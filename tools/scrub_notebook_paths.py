#!/usr/bin/env python3
"""Keep personal filesystem paths out of committed notebooks.

The repository is public. Notebook *outputs* are stored inside the .ipynb, so a
cell that prints a path writes that path into the file and the next commit
publishes it. Scrubbing the sources once does not hold: the next re-run puts the
paths back.

Two different problems, handled two different ways:

* **Outputs** are rewritten in place. They are display text, so substituting a
  placeholder changes nothing that runs, and the reader still sees the shape of
  the path.
* **Sources** are only reported, never rewritten. Editing code automatically is
  how a path constant once became the literal string '$PEPBIND3D_DATA/...',
  which no shell expands and Python happily treats as a directory name. A
  source hit fails the commit and the fix is to import from `paths.py`.

Usage:
    scrub_notebook_paths.py --check  nb.ipynb [...]   # report only, exit 1 on a hit
    scrub_notebook_paths.py          nb.ipynb [...]   # scrub outputs, exit 1 on a source hit
    scrub_notebook_paths.py --self-test
"""
import argparse
import json
import re
import sys
from pathlib import Path

# Longest first: DATA_ROOT sits under HOME, so HOME must not match it first.
# Values come from paths.py when it is importable, so this file holds no
# second copy of the roots to drift out of step.
def _substitutions():
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from paths import CLUSTER_ROOT, DATA_ROOT, MHC_DB_ROOT
        subs = [(str(DATA_ROOT), "<DATA_ROOT>"),
                (str(MHC_DB_ROOT), "<MHC_DB_ROOT>"),
                (str(CLUSTER_ROOT), "<CLUSTER_ROOT>")]
    except Exception:
        subs = []
    # Whatever the roots are, the home directory itself must never ship. Both
    # spellings resolve to the same place on this cluster.
    subs += [("/dors/meilerlab/home/" + _user(), "<HOME>"),
             ("/home/" + _user(), "<HOME>")]
    return sorted(subs, key=lambda kv: -len(kv[0]))


def _user():
    """The account whose paths must not be published."""
    return Path.home().name


def scrub_text(text, subs):
    for real, placeholder in subs:
        text = text.replace(real, placeholder)
    return text


def _walk_outputs(nb):
    for cell in nb.get("cells", []):
        for out in cell.get("outputs", []):
            yield out


def process(path, subs, check_only):
    """Return (n_output_hits, [source hit descriptions]). Writes unless check_only."""
    nb = json.loads(Path(path).read_text())
    pattern = re.compile("|".join(re.escape(r) for r, _ in subs)) if subs else None

    source_hits = []
    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        src = "".join(cell.get("source", []))
        if pattern and pattern.search(src):
            hit = pattern.search(src).group(0)
            source_hits.append(f"{path}: code cell {i} hardcodes {hit}")

    before = json.dumps(nb.get("cells", []))
    for out in _walk_outputs(nb):
        for key in ("text", "traceback"):
            if key in out:
                v = out[key]
                out[key] = ([scrub_text(s, subs) for s in v]
                            if isinstance(v, list) else scrub_text(v, subs))
        data = out.get("data", {})
        for mime, v in list(data.items()):
            if mime.startswith("image/"):          # base64 payloads: leave alone
                continue
            data[mime] = ([scrub_text(s, subs) for s in v]
                          if isinstance(v, list) else
                          scrub_text(v, subs) if isinstance(v, str) else v)
    after = json.dumps(nb.get("cells", []))

    n_out = 0 if before == after else 1
    if n_out and not check_only:
        Path(path).write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    return n_out, source_hits


def self_test():
    """Synthetic notebook with a known answer: one dirty output, one dirty
    source, one image that must survive byte-for-byte."""
    import tempfile

    home = str(Path.home())
    data = home + "/main_project/data/IEDB_data_clean"
    nb = {"cells": [
        {"cell_type": "code", "source": ["print(OUT)\n"],
         "outputs": [{"output_type": "stream", "name": "stdout",
                      "text": ["wrote " + data + "/metadata.csv\n"]}]},
        {"cell_type": "code", "source": ["P = '" + data + "'\n"], "outputs": []},
        {"cell_type": "code", "source": ["plot()\n"],
         "outputs": [{"output_type": "display_data",
                      "data": {"image/png": "iVBORw0KGgoAAAANS"}}]},
        {"cell_type": "markdown", "source": ["prose\n"]},
    ], "metadata": {}, "nbformat": 4, "nbformat_minor": 5}

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "t.ipynb"
        p.write_text(json.dumps(nb))
        subs = _substitutions()

        n_out, hits = process(p, subs, check_only=True)
        assert n_out == 1, f"check mode should see the dirty output, got {n_out}"
        assert len(hits) == 1, f"expected 1 source hit, got {hits}"
        assert json.loads(p.read_text()) == nb, "--check must not write"

        n_out, hits = process(p, subs, check_only=False)
        got = json.loads(p.read_text())
        text = "".join(got["cells"][0]["outputs"][0]["text"])
        assert home not in text, f"home survived in output: {text!r}"
        assert "<" in text, f"no placeholder substituted: {text!r}"
        assert got["cells"][2]["outputs"][0]["data"]["image/png"] == "iVBORw0KGgoAAAANS", \
            "image payload was altered"
        assert data in "".join(got["cells"][1]["source"]), \
            "source was rewritten; it must only be reported"
        assert len(hits) == 1, "source hit must still be reported after scrubbing"

        # Idempotent: a second pass finds nothing to do.
        n_out2, _ = process(p, subs, check_only=False)
        assert n_out2 == 0, "second pass should be a no-op"

    print("self-test passed")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("notebooks", nargs="*", type=Path)
    ap.add_argument("--check", action="store_true",
                    help="report only, do not rewrite")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if not args.notebooks:
        ap.error("give at least one notebook, or --self-test")

    subs = _substitutions()
    scrubbed, source_hits = [], []
    for nb in args.notebooks:
        n_out, hits = process(nb, subs, args.check)
        if n_out:
            scrubbed.append(str(nb))
        source_hits += hits

    for s in scrubbed:
        print(f"{'would scrub' if args.check else 'scrubbed'} paths from outputs: {s}")
    for h in source_hits:
        print(f"ERROR {h}", file=sys.stderr)
    if source_hits:
        print("\nImport the root from paths.py instead of writing the path out.",
              file=sys.stderr)
        return 1
    return 1 if (args.check and scrubbed) else 0


if __name__ == "__main__":
    sys.exit(main())
