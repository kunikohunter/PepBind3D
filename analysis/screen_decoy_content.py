"""
Screen released decoys for CONTENT defects, not just decoy count.

Three defects reach the release because nothing validated what is inside a
decoy, only how many decoys there are:

  A. spurious extra chain, remnant too small to touch the receptor
     -> 3 chains, ALL interface terms exactly 0.000 (score-visible)
  B. spurious extra chain, remnant large enough to form an interface
     -> 3 chains, plausible-looking I_sc computed for the WRONG chain
        (INVISIBLE to any score-based screen)
  C. truncated peptide
     -> 2 chains, but chain B holds only part of the curated sequence
        (INVISIBLE to any score-based screen; the structure does not contain
        the peptide it claims)

All three come from threading a shorter query peptide onto a longer template:
surplus template residues are left behind as their own chain, and FlexPepDock
then refines whichever chain it takes to be the peptide. See
the threading step
for that run.

The reliable test is chain composition plus sequence identity against the
curated peptide -- exactly what ACCRE's stage_redock.py / stage_regen.py assert
before docking, and what release/convert_to_silent.py now enforces before
writing a silent file. This script applies it retrospectively to what is
already on disk.

Reads ONE decoy per pair (the artifact is created at threading, so it is
present in every decoy of an affected pair) and only its CA records.

Usage:
    python3 screen_decoy_content.py --out-dir <dir> [--workers 4] [--limit N]
    python3 screen_decoy_content.py --self-test

NB only pairs whose decoy PDB tree is present locally can be screened. On
Tungsten that is the v1 batch; the v2 decoy PDBs were deleted after the silents
were built, so v2 must be screened where its tree lives.
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

# one implementation, shared with the release gate in
# release/convert_to_silent.py -- see release/decoy_content.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "release"))
from decoy_content import chain_sequences, classify, self_test  # noqa: E402

BASE = Path("<HOME>/main_project/data/IEDB_data_clean")
PDB_ROOT = BASE / "pdb"
SCORE_SUMMARIES = [
    BASE / "IEDB_validation" / "scores_out" / "score_summary.csv",
    BASE / "IEDB_validation" / "scores_out_v2" / "score_summary.csv",
]

def _one(args):
    allele_dir, peptide = args
    d = PDB_ROOT / allele_dir / peptide
    pdbs = sorted(d.glob(f"{peptide}_input_[0-9][0-9][0-9][0-9].pdb"))
    if not pdbs:
        return {"allele_dir": allele_dir, "peptide": peptide,
                "verdict": "no_pdbs", "detail": str(d)}
    try:
        seqs = chain_sequences(pdbs[0])
    except OSError as e:
        return {"allele_dir": allele_dir, "peptide": peptide,
                "verdict": "unreadable", "detail": str(e)}
    verdict, detail = classify(seqs, peptide)
    return {"allele_dir": allele_dir, "peptide": peptide, "decoy": pdbs[0].name,
            "verdict": verdict, "detail": detail}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir")
    ap.add_argument("--workers", type=int, default=4,
                    help="parallel readers. Keep modest: this is a shared "
                         "filesystem and an aggressive walk degrades the machine.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test(); return
    if not args.out_dir:
        raise SystemExit("--out-dir is required")
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    self_test(); print()

    s = pd.concat([pd.read_csv(p) for p in SCORE_SUMMARIES], ignore_index=True)
    jobs = list(zip(s["allele_dir"], s["peptide"]))
    if args.limit:
        jobs = jobs[:args.limit]
    print(f"screening {len(jobs):,} pairs with {args.workers} workers", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, r in enumerate(ex.map(_one, jobs, chunksize=64), 1):
            rows.append(r)
            if i % 5000 == 0:
                bad = sum(1 for x in rows if x["verdict"] not in ("ok", "no_pdbs"))
                print(f"  {i:,}/{len(jobs):,}  defects so far: {bad}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out / "decoy_content_screen.csv", index=False)
    bad = df[~df["verdict"].isin(["ok", "no_pdbs"])]
    bad.to_csv(out / "decoy_content_defects.csv", index=False)

    print("\nverdicts:")
    for v, n in df["verdict"].value_counts().items():
        print(f"  {v:20s} {n:,}")
    screened = int((df["verdict"] != "no_pdbs").sum())
    if screened:
        print(f"\ndefective: {len(bad):,} of {screened:,} screened "
              f"({100*len(bad)/screened:.3f}%)")
    print(f"\nwrote decoy_content_screen.csv and decoy_content_defects.csv to {out}")


if __name__ == "__main__":
    main()
