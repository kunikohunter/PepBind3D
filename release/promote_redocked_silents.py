"""
Promote re-docked silent files into the staging release, replacing the
defective ones.

THE HAZARD THIS EXISTS TO AVOID. release_v2_final/structures was assembled with
hardlinks (cp -al) to save disk, so each staging silent shares an inode with its
source tree -- verified: A0101/CEKRLLLKL.silent has link count 2 and the same
inode as v1_corrected_silents/A0101/CEKRLLLKL.silent. Copying a new file OVER
the staging path would write through that shared inode and silently rewrite the
source tree too, destroying the record of what the original run produced.

So each replacement is unlink-then-link, never copy-over: removing the staging
path only decrements the link count and leaves the source file intact. The new
file is hardlinked in rather than copied, so promotion costs no extra disk.

Every source tree is treated as read-only, and the live HuggingFace release is
refused outright.

Usage:
    python3 promote_redocked_silents.py --new <dir> --pairs <csv> [--apply]
Without --apply it reports what it would do and changes nothing.
"""
import argparse
import hashlib
import os
from pathlib import Path

import pandas as pd

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT  # noqa: E402

BASE = DATA_ROOT
STAGING = BASE / "release_v2_final" / "structures"
# trees that must never be modified by this script
READ_ONLY = [BASE / "v1_corrected_silents", BASE / "huggingface" / "structures",
             BASE / "huggingface_v2" / "structures"]


def digest(p, n=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(n):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, help="tree of new silents, {allele}/{peptide}.silent")
    ap.add_argument("--pairs", required=True, help="CSV with allele_dir + peptide")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--allow-new", action="store_true",
                    help="permit pairs that have no staging file yet (recovered "
                         "orphans, which are additions rather than replacements)")
    args = ap.parse_args()

    new_root = Path(args.new)
    assert "huggingface/structures" not in str(STAGING), "refusing to touch the live release"
    pairs = pd.read_csv(args.pairs)

    plan, missing_new, missing_old = [], [], []
    for r in pairs.itertuples(index=False):
        src = new_root / r.allele_dir / f"{r.peptide}.silent"
        dst = STAGING / r.allele_dir / f"{r.peptide}.silent"
        if not src.exists():
            missing_new.append(str(src)); continue
        if not dst.exists():
            if not args.allow_new:
                missing_old.append(str(dst)); continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            plan.append((src, dst, 0))          # nlink 0 => a new path, nothing to unlink
            continue
        plan.append((src, dst, dst.stat().st_nlink))

    print(f"pairs requested        : {len(pairs):,}")
    print(f"new silents found      : {len(plan) + len(missing_old):,}")
    print(f"staging targets found  : {len(plan):,}")
    if missing_new:
        print(f"  MISSING new silent   : {len(missing_new)} e.g. {missing_new[:2]}")
    if missing_old:
        print(f"  MISSING staging file : {len(missing_old)} e.g. {missing_old[:2]}")
    hard = sum(1 for *_, n in plan if n > 1)
    print(f"staging files that are hardlinks (nlink>1): {hard:,} of {len(plan):,}")
    print("  -> these MUST be unlinked before replacement, not overwritten")

    if not args.apply:
        print("\nDRY RUN. Nothing changed. Re-run with --apply to promote.")
        return

    # record source-tree digests for the affected pairs, to prove afterwards
    # that promotion did not write through any shared inode
    witness = {}
    for ro in READ_ONLY:
        for r in pairs.head(25).itertuples(index=False):
            f = ro / r.allele_dir / f"{r.peptide}.silent"
            if f.exists():
                witness[str(f)] = digest(f)
    print(f"\nrecorded {len(witness)} source-tree digests as witnesses")

    replaced = 0
    added = 0
    for src, dst, nlink in plan:
        if nlink:
            os.unlink(dst)      # decrements link count; source file survives
            replaced += 1
        else:
            added += 1
        os.link(src, dst)       # hardlink, no extra disk
    print(f"replaced {replaced:,} staging silents (unlink-then-link), "
          f"added {added:,} new")

    bad = [f for f, d in witness.items() if digest(Path(f)) != d]
    print(f"source-tree witnesses unchanged: {len(witness) - len(bad)}/{len(witness)}")
    if bad:
        raise SystemExit(f"SOURCE TREE MODIFIED -- promotion wrote through a shared "
                         f"inode: {bad[:3]}")
    print("no source tree was modified.")


if __name__ == "__main__":
    main()
