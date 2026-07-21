"""Extract Rosetta score terms from per-pair FlexPepDock score.sc files.

The released silent files carry only the ref2015 energy terms, but the
local pdb/ tree retains the full FlexPepDock score set written by
FlexPepDocking (-pep_refine), including reweighted_sc and I_sc.

Expected layout:
    {PDB_ROOT}/{allele_dir}/{peptide}/score.sc

Outputs (written to OUT_DIR):
    per_decoy_scores.csv  one row per decoy, all score terms
    score_summary.csv     one row per peptide-allele pair, best + mean

Usage:
    python parse_scorefiles.py <PDB_ROOT> <OUT_DIR> [--metrics a,b,c]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_METRICS = ("reweighted_sc", "I_sc", "total_score")


def parse_scorefile(path: Path) -> pd.DataFrame:
    """Parse one Rosetta score.sc into a DataFrame.

    Line 1 is 'SEQUENCE:'; the header is the first SCORE: line whose second
    token is non-numeric. Later SCORE: lines are data. Handles repeated
    headers (from -overwrite / appended runs) and short/malformed rows.
    """
    header, rows = None, []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith("SCORE:"):
                continue
            fields = line.split()[1:]
            if not fields:
                continue
            if header is None:
                if not _is_number(fields[0]):
                    header = fields
                continue
            if fields[0] == header[0]:      # repeated header block
                continue
            if len(fields) != len(header):  # truncated / interrupted write
                continue
            rows.append(fields)

    if header is None or not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows, columns=header)
    for c in df.columns:
        if c != "description":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb_root", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--metrics", default=",".join(DEFAULT_METRICS),
                    help="comma-separated score terms to summarize")
    ap.add_argument("--expect-decoys", type=int, default=25)
    args = ap.parse_args()

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    scorefiles = sorted(args.pdb_root.glob("*/*/score.sc"))
    print(f"Found {len(scorefiles):,} score.sc files under {args.pdb_root}")
    if not scorefiles:
        raise SystemExit("No score.sc files matched {allele_dir}/{peptide}/score.sc")

    frames, empty = [], []
    for i, sf in enumerate(scorefiles, 1):
        df = parse_scorefile(sf)
        if df.empty:
            empty.append(str(sf))
            continue
        df["allele_dir"] = sf.parent.parent.name
        df["peptide"] = sf.parent.name
        frames.append(df)
        if i % 2000 == 0:
            print(f"  {i:,}/{len(scorefiles):,}")

    if not frames:
        raise SystemExit("No SCORE: data rows parsed.")

    decoys = pd.concat(frames, ignore_index=True)
    print(f"\nParsed {len(decoys):,} decoys from {len(frames):,} files")
    if empty:
        print(f"WARNING: {len(empty)} score.sc files had no data rows; first few:")
        for e in empty[:5]:
            print("   ", e)

    missing = [m for m in metrics if m not in decoys.columns]
    if missing:
        print(f"\nWARNING: requested metrics absent: {missing}")
        metrics = [m for m in metrics if m in decoys.columns]
    if not metrics:
        raise SystemExit("None of the requested metrics are present.")
    print(f"Summarizing: {metrics}")

    decoys.to_csv(args.out_dir / "per_decoy_scores.csv", index=False)

    g = decoys.groupby(["allele_dir", "peptide"], sort=True)
    summary = pd.DataFrame({"n_decoys": g.size()})
    for m in metrics:
        summary[f"{m}_best"] = g[m].min()    # Rosetta: lower is better
        summary[f"{m}_mean"] = g[m].mean()
    summary = summary.reset_index()

    out = args.out_dir / "score_summary.csv"
    summary.to_csv(out, index=False)
    print(f"\nWrote {len(summary):,} pair summaries -> {out}")
    print(summary.head().to_string(index=False))

    bad = summary[summary["n_decoys"] != args.expect_decoys]
    if len(bad):
        print(f"\nWARNING: {len(bad):,} pairs != {args.expect_decoys} decoys")
        print(bad.head(10).to_string(index=False))
    else:
        print(f"\nAll pairs have exactly {args.expect_decoys} decoys.")

    print("\nMetric correlations across decoys (sanity check):")
    print(decoys[metrics].corr(method="spearman").round(3).to_string())


if __name__ == "__main__":
    main()
