"""Rebuild the released metadata.csv with the full FlexPepDock score set.

The originally released metadata.csv carries only total_score summaries
(rosetta_best_score, rosetta_mean_score). This adds the interface (I_sc)
and reweighted (reweighted_sc) score summaries parsed from the score.sc
files, so the public file reproduces the analyses in the manuscript.

Adds columns: I_sc_best, I_sc_mean, reweighted_sc_best, reweighted_sc_mean,
              total_score_best, total_score_mean
Keeps rosetta_best_score / rosetta_mean_score as aliases for back-compat.

Usage:
    python rebuild_metadata.py <old_metadata.csv> <score_summary.csv> <out_metadata.csv>
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd


def allele_to_dir(a: str) -> str:
    s = a[4:] if a.startswith("HLA-") else a
    return s.replace("*", "").replace(":", "")


def main() -> None:
    old_fn, scores_fn, out_fn = map(Path, sys.argv[1:4])

    df = pd.read_csv(old_fn, low_memory=False)
    scores = pd.read_csv(scores_fn)
    print(f"metadata rows: {len(df):,}")
    print(f"score_summary pairs: {len(scores):,}")
    print(f"score_summary columns: {list(scores.columns)}")

    df["allele_dir"] = df["allele"].map(allele_to_dir)

    score_cols = [c for c in scores.columns if c not in ("allele_dir", "peptide")]
    # drop any of these if a previous run added them, to stay idempotent
    df = df.drop(columns=[c for c in score_cols if c in df.columns], errors="ignore")

    before = len(df)
    df = df.merge(scores[["allele_dir", "peptide"] + score_cols],
                  on=["allele_dir", "peptide"], how="left", validate="many_to_one")
    assert len(df) == before, "merge changed row count"

    # Report coverage
    for c in ("I_sc_best", "reweighted_sc_best", "total_score_best"):
        if c in df.columns:
            n_missing = df[c].isna().sum()
            print(f"  {c}: {n_missing:,} rows without a value "
                  f"({n_missing/len(df)*100:.2f}%)")

    # Keep legacy aliases pointing at total_score (what they always were)
    if "total_score_best" in df.columns:
        df["rosetta_best_score"] = df["total_score_best"]
        df["rosetta_mean_score"] = df["total_score_mean"]

    df = df.drop(columns=["allele_dir"])
    df.to_csv(out_fn, index=False)
    print(f"\nWrote {out_fn} ({len(df):,} rows, {len(df.columns)} columns)")
    print("Columns:", list(df.columns))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main()
