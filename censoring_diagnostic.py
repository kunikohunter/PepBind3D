"""Diagnose the censoring filter: what values pile up, and how much
the correlation depends on excluding them.

Usage:
    python censoring_diagnostic.py <HF_DIR> <SCORES_DIR>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

IC50_CEILINGS = (20000.0, 70000.0)
KD_CEILINGS   = (5000.0, 20000.0)
METRIC = "I_sc_best"


def allele_to_dir(a: str) -> str:
    s = a[4:] if a.startswith("HLA-") else a
    return s.replace("*", "").replace(":", "")


def main() -> None:
    hf_dir, scores_dir = Path(sys.argv[1]), Path(sys.argv[2])

    df = pd.read_csv(hf_dir / "metadata.csv", low_memory=False)
    scores = pd.read_csv(scores_dir / "score_summary.csv")
    if "flagged" in df.columns:
        df = df[~df["flagged"].astype(bool)].copy()
    df["allele_dir"] = df["allele"].map(allele_to_dir)
    df = df.merge(scores, on=["allele_dir", "peptide"], how="left",
                  validate="many_to_one")

    for assay, ceilings in (("IC50", IC50_CEILINGS), ("Kd", KD_CEILINGS)):
        sub = df[df["measurement_type"] == assay]
        v = sub["measurement_value"]
        print("=" * 68)
        print(f"{assay}")
        print("=" * 68)
        print(f"  rows                       : {len(sub):,}")
        print(f"  with a value               : {v.notna().sum():,}")
        print(f"  with value AND {METRIC:<12}: "
              f"{(v.notna() & sub[METRIC].notna()).sum():,}")

        print(f"\n  20 most common values (a spike = a detection limit):")
        top = v.value_counts().head(20)
        for val, n in top.items():
            flag = "  <-- currently treated as censored" if (
                val in ceilings or val >= max(ceilings)) else ""
            print(f"    {val:>12,.1f}  n={n:>6,}  ({n/len(sub)*100:>5.2f}%){flag}")

        cens = v.isin(ceilings) | (v >= max(ceilings))
        print(f"\n  currently excluded as censored: {cens.sum():,} "
              f"({cens.sum()/len(sub)*100:.1f}%)")
        print(f"  at/above max ceiling only     : {(v >= max(ceilings)).sum():,}")
        print(f"  exactly equal to a ceiling    : {v.isin(ceilings).sum():,}")

        # Sensitivity: correlation with vs without censored values
        d = sub.dropna(subset=["measurement_value", METRIC]).copy()
        d["log_value"] = np.log10(d["measurement_value"])
        d_cens = d["measurement_value"].isin(ceilings) | \
                 (d["measurement_value"] >= max(ceilings))

        print(f"\n  Spearman {METRIC} vs log10({assay}):")
        for label, dd in (("excluding censored (current)", d[~d_cens]),
                          ("including everything        ", d)):
            if len(dd) > 3:
                r = stats.spearmanr(dd[METRIC], dd["log_value"])
                print(f"    {label}: rho={r.statistic:+.3f}  "
                      f"p={r.pvalue:.2g}  n={len(dd):,}")

        # How many distinct values sit at the very top of the range?
        hi = v[v >= v.quantile(0.98)]
        print(f"\n  top 2% of values span {hi.min():,.0f} to {hi.max():,.0f} nM "
              f"across {hi.nunique():,} distinct values")
        print()


if __name__ == "__main__":
    main()
