"""Compare score metrics against experimental affinity before committing to one.

Merges score_summary.csv onto metadata.csv and reports pooled and per-allele
Spearman correlations for every metric x {best, mean} x {IC50, Kd}.

Usage:
    python quick_metric_comparison.py <HF_DIR> <SCORES_DIR> [OUT_DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

IC50_CEILINGS = (20000.0, 50000.0, 70000.0)
KD_CEILINGS   = (5000.0, 10000.0, 20000.0)
MIN_N_ALLELE  = 10


def allele_to_dir(a: str) -> str:
    s = a[4:] if a.startswith("HLA-") else a
    return s.replace("*", "").replace(":", "")


def is_censored(v: pd.Series, ceilings) -> pd.Series:
    return v.isin(ceilings) | (v >= max(ceilings))


def main() -> None:
    hf_dir     = Path(sys.argv[1])
    scores_dir = Path(sys.argv[2])
    out_dir    = Path(sys.argv[3]) if len(sys.argv) > 3 else scores_dir

    df = pd.read_csv(hf_dir / "metadata.csv", low_memory=False)
    scores = pd.read_csv(scores_dir / "score_summary.csv")
    print(f"metadata rows      : {len(df):,}")
    print(f"score_summary pairs: {len(scores):,}")

    if "flagged" in df.columns:
        df = df[~df["flagged"].astype(bool)].copy()
        print(f"after dropping flagged: {len(df):,}")

    df["allele_dir"] = df["allele"].map(allele_to_dir)
    before = len(df)
    df = df.merge(scores, on=["allele_dir", "peptide"], how="left",
                  validate="many_to_one")
    assert len(df) == before, "merge changed row count"

    metrics = sorted({c.rsplit("_", 1)[0] for c in scores.columns
                      if c.endswith(("_best", "_mean"))})
    print(f"metrics found      : {metrics}")

    unmatched = df[f"{metrics[0]}_best"].isna().sum()
    print(f"rows without score : {unmatched:,} ({unmatched/len(df)*100:.2f}%)\n")

    rows, per_allele_rows = [], []
    for assay, ceil in (("IC50", IC50_CEILINGS), ("Kd", KD_CEILINGS)):
        sub = df[df["measurement_type"] == assay].copy()
        sub = sub.dropna(subset=["measurement_value"])
        sub = sub[~is_censored(sub["measurement_value"], ceil)]
        sub["log_value"] = np.log10(sub["measurement_value"])

        for m in metrics:
            for agg in ("best", "mean"):
                col = f"{m}_{agg}"
                d = sub.dropna(subset=[col])
                if len(d) < 3:
                    continue
                r = stats.spearmanr(d[col], d["log_value"])
                rows.append({"assay": assay, "metric": m, "agg": agg,
                             "rho": r.statistic, "p": r.pvalue, "n": len(d)})

                for allele, g in d.groupby("allele"):
                    if len(g) < MIN_N_ALLELE:
                        continue
                    ra = stats.spearmanr(g[col], g["log_value"])
                    per_allele_rows.append({
                        "assay": assay, "metric": m, "agg": agg,
                        "allele": allele, "rho": ra.statistic,
                        "p": ra.pvalue, "n": len(g)})

    pooled = pd.DataFrame(rows)
    per_allele = pd.DataFrame(per_allele_rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    pooled.to_csv(out_dir / "metric_comparison_pooled.csv", index=False)
    per_allele.to_csv(out_dir / "metric_comparison_per_allele.csv", index=False)

    print("=" * 72)
    print("POOLED (lower score should track lower affinity value -> positive rho)")
    print("=" * 72)
    print(pooled.assign(rho=pooled["rho"].round(3),
                        p=pooled["p"].map(lambda x: f"{x:.2g}"))
                .to_string(index=False))

    print()
    print("=" * 72)
    print(f"PER-ALLELE SUMMARY (alleles with n >= {MIN_N_ALLELE})")
    print("=" * 72)
    summ = (per_allele.groupby(["assay", "metric", "agg"])
            .apply(lambda g: pd.Series({
                "n_alleles": len(g),
                "median_rho": round(g["rho"].median(), 3),
                "n_sig": int((g["p"] < 0.05).sum()),
                "n_sig_pos": int(((g["p"] < 0.05) & (g["rho"] > 0)).sum()),
                "n_sig_neg": int(((g["p"] < 0.05) & (g["rho"] < 0)).sum()),
            }), include_groups=False)
            .reset_index())
    print(summ.to_string(index=False))

    print()
    print("Best metric by pooled |rho| per assay:")
    for assay, g in pooled.groupby("assay"):
        b = g.loc[g["rho"].abs().idxmax()]
        print(f"  {assay}: {b['metric']}_{b['agg']}  rho={b['rho']:.3f}  "
              f"p={b['p']:.2g}  n={b['n']:,}")
    print("\nBest metric by median per-allele rho:")
    for assay, g in summ.groupby("assay"):
        b = g.loc[g["median_rho"].idxmax()]
        print(f"  {assay}: {b['metric']}_{b['agg']}  median rho={b['median_rho']}  "
              f"({b['n_sig']}/{b['n_alleles']} significant)")


if __name__ == "__main__":
    main()
