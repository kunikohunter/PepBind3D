"""
Recompute the censored-vs-quantitative separation and the score-affinity
Spearman correlations on the FULL merged v2 dataset (112,378 pairs), instead of
v1's 49,268. Merges v1 + v2 per-pair score summaries, joins to the reconciled
merged metadata, and reuses the validated analysis functions from
censored_vs_quantitative_auroc.py (same censoring definition, same metric).

Usage: python3 recompute_affinity_112k.py --out-dir <dir>
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from censored_vs_quantitative_auroc import (
    is_censored, auroc_and_effect_size,
    IC50_CEILINGS, KD_CEILINGS, IC50_FLOOR, KD_FLOOR, STRONG_BINDER_THRESHOLD_NM,
    METRICS, PRIMARY_METRIC,
)

METADATA = Path("<HOME>/main_project/data/IEDB_data_clean/release_v2_final/metadata.csv")
V1_SCORES = Path("<HOME>/main_project/data/IEDB_data_clean/IEDB_validation/scores_out/score_summary.csv")
V2_SCORES = Path("<HOME>/main_project/data/IEDB_data_clean/IEDB_validation/scores_out_v2/score_summary.csv")

SCORE_COLS = [f"{m}_{agg}" for m in METRICS for agg in ("best", "mean")]


def build_scored_metadata():
    md = pd.read_csv(METADATA, low_memory=False)
    s1 = pd.read_csv(V1_SCORES); s2 = pd.read_csv(V2_SCORES)
    scores = pd.concat([s1, s2], ignore_index=True).drop_duplicates(subset=["allele_dir", "peptide"])
    keep = ["allele_dir", "peptide"] + [c for c in SCORE_COLS if c in scores.columns]
    merged = md.merge(scores[keep], left_on=["allele_compact", "peptide"],
                      right_on=["allele_dir", "peptide"], how="left")
    n_missing = merged[f"{PRIMARY_METRIC}_best"].isna().sum()
    print(f"merged metadata rows: {len(merged):,}; pairs missing a score: "
          f"{n_missing:,} ({100*n_missing/len(merged):.2f}%)", flush=True)
    return merged


def run_censored(df, out_dir):
    rows = []
    for assay, ceil, floor in [("IC50", IC50_CEILINGS, IC50_FLOOR), ("Kd", KD_CEILINGS, KD_FLOOR)]:
        sub = df[df["measurement_type"] == assay]
        for metric in METRICS:
            for agg in ("best", "mean"):
                col = f"{metric}_{agg}"
                if col not in sub: continue
                d = sub.dropna(subset=[col, "measurement_value"])
                cens = is_censored(d["measurement_value"], ceil, floor)
                r = auroc_and_effect_size(d[col], ~cens, n_boot=500)
                if r: rows.append({"assay": assay, "comparison": "pooled", "metric": col, **r})
                strong = (~cens) & (d["measurement_value"] < STRONG_BINDER_THRESHOLD_NM)
                dd = d[strong | cens]; cens2 = is_censored(dd["measurement_value"], ceil, floor)
                r2 = auroc_and_effect_size(dd[col], ~cens2, n_boot=500)
                if r2: rows.append({"assay": assay, "comparison": "strong_binder_only", "metric": col, **r2})
    pd.DataFrame(rows).to_csv(out_dir / "censored_vs_quantitative_112k.csv", index=False)
    return rows


def run_affinity(df, out_dir):
    out = {}
    for assay, ceil, floor in [("IC50", IC50_CEILINGS, IC50_FLOOR), ("Kd", KD_CEILINGS, KD_FLOOR)]:
        sub = df[df["measurement_type"] == assay].dropna(subset=[f"{PRIMARY_METRIC}_best", "measurement_value"])
        cens = is_censored(sub["measurement_value"], ceil, floor)
        q = sub[~cens].copy()
        q = q[q["measurement_value"] > 0]
        logv = np.log10(q["measurement_value"])
        out[assay] = {}
        for metric in METRICS:
            col = f"{metric}_best"
            if col not in q: continue
            rho = spearmanr(q[col], logv)
            out[assay][col] = {"rho": float(rho.statistic), "p": float(rho.pvalue), "n": int(len(q))}
    with open(out_dir / "score_affinity_112k.json", "w") as f:
        json.dump(out, f, indent=2)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(); out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    df = build_scored_metadata()
    cens = run_censored(df, out)
    aff = run_affinity(df, out)
    print("\n=== censored-vs-quantitative (I_sc_best, pooled) ===")
    for r in cens:
        if r["metric"] == "I_sc_best" and r["comparison"] == "pooled":
            print(f"  {r['assay']}: AUROC={r['auroc']:.3f} [{r['auroc_ci_lo']:.3f},{r['auroc_ci_hi']:.3f}] "
                  f"n={r['n']:,} (quant={r['n_quantitative']:,}, cens={r['n_censored']:,})")
    print("\n=== score-affinity (I_sc_best, quantitative) ===")
    for assay, d in aff.items():
        if "I_sc_best" in d:
            print(f"  {assay}: rho={d['I_sc_best']['rho']:.3f} p={d['I_sc_best']['p']:.2g} n={d['I_sc_best']['n']:,}")
    print(f"\nwrote outputs to {out}")


if __name__ == "__main__":
    main()
