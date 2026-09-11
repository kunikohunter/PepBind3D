"""Sensitivity of the score-affinity correlation to the censoring rule.

Runs several censoring policies side by side so the choice can be made on
evidence rather than inherited constants. Reports pooled Spearman rho, n,
and the per-allele summary for each policy.

Usage:
    python censoring_sensitivity.py <HF_DIR> <SCORES_DIR> [OUT_DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

METRIC       = "I_sc_best"
MIN_N_ALLELE = 10

# Each policy: name -> {assay: (exact_values_to_drop, floor_or_None)}
# 'floor' drops everything >= that value, in addition to the exact matches.
POLICIES = {
    "none (keep everything)": {
        "IC50": ((), None),
        "Kd":   ((), None),
    },
    "two spikes only": {
        "IC50": ((20000.0, 70000.0), None),
        "Kd":   ((5000.0, 20000.0), None),
    },
    "manuscript (3 spikes, exact)": {
        "IC50": ((20000.0, 50000.0, 70000.0), None),
        "Kd":   ((5000.0, 10000.0, 20000.0), None),
    },
    "two spikes + tail cutoff": {
        "IC50": ((20000.0, 70000.0), 70000.0),
        "Kd":   ((5000.0, 20000.0), 20000.0),
    },
    "union (3 spikes + tail)": {
        "IC50": ((20000.0, 50000.0, 70000.0, 100000.0), 70000.0),
        "Kd":   ((5000.0, 10000.0, 20000.0), 20000.0),
    },
}


def allele_to_dir(a: str) -> str:
    s = a[4:] if a.startswith("HLA-") else a
    return s.replace("*", "").replace(":", "")


def censored_mask(v: pd.Series, exact, floor) -> pd.Series:
    m = v.isin(exact) if exact else pd.Series(False, index=v.index)
    if floor is not None:
        m = m | (v >= floor)
    return m


def main() -> None:
    hf_dir, scores_dir = Path(sys.argv[1]), Path(sys.argv[2])
    out_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else scores_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(hf_dir / "metadata.csv", low_memory=False)
    scores = pd.read_csv(scores_dir / "score_summary.csv")
    if "flagged" in df.columns:
        df = df[~df["flagged"].astype(bool)].copy()
    df["allele_dir"] = df["allele"].map(allele_to_dir)
    df = df.merge(scores, on=["allele_dir", "peptide"], how="left",
                  validate="many_to_one")
    df = df.dropna(subset=["measurement_value", METRIC]).copy()
    df["log_value"] = np.log10(df["measurement_value"])

    rows = []
    for policy, spec in POLICIES.items():
        for assay in ("IC50", "Kd"):
            exact, floor = spec[assay]
            sub = df[df["measurement_type"] == assay]
            keep = sub[~censored_mask(sub["measurement_value"], exact, floor)]
            if len(keep) < 3:
                continue

            r = stats.spearmanr(keep[METRIC], keep["log_value"])

            # per-allele
            rhos, n_sig = [], 0
            for _, g in keep.groupby("allele"):
                if len(g) < MIN_N_ALLELE:
                    continue
                ra = stats.spearmanr(g[METRIC], g["log_value"])
                rhos.append(ra.statistic)
                n_sig += int(ra.pvalue < 0.05)

            rows.append({
                "policy":     policy,
                "assay":      assay,
                "n":          len(keep),
                "pct_kept":   round(100 * len(keep) / len(sub), 1),
                "rho":        round(r.statistic, 3),
                "p":          f"{r.pvalue:.1e}",
                "n_alleles":  len(rhos),
                "median_rho": round(float(np.median(rhos)), 3) if rhos else np.nan,
                "n_sig":      n_sig,
            })

    res = pd.DataFrame(rows)
    res.to_csv(out_dir / "censoring_sensitivity.csv", index=False)

    for assay in ("IC50", "Kd"):
        print("=" * 92)
        print(f"{assay}   (metric = {METRIC})")
        print("=" * 92)
        print(res[res["assay"] == assay]
              .drop(columns="assay")
              .to_string(index=False))
        print()

    print(f"Wrote {out_dir / 'censoring_sensitivity.csv'}")
    print("\nInterpretation: if rho is stable across policies, the censoring")
    print("choice is not driving the result and any defensible rule will do.")


if __name__ == "__main__":
    main()
