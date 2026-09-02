"""
Phase 1 item 1.1 (REVISION_PLAN.md): censored vs quantitative separation.

~20K pairs at assay detection limits were dropped from the Validation 2
correlation analysis (02_score_affinity_validation.ipynb). This script asks
the complementary question: can the bundled Rosetta scores (I_sc,
reweighted_sc, total_score) tell a censored (non-binder / weak-binder-at-best)
measurement apart from a quantitative (real, measured) one? Reported as
AUROC (pooled + per-allele) and a rank-biserial effect size, plus a
strong-binder-only variant (IC50/KD < 500 nM vs censored) per the plan's own
suggestion.

Uses the exact same censoring definition, ceilings, and score-column
convention as 02_score_affinity_validation.ipynb (kept in sync deliberately
so "censored" means the same thing everywhere in the paper).

Usage:
    $ python3 censored_vs_quantitative_auroc.py --out-dir <output dir>

Self-test:
    $ python3 censored_vs_quantitative_auroc.py --self-test
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

HF_DIR = Path("/home/huntek1/main_project/data/IEDB_data_clean/huggingface")

COL_ALLELE = "allele"
COL_MEAS_TYPE = "measurement_type"
COL_MEAS_VALUE = "measurement_value"
COL_FLAGGED = "flagged"

METRICS = ["I_sc", "reweighted_sc", "total_score"]
PRIMARY_METRIC = "I_sc"

# Identical to 02_score_affinity_validation.ipynb, deliberately.
IC50_CEILINGS = (20000.0, 50000.0, 70000.0)
KD_CEILINGS = (5000.0, 10000.0, 20000.0)
IC50_FLOOR = 70000.0
KD_FLOOR = 20000.0
STRONG_BINDER_THRESHOLD_NM = 500.0


def is_censored(values, ceilings, floor):
    return values.isin(ceilings) | (values >= floor)


def auroc_and_effect_size(score, is_quantitative):
    """AUROC for discriminating quantitative (favorable/lower score expected)
    from censored (unfavorable/higher score expected) pairs, plus the
    rank-biserial effect size (r = 2*AUROC - 1, equivalent to Cliff's delta).
    Lower Rosetta REU = more favorable, so we score on -value: higher -value
    (i.e. lower raw score) should predict "quantitative" (a real binder)."""
    y = is_quantitative.astype(int).values
    if y.sum() == 0 or (1 - y).sum() == 0:
        return None
    auc = roc_auc_score(y, -score.values)
    return {"auroc": float(auc), "rank_biserial": float(2 * auc - 1), "n": int(len(y)),
            "n_quantitative": int(y.sum()), "n_censored": int((1 - y).sum())}


def analyze_assay(df, assay_label, meas_type, ceilings, floor, strong_binder_thr):
    sub = df[df[COL_MEAS_TYPE] == meas_type].copy()
    censored = is_censored(sub[COL_MEAS_VALUE], ceilings, floor)

    results = {"assay": assay_label, "pooled": {}, "pooled_strong_binder_only": {}, "per_allele": []}

    for metric in METRICS:
        for agg in ("best", "mean"):
            col = f"{metric}_{agg}"
            d = sub.dropna(subset=[col, COL_MEAS_VALUE])
            cens_d = is_censored(d[COL_MEAS_VALUE], ceilings, floor)
            r = auroc_and_effect_size(d[col], ~cens_d)
            if r is not None:
                results["pooled"][f"{metric}_{agg}"] = r

            # strong-binder-only variant: quantitative group restricted to
            # value < STRONG_BINDER_THRESHOLD_NM, vs the same censored group
            strong = (~cens_d) & (d[COL_MEAS_VALUE] < strong_binder_thr)
            keep = strong | cens_d
            dd = d[keep]
            cens_dd = is_censored(dd[COL_MEAS_VALUE], ceilings, floor)
            r2 = auroc_and_effect_size(dd[col], ~cens_dd)
            if r2 is not None:
                results["pooled_strong_binder_only"][f"{metric}_{agg}"] = r2

    # per-allele, primary metric only (best), n>=10 per group, matching the
    # per-allele Spearman convention already used in the paper. Outside the
    # metric loop deliberately -- only needs to run once, using PRIMARY_METRIC.
    col = f"{PRIMARY_METRIC}_best"
    for allele, g in sub.groupby(COL_ALLELE):
        g2 = g.dropna(subset=[col, COL_MEAS_VALUE])
        cens_g = is_censored(g2[COL_MEAS_VALUE], ceilings, floor)
        if (~cens_g).sum() < 10 or cens_g.sum() < 10:
            continue
        r = auroc_and_effect_size(g2[col], ~cens_g)
        if r is not None:
            results["per_allele"].append({"allele": allele, "metric": f"{PRIMARY_METRIC}_best", **r})

    return results


def self_test():
    """Synthetic data with a known analytic AUROC. For two normals with
    means mu0 (censored, unfavorable/high score) and mu1 (quantitative,
    favorable/low score) and equal variance sigma, the true AUROC for
    discriminating on -score is Phi((mu0 - mu1) / (sigma*sqrt(2)))."""
    rng = np.random.default_rng(0)
    n = 5000
    sigma = 3.0
    mu_censored, mu_quant = 5.0, -2.0  # censored scores are worse (higher)
    censored_scores = rng.normal(mu_censored, sigma, size=n)
    quant_scores = rng.normal(mu_quant, sigma, size=n)

    expected_auc = norm.cdf((mu_censored - mu_quant) / (sigma * np.sqrt(2)))

    df = pd.DataFrame({
        "score": np.concatenate([censored_scores, quant_scores]),
        "is_quantitative": np.array([False] * n + [True] * n),
    })
    r = auroc_and_effect_size(df["score"], df["is_quantitative"])
    assert abs(r["auroc"] - expected_auc) < 0.02, (
        f"self-test FAILED: computed AUROC {r['auroc']:.4f} vs analytic {expected_auc:.4f}"
    )

    # Null case: identical distributions -> AUROC ~ 0.5
    same_a = rng.normal(0, 1, size=n)
    same_b = rng.normal(0, 1, size=n)
    df_null = pd.DataFrame({
        "score": np.concatenate([same_a, same_b]),
        "is_quantitative": np.array([False] * n + [True] * n),
    })
    r_null = auroc_and_effect_size(df_null["score"], df_null["is_quantitative"])
    assert abs(r_null["auroc"] - 0.5) < 0.02, (
        f"self-test FAILED: null-case AUROC {r_null['auroc']:.4f} should be ~0.5"
    )

    print(f"Self-test PASSED: computed AUROC {r['auroc']:.4f} matches analytic {expected_auc:.4f} "
          f"(known-separation case); null case AUROC {r_null['auroc']:.4f} ~ 0.5 as expected.")


def plot_distributions(df, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    for ax, (assay_label, meas_type, ceilings, floor) in zip(
        axes, [("IC50", "IC50", IC50_CEILINGS, IC50_FLOOR), ("KD", "Kd", KD_CEILINGS, KD_FLOOR)]
    ):
        sub = df[df[COL_MEAS_TYPE] == meas_type].dropna(subset=[f"{PRIMARY_METRIC}_best", COL_MEAS_VALUE])
        cens = is_censored(sub[COL_MEAS_VALUE], ceilings, floor)
        data = [sub.loc[~cens, f"{PRIMARY_METRIC}_best"], sub.loc[cens, f"{PRIMARY_METRIC}_best"]]
        bp = ax.boxplot(data, labels=[f"quantitative\n(n={len(data[0]):,})", f"censored\n(n={len(data[1]):,})"],
                         showfliers=False, patch_artist=True)
        for patch, color in zip(bp["boxes"], ["#4477AA", "#CC3311"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax.set_title(assay_label, fontsize=10)
        ax.set_ylabel(f"{PRIMARY_METRIC}_best (REU)")
    fig.suptitle("I_sc_best: censored vs quantitative measurements", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_dir / "censored_vs_quantitative_boxplots.pdf")
    fig.savefig(out_dir / "censored_vs_quantitative_boxplots.png", dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=str, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    if args.out_dir is None:
        raise SystemExit("--out-dir is required (never write outputs into the scripts directory)")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(HF_DIR / "metadata.csv", low_memory=False)
    df = df[~df[COL_FLAGGED].astype(bool)].copy()
    print(f"Loaded {len(df):,} unflagged rows from {HF_DIR / 'metadata.csv'}")

    all_results = {}
    for assay_label, meas_type, ceilings, floor in [
        ("IC50", "IC50", IC50_CEILINGS, IC50_FLOOR),
        ("KD", "Kd", KD_CEILINGS, KD_FLOOR),
    ]:
        res = analyze_assay(df, assay_label, meas_type, ceilings, floor, STRONG_BINDER_THRESHOLD_NM)
        all_results[assay_label] = res

        print(f"\n=== {assay_label} pooled AUROC (censored vs quantitative) ===")
        for metric, r in res["pooled"].items():
            print(f"  {metric:<20} AUROC={r['auroc']:.3f}  rank-biserial={r['rank_biserial']:+.3f}  "
                  f"n={r['n']:,} (quant={r['n_quantitative']:,}, censored={r['n_censored']:,})")

        print(f"\n=== {assay_label} pooled AUROC, strong-binder-only "
              f"(value < {STRONG_BINDER_THRESHOLD_NM:.0f} nM) vs censored ===")
        for metric, r in res["pooled_strong_binder_only"].items():
            print(f"  {metric:<20} AUROC={r['auroc']:.3f}  rank-biserial={r['rank_biserial']:+.3f}  "
                  f"n={r['n']:,} (quant={r['n_quantitative']:,}, censored={r['n_censored']:,})")

        print(f"\n{assay_label}: {len(res['per_allele'])} alleles with n>=10 in both groups "
              f"(primary metric {PRIMARY_METRIC}_best)")

    with open(out_dir / "censored_vs_quantitative_auroc.json", "w") as f:
        json.dump(all_results, f, indent=2)

    per_allele_rows = []
    for assay_label, res in all_results.items():
        for row in res["per_allele"]:
            per_allele_rows.append({"assay": assay_label, **row})
    pd.DataFrame(per_allele_rows).to_csv(out_dir / "censored_vs_quantitative_per_allele.csv", index=False)

    plot_distributions(df, out_dir)
    print(f"\nOutputs written to {out_dir}")


if __name__ == "__main__":
    main()
