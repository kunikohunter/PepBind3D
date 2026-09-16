"""Does pooling the three Kd assays distort the per-allele correlations?

Figure 1 splits Kd by assay method because the three dissociation-constant
assays differ in censoring (the competitive radioligand subset is ~1% censored,
the two fluorescence subsets 71-74%) and sit about one log unit apart. Figure 3C
then pools them: one Spearman rho per allele over all of that allele's Kd
measurements. That is only legitimate if pooling does not move the per-allele
estimate, so this measures it rather than asserting it.

The comparison is, for every allele that genuinely mixes assays (two or more
methods with n >= 10):

    pooled rho   against   the n-weighted mean of the within-assay rho

Why it comes out small: the censored values are excluded before any correlation
is computed, and censoring is the thing that differs most between these assays.
What remains is a location shift, and Spearman ranks within an allele are
largely insensitive to one.

This does NOT say the assays are interchangeable. It says that for this
statistic, on these alleles, pooling and stratifying give the same answer.
Anyone modelling absolute Kd values still has to stratify: see the dataset card.

Usage:
    python3 kd_assay_pooling_check.py --out-dir <dir>
Self-test:
    python3 kd_assay_pooling_check.py --self-test
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import RELEASE_METADATA  # noqa: E402

# Assay detection ceilings for Kd, and the floor at or above which a value is
# also a ceiling report. Same rule as every other correlation in the paper.
KD_CENSORED = {5000, 10000, 20000}
KD_FLOOR = 20000
MIN_N = 10          # minimum measurements for an allele to get a rho
SCORE = "I_sc_best"


def quantitative_kd(df: pd.DataFrame) -> pd.DataFrame:
    """Uncensored Kd rows, with a log10 column, flagged records dropped."""
    df = df[~df["flagged"].astype(bool)]
    kd = df[df["measurement_type"].astype(str).str.lower().str.fullmatch("kd")].copy()
    keep = ~(kd["measurement_value"].isin(KD_CENSORED)
             | (kd["measurement_value"] >= KD_FLOOR))
    kd = kd[keep].copy()
    kd["log_value"] = np.log10(kd["measurement_value"])
    return kd


def compare(kd: pd.DataFrame, min_n: int = MIN_N) -> tuple:
    """Per-allele pooled rho against the within-assay weighted mean."""
    rows = []
    for allele, g in kd.groupby("allele"):
        if len(g) < min_n:
            continue
        counts = g["assay_method"].value_counts()
        pooled, _ = spearmanr(g[SCORE], g["log_value"])
        per = {}
        for method, n in counts.items():
            if n >= min_n:
                sub = g[g["assay_method"] == method]
                r, _ = spearmanr(sub[SCORE], sub["log_value"])
                per[method] = {"rho": float(r), "n": int(n)}
        row = {
            "allele": allele,
            "n": int(len(g)),
            "n_methods_over_min": int(len(per)),
            "dominant_fraction": float(counts.iloc[0] / len(g)),
            "pooled_rho": float(pooled),
        }
        if len(per) >= 2:
            rhos = [v["rho"] for v in per.values()]
            ns = [v["n"] for v in per.values()]
            row["within_assay_weighted_mean"] = float(np.average(rhos, weights=ns))
            row["difference"] = row["pooled_rho"] - row["within_assay_weighted_mean"]
        rows.append(row)
    df = pd.DataFrame(rows)
    mixed = df[df["n_methods_over_min"] >= 2]
    summary = {
        "n_alleles": int(len(df)),
        "n_single_assay": int((df["dominant_fraction"] >= 0.95).sum()),
        "n_mixed": int(len(mixed)),
        "median_abs_difference": float(mixed["difference"].abs().median()),
        "max_abs_difference": float(mixed["difference"].abs().max()),
        "min_n": min_n,
    }
    return summary, df


def self_test() -> None:
    """Two alleles with a known answer.

    `clean` mixes two assays that agree, so pooling must reproduce them.
    `offset` mixes two assays whose values sit a decade apart but rank the same
    way within each: pooling must still recover a positive rho, which is the
    property the real check relies on.
    """
    rng = np.random.default_rng(0)
    n = 60
    score = np.linspace(-80, -55, n)

    def frame(allele, method, values):
        return pd.DataFrame({
            "allele": allele, "assay_method": method, "flagged": False,
            "measurement_type": "Kd", "measurement_value": values,
            SCORE: score,
        })

    val = 10 ** (2 + 0.02 * (score + 80) + rng.normal(0, 0.05, n))
    df = pd.concat([
        frame("A*01:01", "m1", val),
        frame("A*01:01", "m2", val),
        frame("B*01:01", "m1", val),
        frame("B*01:01", "m2", val * 10),     # same ranks, one decade higher
    ], ignore_index=True)

    kd = quantitative_kd(df)
    assert len(kd) == 4 * n, len(kd)
    summary, res = compare(kd)
    assert summary["n_alleles"] == 2, summary
    assert summary["n_mixed"] == 2, summary

    clean = res[res.allele == "A*01:01"].iloc[0]
    assert abs(clean["difference"]) < 1e-9, clean["difference"]

    offset = res[res.allele == "B*01:01"].iloc[0]
    # Clearly positive, but visibly attenuated from the ~0.99 each assay gives
    # on its own: interleaving two offset distributions costs rank agreement.
    assert 0.4 < offset["pooled_rho"] < 0.6, offset["pooled_rho"]
    # A pure location shift does move the pooled estimate, which is exactly why
    # the real data has to be measured rather than argued about.
    assert abs(offset["difference"]) > 0.01, offset["difference"]
    print("Self-test PASSED: identical assays pool exactly; a one-decade offset "
          "keeps the pooled correlation positive but does shift it, so the "
          "effect on real data is an empirical question.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.out_dir:
        raise SystemExit("--out-dir is required")

    self_test()
    print()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    kd = quantitative_kd(pd.read_csv(RELEASE_METADATA, low_memory=False))
    summary, res = compare(kd)

    print(f"quantitative Kd measurements : {len(kd):,}")
    print(f"alleles with n >= {summary['min_n']}          : {summary['n_alleles']}")
    print(f"  effectively single-assay   : {summary['n_single_assay']}")
    print(f"  genuinely mixed            : {summary['n_mixed']}")
    print(f"\npooled rho vs within-assay weighted mean, over the mixed alleles:")
    print(f"  median |difference|        : {summary['median_abs_difference']:.3f}")
    print(f"  largest |difference|       : {summary['max_abs_difference']:.3f}")

    res.to_csv(out / "kd_assay_pooling_per_allele.csv", index=False)
    with open(out / "kd_assay_pooling.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nwrote kd_assay_pooling.{{csv,json}} to {out}")


if __name__ == "__main__":
    main()
