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


MIN_N_PER_ALLELE = 10          # Table S3 / Figure 3C convention
WELL_SAMPLED_N = 50            # the n>=50 labelling convention in Figure 3C


def run_per_allele(df, out_dir):
    """Per-allele Spearman between best-decoy I_sc and log affinity -- the
    Supplementary Table S3 / Figure 3C content, regenerated for 95 alleles.

    Same convention as the submitted paper: one row per (allele, assay) with at
    least MIN_N_PER_ALLELE quantitative measurements, censored values excluded.
    Also returns the summary statistics the Results text quotes (median rho,
    IQR, how many alleles are significant, and the strongest alleles among the
    well-sampled ones) so those sentences do not have to be read off a figure."""
    rows = []
    for assay, ceil, floor in [("IC50", IC50_CEILINGS, IC50_FLOOR),
                               ("Kd", KD_CEILINGS, KD_FLOOR)]:
        sub = df[df["measurement_type"] == assay]
        cens = is_censored(sub["measurement_value"], ceil, floor)
        q = sub[(~cens) & (sub["measurement_value"] > 0)].dropna(
            subset=[f"{PRIMARY_METRIC}_best"])
        for allele, g in q.groupby("allele"):
            if len(g) < MIN_N_PER_ALLELE:
                continue
            r = spearmanr(g[f"{PRIMARY_METRIC}_best"], np.log10(g["measurement_value"]))
            rows.append({"assay": assay, "allele": allele, "n": int(len(g)),
                         "rho": float(r.statistic), "p": float(r.pvalue),
                         "significant_p05": bool(r.pvalue < 0.05)})
    per_allele = pd.DataFrame(rows).sort_values(["assay", "n"], ascending=[True, False])

    # Holm-Bonferroni within each assay. Needed for one specific claim: the
    # submitted paper says no allele shows a significant negative correlation.
    # On 112k one does at nominal alpha (B*15:42, KD, n=16, rho=-0.50,
    # p=0.047), so the claim has to be stated against corrected p-values, where
    # it still holds.
    per_allele["p_holm"] = np.nan
    for assay, g in per_allele.groupby("assay"):
        g = g.sort_values("p")
        m = len(g)
        adj, running = [], 0.0
        for rank, p in enumerate(g["p"].values):
            running = max(running, min(1.0, p * (m - rank)))
            adj.append(running)          # enforce monotonicity, as Holm requires
        per_allele.loc[g.index, "p_holm"] = adj
    per_allele["significant_holm"] = per_allele["p_holm"] < 0.05
    per_allele.to_csv(out_dir / "per_allele_spearman_112k.csv", index=False)

    summary = {}
    print("\n=== per-allele Spearman (Table S3 / Figure 3C), 112k ===")
    for assay, g in per_allele.groupby("assay"):
        ws = g[g["n"] >= WELL_SAMPLED_N]
        top = ws.nlargest(3, "rho")[["allele", "rho", "p", "n"]].to_dict("records")
        bot = ws.nsmallest(1, "rho")[["allele", "rho", "p", "n"]].to_dict("records")
        summary[assay] = {
            "n_alleles": int(len(g)),
            "n_significant": int(g["significant_p05"].sum()),
            "median_rho": float(g["rho"].median()),
            "iqr": [float(g["rho"].quantile(.25)), float(g["rho"].quantile(.75))],
            "n_well_sampled": int(len(ws)),
            "strongest_positive": top,
            "most_negative": bot,
            "n_significant_holm": int(g["significant_holm"].sum()),
            "n_significant_negative": int(((g["rho"] < 0) & g["significant_p05"]).sum()),
            "n_significant_negative_holm": int(
                ((g["rho"] < 0) & g["significant_holm"]).sum()),
            "nominal_negative_alleles": g.loc[
                (g["rho"] < 0) & g["significant_p05"],
                ["allele", "n", "rho", "p", "p_holm"]].to_dict("records"),
        }
        s = summary[assay]
        print(f"  {assay}: {s['n_alleles']} alleles (n>={MIN_N_PER_ALLELE}), "
              f"median rho={s['median_rho']:.2f} "
              f"IQR {s['iqr'][0]:.2f}-{s['iqr'][1]:.2f}, "
              f"{s['n_significant']} significant at p<0.05 "
              f"({s['n_significant_holm']} after Holm)")
        print(f"    strongest (n>={WELL_SAMPLED_N}): " + ", ".join(
            f"{t['allele']} rho={t['rho']:.2f} n={t['n']}" for t in top))
        print(f"    most negative: " + ", ".join(
            f"{t['allele']} rho={t['rho']:.2f} p={t['p']:.2g} n={t['n']}" for t in bot)
            + f"; nominally significant negatives: {s['n_significant_negative']}"
              f", after Holm: {s['n_significant_negative_holm']}")
        for a in s["nominal_negative_alleles"]:
            print(f"      {a['allele']} n={a['n']} rho={a['rho']:.2f} "
                  f"p={a['p']:.3f} p_holm={a['p_holm']:.2f}")
    with open(out_dir / "per_allele_spearman_112k_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    return per_allele, summary


def run_composition(df, out_dir):
    """Supplementary Table S5: per-allele composition -- unique peptides, pairs,
    IC50/KD measurement counts, peptide-length range, and share of the dataset."""
    pairs = df[["allele", "allele_compact", "peptide", "peptide_length"]].drop_duplicates(
        subset=["allele_compact", "peptide"])
    total_pairs = len(pairs)
    g = pairs.groupby("allele")
    comp = pd.DataFrame({
        "n_unique_peptides": g["peptide"].nunique(),
        "n_pairs": g.size(),
        "peptide_len_min": g["peptide_length"].min(),
        "peptide_len_max": g["peptide_length"].max(),
    })
    counts = (df.groupby(["allele", "measurement_type"]).size()
                .unstack(fill_value=0).rename(columns={"IC50": "n_IC50", "Kd": "n_KD"}))
    comp = comp.join(counts, how="left").fillna({"n_IC50": 0, "n_KD": 0})
    comp["pct_of_pairs"] = 100 * comp["n_pairs"] / total_pairs
    comp = comp.sort_values("n_pairs", ascending=False).reset_index()
    comp.to_csv(out_dir / "per_allele_composition_112k.csv", index=False)

    top = comp.iloc[0]
    few = comp[comp["n_pairs"] < 5]["allele"].tolist()
    single_assay = comp[(comp["n_IC50"] == 0) | (comp["n_KD"] == 0)]
    print("\n=== per-allele composition (Table S5), 112k ===")
    print(f"  {len(comp)} alleles, {total_pairs:,} pairs")
    print(f"  most represented: {top['allele']} with {int(top['n_pairs']):,} pairs "
          f"({top['pct_of_pairs']:.1f}%)")
    print(f"  alleles with only one assay type: {len(single_assay)}")
    print(f"  alleles with fewer than 5 pairs ({len(few)}): {', '.join(few) if few else 'none'}")
    return comp


def run_affinity_by_kd_label(df, kd_labels_csv, out_dir):
    """Split the KD score-affinity correlation by the ORIGINAL IEDB assay-response
    label recovered by kd_label_pooling.py.

    Necessary because the three labels pooled into "Kd" are not one population:
    the two true-KD labels are 71-74% censored, while "dissociation constant KD
    (~IC50)" is 1% censored and sits ~1 log unit weaker. A single pooled rho over
    all 48,311 quantitative KD rows therefore averages over two different assay
    populations, and the response letter should report the split."""
    rec = pd.read_csv(kd_labels_csv, low_memory=False).dropna(subset=["raw_label"])
    key = ["allele_iedb", "peptide", "measurement_value"]
    for d in (rec, df):
        d["allele_iedb"] = d["allele_iedb"].astype(str).str.strip()
        d["peptide"] = d["peptide"].astype(str).str.strip()
        d["measurement_value"] = pd.to_numeric(d["measurement_value"], errors="coerce")
    rec = rec.drop_duplicates(subset=key)[key + ["raw_label"]]

    kd = df[df["measurement_type"] == "Kd"].merge(rec, on=key, how="left")
    cens = is_censored(kd["measurement_value"], KD_CEILINGS, KD_FLOOR)
    q = kd[(~cens) & (kd["measurement_value"] > 0)].dropna(subset=[f"{PRIMARY_METRIC}_best"])
    rows = []
    for label, sub in q.groupby("raw_label"):
        if len(sub) < 10:
            continue
        r = spearmanr(sub[f"{PRIMARY_METRIC}_best"], np.log10(sub["measurement_value"]))
        rows.append({"raw_label": label, "n": int(len(sub)),
                     "rho": float(r.statistic), "p": float(r.pvalue)})
    res = pd.DataFrame(rows).sort_values("n", ascending=False)
    res.to_csv(out_dir / "score_affinity_112k_by_kd_label.csv", index=False)
    n_unlabelled = int(q["raw_label"].isna().sum())
    print(f"\n=== KD score-affinity split by original IEDB label "
          f"({PRIMARY_METRIC}_best, quantitative; {n_unlabelled:,} rows unlabelled) ===")
    print(res.to_string(index=False))
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True)
    ap.add_argument("--kd-labels", default=None,
                    help="kd_label_recovery.csv from kd_label_pooling.py; if given, "
                         "the KD correlation is also reported split by original label.")
    args = ap.parse_args(); out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    df = build_scored_metadata()
    cens = run_censored(df, out)
    aff = run_affinity(df, out)
    run_per_allele(df, out)
    run_composition(df, out)
    if args.kd_labels:
        run_affinity_by_kd_label(df, args.kd_labels, out)
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
