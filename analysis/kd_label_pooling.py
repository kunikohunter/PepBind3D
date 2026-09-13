"""
Phase 1 item 1.3 (REVISION_PLAN.md): KD label pooling.

Three IEDB assay-response labels were *intended* to be normalized to a single
"dissociation constant (KD)" category during curation (IEDBTestPipeline.py:136):
  - "dissociation constant KD (~EC50)"
  - "dissociation constant KD"
  - "dissociation constant KD (~IC50)"  (IEDBTestPipeline.py:136 has a typo,
    "dissociation constant (~IC50)" missing "KD" -- that string never occurs
    in the raw IEDB data, so this label's rows never match the pipeline's
    normalization and are dropped from the curated dataset entirely, not
    pooled. See the 2026-09-04 note by RAW_LABELS below.)

metadata.csv only stores the post-normalization label, so the original
three-way split has to be recovered from the raw IEDB bulk download and
joined back onto the curated KD rows. This script does that recovery,
reports per-label medians/IQRs, and runs a Kruskal-Wallis test (+ pairwise
KS tests) on log10(measurement_value) across the three groups.

Usage:
    $ python3 kd_label_pooling.py --out-dir <output dir>

Self-test:
    $ python3 kd_label_pooling.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, ks_2samp

# NB: the working-tree copy at data/IEDB_data_clean/metadata.csv is a stale
# (Mar 2026) pre-rebuild file missing the six per-metric score columns; the
# huggingface-staged copy is the one all validation notebooks actually read
# (HF_DIR/metadata.csv) and is what's used here. Row-level curation (allele,
# peptide, measurement_type/value, flagged) was confirmed byte-identical
# between the two files, so this only matters for analyses touching scores.
#
# 2026-09-04 fix: HF_DIR/metadata.csv is not present on this filesystem (see
# 08_revision_figures.ipynb, same finding). The file is public and un-gated
# on HuggingFace (kunikohunter/PepBind3D), so METADATA_FN now falls back to
# a local, read-only download of it (fetched via huggingface_hub.hf_hub_download,
# never written into any huggingface/ or released-dataset directory) if the
# authors' hardcoded cluster path isn't found. RAW_IEDB_FN is unchanged and
# was found in place.
#
# 2026-09-13: the release is now the merged v2 set (95 alleles, 112,378 pairs,
# 118,751 measurement rows, HLA-C included), so the default metadata is
# release_v2_final/metadata.csv -- 97,574 unflagged KD rows against v1's
# 35,851. The v1-scope run stays reproducible via --metadata.
METADATA_FN = "<HOME>/main_project/data/IEDB_data_clean/release_v2_final/metadata.csv"
METADATA_V1_FN = "<HOME>/main_project/data/IEDB_data_clean/huggingface/metadata.csv"
RAW_IEDB_FN = "<HOME>/Data/MHC_database/build/mhc_ligand_full.csv"

# 2026-09-04 fix: the third label was "dissociation constant (~IC50)" (missing
# "KD"), copied from the same typo in IEDBTestPipeline.py:136's normalization
# dict. Scanning the raw IEDB bulk file's actual "Assay Response measured"
# values shows the real label is "dissociation constant KD (~IC50)" -- the
# typo'd string never occurs in the raw data. This means the production
# curation pipeline's `.replace()` for that label silently never fires: rows
# with this raw label are left as "dissociation constant KD (~IC50)", which
# isn't in `valid_responses`, so they are dropped entirely rather than pooled
# into the "Kd" category. Corrected here so the join recovers the true label
# set; see the run report for what this means for the pooling question.
RAW_LABELS = [
    "dissociation constant KD (~EC50)",
    "dissociation constant KD",
    "dissociation constant KD (~IC50)",
]


def flatten_columns(df):
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = [" ".join(col).strip() for col in df.columns]
    return df


def normalize_pmid(s):
    """PMIDs arrive as floats whenever the column holds any NaN, so the same
    reference is "22508927.0" in one file and "22508927" in the other. Strip a
    trailing ".0" on both sides or the join silently matches nothing."""
    return (s.astype(str).str.strip()
             .str.replace(r"\.0$", "", regex=True))


def load_curated_kd_rows(metadata_fn):
    df = pd.read_csv(metadata_fn, low_memory=False)
    kd = df[df["measurement_type"] == "Kd"].copy()
    # The merged v2 metadata has no `flagged` column: flagged rows were dropped
    # during reconciliation rather than carried with a marker. Only filter when
    # the column is actually present (v1 huggingface/metadata.csv).
    if "flagged" in kd.columns:
        kd = kd[kd["flagged"] == False]  # noqa: E712
    kd["measurement_value"] = pd.to_numeric(kd["measurement_value"], errors="coerce")
    kd = kd.dropna(subset=["measurement_value"])
    cols = ["allele_iedb", "peptide", "measurement_value", "pubmed_id"]
    # merged v2 metadata only; lets the report attribute each label to the
    # release batch that produced it (see the typo note in main()).
    if "source_version" in kd.columns:
        cols.append("source_version")
    return kd[cols].copy()


def scan_raw_labels(raw_fn, alleles, chunksz=500_000, verbose=True):
    """Stream the raw IEDB bulk download, keep rows whose raw assay-response
    label is one of the three KD variants AND whose allele appears in the
    curated release, and return the minimal columns needed to join back onto
    metadata.csv.

    `alleles` is the set of curated allele keys in filesystem-safe form
    (e.g. "HLA-B_35_03"). Filtering on release membership rather than on an
    "HLA-A*"/"HLA-B*" substring is what makes this correct for the merged v2
    release, which added HLA-C: a hardcoded locus filter would drop every
    HLA-C row before the join and report them as unmatched.
    """
    keep = []
    with pd.read_csv(raw_fn, header=[0, 1], chunksize=chunksz, low_memory=False) as reader:
        for i, chunk in enumerate(reader):
            chunk = flatten_columns(chunk)
            allele_key = (chunk["MHC Restriction Name"].astype(str)
                          .str.replace(":", "_", regex=False)
                          .str.replace("*", "_", regex=False))
            in_release = allele_key.isin(alleles)
            is_kd_variant = chunk["Assay Response measured"].isin(RAW_LABELS)
            sub = chunk.loc[in_release & is_kd_variant, [
                "MHC Restriction Name",
                "Epitope Name",
                "Assay Quantitative measurement",
                "Reference PMID",
                "Assay Response measured",
            ]].copy()
            keep.append(sub)
            if verbose:
                print(f"  chunk {i+1}: {len(sub)} matching rows (cumulative {sum(len(k) for k in keep)})",
                      file=sys.stderr)
    raw = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    raw = raw.rename(columns={
        "MHC Restriction Name": "allele_iedb",
        "Epitope Name": "peptide",
        "Assay Quantitative measurement": "measurement_value",
        "Reference PMID": "pubmed_id",
        "Assay Response measured": "raw_label",
    })
    # metadata.csv stores alleles in the filesystem-safe form used for per-allele
    # output directories (IEDBTestPipeline.py: a.replace(":","_").replace("*","_")),
    # not the raw IEDB "HLA-B*35:03" form. Apply the same transform here so the
    # join key matches.
    raw["allele_iedb"] = raw["allele_iedb"].str.replace(":", "_", regex=False).str.replace("*", "_", regex=False)
    raw["measurement_value"] = pd.to_numeric(raw["measurement_value"], errors="coerce")
    return raw.dropna(subset=["measurement_value"])


def join_labels(curated_kd, raw_kd):
    """Join curated KD rows to their original raw label on
    (allele_iedb, peptide, measurement_value, pubmed_id). Rows with an
    ambiguous (many-to-one) or missing match are reported, not silently
    dropped."""
    for df in (curated_kd, raw_kd):
        df["pubmed_id"] = normalize_pmid(df["pubmed_id"])
        df["peptide"] = df["peptide"].astype(str).str.strip()
        df["allele_iedb"] = df["allele_iedb"].astype(str).str.strip()

    raw_dedup = raw_kd.drop_duplicates(subset=["allele_iedb", "peptide", "measurement_value", "pubmed_id"])
    dup_keys = raw_kd.duplicated(subset=["allele_iedb", "peptide", "measurement_value", "pubmed_id"], keep=False)
    n_ambiguous_raw_groups = raw_kd[dup_keys].drop_duplicates(
        subset=["allele_iedb", "peptide", "measurement_value", "pubmed_id"]
    ).shape[0]

    merged = curated_kd.merge(
        raw_dedup[["allele_iedb", "peptide", "measurement_value", "pubmed_id", "raw_label"]],
        on=["allele_iedb", "peptide", "measurement_value", "pubmed_id"],
        how="left",
    )
    return merged, n_ambiguous_raw_groups


def summarize(merged):
    matched = merged.dropna(subset=["raw_label"])
    rows = []
    for label in RAW_LABELS:
        vals = matched.loc[matched["raw_label"] == label, "measurement_value"].astype(float)
        log_vals = np.log10(vals[vals > 0])
        rows.append({
            "raw_label": label,
            "n": int(len(log_vals)),
            "median_nM": float(vals.median()) if len(vals) else float("nan"),
            "log10_median": float(log_vals.median()) if len(log_vals) else float("nan"),
            "log10_q1": float(log_vals.quantile(0.25)) if len(log_vals) else float("nan"),
            "log10_q3": float(log_vals.quantile(0.75)) if len(log_vals) else float("nan"),
        })
    return pd.DataFrame(rows), matched


# KD assay detection ceilings (same values as the correlation analyses; see
# CLAUDE.md "Conventions"). Needed here because the two true-KD labels are
# 71-74% censored while the (~IC50) label is 1% censored, so a pooling test on
# raw values measures the difference in censoring rate, not a difference in
# reported affinity. Both versions of the test are reported.
KD_CEILINGS = (5000.0, 10000.0, 20000.0)


def stratify_by_censoring(matched):
    """Per-label censoring rate and quantitative-only log10 distribution.

    This is the comparison that answers the pooling question: if the labels
    differ only in how often they hit the ceiling, pooling is defensible; if
    they still differ once censored values are removed, they are measuring
    different things and must not be pooled silently."""
    m = matched.copy()
    m["measurement_value"] = m["measurement_value"].astype(float)
    m["censored"] = m["measurement_value"].isin(KD_CEILINGS)
    rows = []
    for label in RAW_LABELS:
        sub = m[m["raw_label"] == label]
        q = sub[(~sub["censored"]) & (sub["measurement_value"] > 0)]
        lg = np.log10(q["measurement_value"])
        rows.append({
            "raw_label": label,
            "n": int(len(sub)),
            "n_censored": int(sub["censored"].sum()),
            "censored_frac": float(sub["censored"].mean()) if len(sub) else float("nan"),
            "n_quantitative": int(len(lg)),
            "quant_log10_median": float(lg.median()) if len(lg) else float("nan"),
            "quant_log10_q1": float(lg.quantile(0.25)) if len(lg) else float("nan"),
            "quant_log10_q3": float(lg.quantile(0.75)) if len(lg) else float("nan"),
        })
    return pd.DataFrame(rows), m[~m["censored"] & (m["measurement_value"] > 0)]


def run_tests(matched):
    groups = [
        np.log10(matched.loc[matched["raw_label"] == label, "measurement_value"].astype(float))
        for label in RAW_LABELS
    ]
    groups = [g[g.notna()] for g in groups]
    result = {}
    if all(len(g) >= 2 for g in groups):
        stat, p = kruskal(*groups)
        result["kruskal_wallis"] = {"H": float(stat), "p": float(p), "n_groups": len(groups)}
    pairwise = {}
    for i in range(len(RAW_LABELS)):
        for j in range(i + 1, len(RAW_LABELS)):
            if len(groups[i]) >= 2 and len(groups[j]) >= 2:
                stat, p = ks_2samp(groups[i], groups[j])
                pairwise[f"{RAW_LABELS[i]} vs {RAW_LABELS[j]}"] = {"D": float(stat), "p": float(p)}
    result["pairwise_ks"] = pairwise
    return result


def self_test():
    """Construct synthetic data with a known analytic answer: three groups
    where two are drawn from the same distribution and one is shifted by
    1 log10 unit (10x). Kruskal-Wallis and the KS tests must both flag the
    shifted group as significantly different and the two identical groups
    as not significantly different."""
    rng = np.random.default_rng(0)
    same_a = 10 ** rng.normal(2.0, 0.3, size=500)
    same_b = 10 ** rng.normal(2.0, 0.3, size=500)
    shifted = 10 ** rng.normal(3.0, 0.3, size=500)  # 10x shift in nM

    fake = pd.DataFrame({
        "measurement_value": np.concatenate([same_a, same_b, shifted]),
        "raw_label": (
            [RAW_LABELS[0]] * len(same_a) + [RAW_LABELS[1]] * len(same_b) + [RAW_LABELS[2]] * len(shifted)
        ),
    })
    result = run_tests(fake)
    kw_p = result["kruskal_wallis"]["p"]
    assert kw_p < 1e-6, f"self-test FAILED: Kruskal-Wallis should detect the shifted group (p={kw_p})"

    pw = result["pairwise_ks"]
    same_vs_same_p = pw[f"{RAW_LABELS[0]} vs {RAW_LABELS[1]}"]["p"]
    same_vs_shifted_p = pw[f"{RAW_LABELS[0]} vs {RAW_LABELS[2]}"]["p"]
    assert same_vs_same_p > 0.01, f"self-test FAILED: identical groups falsely significant (p={same_vs_same_p})"
    assert same_vs_shifted_p < 1e-6, f"self-test FAILED: shifted group not detected (p={same_vs_shifted_p})"

    print("Self-test PASSED: Kruskal-Wallis and pairwise KS correctly detect a known 10x shift "
          "and correctly find no difference between two identical distributions.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=str, default=None,
                     help="Directory to write outputs (label_recovery.csv, kd_label_summary.csv, "
                          "kd_label_tests.json). Required unless --self-test.")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--metadata", type=str, default=METADATA_FN,
                    help=f"metadata.csv to analyse. Default: the merged v2 release "
                         f"({METADATA_FN}). Pass {METADATA_V1_FN} to reproduce the v1-scope run.")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    if args.out_dir is None:
        raise SystemExit("--out-dir is required (never write outputs into the scripts directory)")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading curated KD rows from {args.metadata} ...", file=sys.stderr)
    curated_kd = load_curated_kd_rows(args.metadata)
    alleles = set(curated_kd["allele_iedb"].astype(str).str.strip())
    print(f"  {len(curated_kd)} curated, unflagged KD rows across {len(alleles)} alleles",
          file=sys.stderr)

    print(f"Scanning raw IEDB bulk file for original KD-variant labels ({RAW_IEDB_FN})...", file=sys.stderr)
    raw_kd = scan_raw_labels(RAW_IEDB_FN, alleles)
    print(f"  {len(raw_kd)} raw rows with a KD-variant label in a released allele", file=sys.stderr)

    merged, n_ambiguous = join_labels(curated_kd, raw_kd)
    n_matched = merged["raw_label"].notna().sum()
    n_total = len(merged)
    print(f"Join: {n_matched}/{n_total} curated KD rows matched to a raw label "
          f"({100*n_matched/n_total:.1f}%); {n_ambiguous} ambiguous raw (allele,peptide,value,pubmed) "
          f"groups collapsed via drop_duplicates before the join.", file=sys.stderr)

    merged.to_csv(out_dir / "kd_label_recovery.csv", index=False)

    summary_df, matched = summarize(merged)
    summary_df.to_csv(out_dir / "kd_label_summary.csv", index=False)
    print(summary_df.to_string(index=False))

    # Which release version contributed each label. This is how the
    # IEDBTestPipeline.py:136 typo shows up in the data: the "(~IC50)" label has
    # zero v1 rows (the .replace() never fired, so those records were dropped)
    # and 42,551 v2 rows (the v2 curation run matches the real IEDB string). The
    # released v1 dataset is missing them; the merged release contains them.
    if "source_version" in matched.columns:
        by_version = matched.groupby(["raw_label", "source_version"]).size().unstack(fill_value=0)
        by_version.to_csv(out_dir / "kd_label_by_source_version.csv")
        print("\nraw_label x source_version:\n" + by_version.to_string())

    censor_df, quant = stratify_by_censoring(matched)
    censor_df.to_csv(out_dir / "kd_label_censoring.csv", index=False)
    print("\nper-label censoring and quantitative-only distribution:")
    print(censor_df.to_string(index=False))

    quant_tests = run_tests(quant)
    test_results = run_tests(matched)
    test_results["quantitative_only"] = quant_tests
    import json
    with open(out_dir / "kd_label_tests.json", "w") as f:
        json.dump({
            "metadata_file": str(args.metadata),
            "n_alleles": len(alleles),
            "n_curated_kd_rows": int(n_total),
            "n_matched": int(n_matched),
            "match_rate": float(n_matched / n_total),
            "n_ambiguous_raw_groups": int(n_ambiguous),
            **test_results,
        }, f, indent=2)
    print(json.dumps(test_results, indent=2))

    plot_distributions(matched, out_dir)
    print(f"\nOutputs written to {out_dir}")


def plot_distributions(matched, out_dir):
    # 2026-09-04: house-style fix so this figure is a drop-in match for
    # Figure 1E (04_figure1_panels.ipynb Panel D, Kd sub-panel) and the
    # revision figures (08_revision_figures.ipynb): utils.set_plot_style(),
    # Scientific-Data mm sizing, 300 dpi PNG + PDF, same #4477AA/#CC3311
    # palette. Previously this used ad hoc inch sizing and 150 dpi PNG.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # NB: `from utils import set_plot_style` pulls in utils/__init__.py, which
    # eagerly imports utils.structure -> Bio.PDB (biopython), not installed in
    # this environment (module load python/3.11.5 scipy-stack/2025a). Rather
    # than add a biopython dependency to a KD-pooling analysis script, inline
    # utils.plotting.set_plot_style()'s exact rcParams so the output still
    # matches house style byte-for-byte.
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "lines.linewidth": 1.5,
    })

    MM = 1 / 25.4
    COL2_WIDTH = 183 * MM  # Scientific Data double-column width, matches 04/08 notebooks

    # Log y-axis (as in 04_figure1_panels.ipynb Panel D/E) -- the raw counts
    # are dominated by the 20,000 nM detection ceiling, so a linear axis
    # hides everything else in the distribution.
    censor_values_nM = [5000, 10000, 20000]
    fig, axes = plt.subplots(1, 3, figsize=(COL2_WIDTH, COL2_WIDTH / 2.6), sharey=True)
    bins = np.arange(0, 5.05, 0.1)  # log10 nM, 0 to ~100 uM
    for ax, label in zip(axes, RAW_LABELS):
        vals = matched.loc[matched["raw_label"] == label, "measurement_value"].astype(float)
        vals = vals[vals > 0]
        n = len(vals)
        if n > 0:
            ax.hist(np.log10(vals), bins=bins, color="#4477AA", edgecolor="white", linewidth=0.4, zorder=2)
        else:
            ax.text(0.5, 0.5, "n = 0\n(dropped in curation,\nnot pooled -- see report)",
                    ha="center", va="center", fontsize=6.5, color="#555555",
                    transform=ax.transAxes)
        for cv in censor_values_nM:
            ax.axvline(np.log10(cv), color="#CC3311", linestyle="--", linewidth=1, alpha=0.7, zorder=10)
        ax.set_title(f"{label}\n(n={n:,})", fontsize=8)
        ax.set_xlabel("log10(KD, nM)")
        ax.set_yscale("log")
        ax.set_ylim(0.8, 40000)
        ax.minorticks_off()
        ax.yaxis.grid(True, color="#dddddd", linewidth=0.5, alpha=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("count (log scale)")
    fig.suptitle("KD label pooling: raw IEDB assay-response labels normalized to 'dissociation constant (KD)'"
                 " in metadata.csv (dashed red = Fig. 1E assay detection ceilings 5,000 / 10,000 / 20,000 nM)",
                 fontsize=7.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out_dir / "kd_label_distributions.pdf")
    fig.savefig(out_dir / "kd_label_distributions.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
