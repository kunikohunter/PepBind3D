"""
Phase 1 item 1.3 (REVISION_PLAN.md): KD label pooling.

Three IEDB assay-response labels were normalized to a single "dissociation
constant (KD)" category during curation (IEDBTestPipeline.py:136):
  - "dissociation constant KD (~EC50)"
  - "dissociation constant KD"
  - "dissociation constant (~IC50)"

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
METADATA_FN = "/home/huntek1/main_project/data/IEDB_data_clean/huggingface/metadata.csv"
RAW_IEDB_FN = "/home/huntek1/Data/MHC_database/build/mhc_ligand_full.csv"

RAW_LABELS = [
    "dissociation constant KD (~EC50)",
    "dissociation constant KD",
    "dissociation constant (~IC50)",
]


def flatten_columns(df):
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = [" ".join(col).strip() for col in df.columns]
    return df


def load_curated_kd_rows(metadata_fn):
    df = pd.read_csv(metadata_fn)
    kd = df[(df["measurement_type"] == "Kd") & (df["flagged"] == False)].copy()  # noqa: E712
    kd["measurement_value"] = pd.to_numeric(kd["measurement_value"], errors="coerce")
    kd = kd.dropna(subset=["measurement_value"])
    return kd[["allele_iedb", "peptide", "measurement_value", "pubmed_id"]].copy()


def scan_raw_labels(raw_fn, chunksz=500_000, verbose=True):
    """Stream the raw IEDB bulk download, keep only HLA-A/B rows whose raw
    assay-response label is one of the three KD variants, and return the
    minimal columns needed to join back onto metadata.csv."""
    keep = []
    with pd.read_csv(raw_fn, header=[0, 1], chunksize=chunksz, low_memory=False) as reader:
        for i, chunk in enumerate(reader):
            chunk = flatten_columns(chunk)
            allele_col = chunk["MHC Restriction Name"]
            is_ab = allele_col.str.contains("HLA-A*", regex=False, na=False) | allele_col.str.contains(
                "HLA-B*", regex=False, na=False
            )
            is_kd_variant = chunk["Assay Response measured"].isin(RAW_LABELS)
            sub = chunk.loc[is_ab & is_kd_variant, [
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
        df["pubmed_id"] = df["pubmed_id"].astype(str).str.strip()
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
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    if args.out_dir is None:
        raise SystemExit("--out-dir is required (never write outputs into the scripts directory)")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading curated KD rows from metadata.csv...", file=sys.stderr)
    curated_kd = load_curated_kd_rows(METADATA_FN)
    print(f"  {len(curated_kd)} curated, unflagged KD rows", file=sys.stderr)

    print(f"Scanning raw IEDB bulk file for original KD-variant labels ({RAW_IEDB_FN})...", file=sys.stderr)
    raw_kd = scan_raw_labels(RAW_IEDB_FN)
    print(f"  {len(raw_kd)} raw HLA-A/B rows with a KD-variant label", file=sys.stderr)

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

    test_results = run_tests(matched)
    import json
    with open(out_dir / "kd_label_tests.json", "w") as f:
        json.dump({
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
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    censor_values_nM = [5000, 10000, 20000]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3), sharey=True)
    bins = np.arange(0, 5.05, 0.1)  # log10 nM, 0 to ~100 uM
    for ax, label in zip(axes, RAW_LABELS):
        vals = matched.loc[matched["raw_label"] == label, "measurement_value"].astype(float)
        vals = vals[vals > 0]
        n = len(vals)
        if n > 0:
            ax.hist(np.log10(vals), bins=bins, color="#4477AA", edgecolor="white", linewidth=0.4)
        for cv in censor_values_nM:
            ax.axvline(np.log10(cv), color="#CC3311", linestyle="--", linewidth=1, alpha=0.7)
        ax.set_title(f"{label}\n(n={n})", fontsize=9)
        ax.set_xlabel("log10(KD, nM)")
    axes[0].set_ylabel("count")
    fig.suptitle("KD label pooling: three raw IEDB assay-response labels normalized to 'dissociation constant (KD)'"
                 " (dashed red = assay detection ceilings 5,000 / 10,000 / 20,000 nM)", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out_dir / "kd_label_distributions.pdf")
    fig.savefig(out_dir / "kd_label_distributions.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
