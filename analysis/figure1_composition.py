"""
Figure 1 composition panels, regenerated on the merged v2 release
(95 alleles / 112,378 complexes / 118,751 measurement rows).

Panels (A = PyMOL ensemble render, unchanged, not produced here):
  B  pairs per allele, log-scaled, colored by locus (HLA-A/B/C)
  C  peptide-length distribution, split by locus
  D  IC50 affinity distribution, log10, with assay-ceiling lines
  E  KD  affinity distribution, log10, with assay-ceiling lines

All counts are read live from release_v2_final/metadata.csv; the numbers a
reader sees are asserted against that file (never hard-coded).

Usage:   python3 figure1_composition.py --out-dir <dir>
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT  # noqa: E402

METADATA = DATA_ROOT / "release_v2_final/metadata.csv"
IC50_CEIL = (20000.0, 50000.0, 70000.0)
KD_CEIL = (5000.0, 10000.0, 20000.0)
LOCUS_COLOR = {"A": "#4477AA", "B": "#EE6677", "C": "#228833"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    md = pd.read_csv(METADATA, low_memory=False)
    pairs = md[["allele_compact", "peptide", "peptide_length"]].drop_duplicates()
    pairs["locus"] = pairs["allele_compact"].str[0]

    # --- assertions: the figure's headline counts must match the file ---
    assert pairs.shape[0] == 112378, f"expected 112,378 pairs, got {pairs.shape[0]}"
    assert pairs["allele_compact"].nunique() == 95, f"expected 95 alleles, got {pairs['allele_compact'].nunique()}"

    per_allele = pairs.groupby("allele_compact").size().sort_values(ascending=False)
    per_allele_locus = pairs.groupby("allele_compact")["locus"].first()

    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], hspace=0.45, wspace=0.25)

    # Panel B: pairs per allele (all 95), log y, colored by locus
    axB = fig.add_subplot(gs[0, :])
    x = np.arange(len(per_allele))
    colors = [LOCUS_COLOR[per_allele_locus[a]] for a in per_allele.index]
    axB.bar(x, per_allele.values, color=colors, edgecolor="white", linewidth=0.3)
    axB.set_yscale("log")
    axB.set_xticks(x)
    axB.set_xticklabels(per_allele.index, rotation=90, fontsize=5)
    axB.set_ylabel("peptide–allele pairs (log)")
    axB.set_title(f"B  Pairs per allele — {len(per_allele)} alleles, "
                  f"{pairs.shape[0]:,} pairs", fontsize=9, loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in LOCUS_COLOR.values()]
    axB.legend(handles, [f"HLA-{k}" for k in LOCUS_COLOR], fontsize=7, frameon=False)

    # Panel C: peptide length distribution by locus
    axC = fig.add_subplot(gs[1, 0])
    lengths = sorted(pairs["peptide_length"].unique())
    width = 0.8 / len(LOCUS_COLOR)
    for i, (loc, col) in enumerate(LOCUS_COLOR.items()):
        sub = pairs[pairs["locus"] == loc]["peptide_length"].value_counts()
        vals = [sub.get(L, 0) for L in lengths]
        axC.bar(np.arange(len(lengths)) + i * width, vals, width, color=col, label=f"HLA-{loc}")
    axC.set_yscale("log")
    axC.set_xticks(np.arange(len(lengths)) + width)
    axC.set_xticklabels(lengths)
    axC.set_xlabel("peptide length"); axC.set_ylabel("pairs (log)")
    axC.set_title("C  Peptide length by locus", fontsize=9, loc="left")
    axC.legend(fontsize=7, frameon=False)

    # Panels D/E: affinity distributions share the lower-right cell (two stacked axes)
    gsDE = gs[1, 1].subgridspec(2, 1, hspace=0.6)
    for row, (mtype, ceil, label) in enumerate(
            [("IC50", IC50_CEIL, "IC50"), ("Kd", KD_CEIL, "KD")]):
        ax = fig.add_subplot(gsDE[row])
        v = pd.to_numeric(md.loc[md["measurement_type"] == mtype, "measurement_value"],
                          errors="coerce").dropna()
        v = v[v > 0]
        ax.hist(np.log10(v), bins=40, color="#666666", edgecolor="white", linewidth=0.3)
        for c in ceil:
            ax.axvline(np.log10(c), color="#CC3311", ls="--", lw=0.8, alpha=0.7)
        ax.set_title(f"{'D' if row==0 else 'E'}  {label}  (n={len(v):,})",
                     fontsize=8, loc="left")
        ax.set_xlabel("log10(nM)"); ax.set_ylabel("count")

    fig.suptitle(f"Figure 1 (regenerated) — merged v2 dataset: "
                 f"{pairs['allele_compact'].nunique()} alleles, {pairs.shape[0]:,} pairs, "
                 f"{len(md):,} measurements", fontsize=10)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure1_composition_v2.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # write the underlying counts too (provenance for the caption)
    per_allele.rename("n_pairs").reset_index().assign(
        locus=lambda d: d["allele_compact"].str[0]).to_csv(
        out / "figure1_pairs_per_allele.csv", index=False)
    print(f"pairs={pairs.shape[0]:,} alleles={pairs['allele_compact'].nunique()} "
          f"loci={pairs['locus'].value_counts().to_dict()}")
    print(f"wrote figure1_composition_v2.pdf/png + figure1_pairs_per_allele.csv to {out}")


if __name__ == "__main__":
    main()
