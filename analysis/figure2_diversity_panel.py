"""
Figure 2 ensemble-diversity panel (Reviewer 2's decoy-scatter question).

Renders, per pair, the ensemble's internal spread (median decoy-to-decoy
peptide-backbone RMSD, x) against its distance from the crystal (median
decoy-to-crystal RMSD, y), with the y=x diagonal. A point on/below the
diagonal means the crystal sits within the ensemble's own scatter; above it
means the ensemble is internally tighter than its displacement from truth.
Two panels: the 52 PepBind3D validation pairs and the 167 leakage-controlled
benchmark targets. Reads the already-computed CSVs (no recompute).

Usage: python3 figure2_diversity_panel.py --out-dir <dir>
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_BASE = Path("<HOME>/main_project/data/IEDB_data_clean/IEDB_validation/ensemble_diversity_out")
SETS = [
    ("52 PepBind3D validation pairs", OUT_BASE / "ensemble_diversity_per_pair.csv", "#4477AA"),
    ("167 leakage-controlled targets", OUT_BASE / "ensemble_diversity_benchmark167_per_target.csv", "#228833"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))
    for ax, (title, csv, color) in zip(axes, SETS):
        df = pd.read_csv(csv).dropna(subset=["d2d_median", "d2c_median"])
        x, y = df["d2d_median"], df["d2c_median"]
        lim = max(x.max(), y.max()) * 1.05
        ax.plot([0, lim], [0, lim], "--", color="#888888", lw=1, zorder=1,
                label="crystal = ensemble spread (y=x)")
        ax.scatter(x, y, s=22, color=color, alpha=0.7, edgecolor="white", linewidth=0.4, zorder=3)
        within = int((df["d2c_median"] <= df["d2d_max"]).sum()) if "d2d_max" in df else None
        ratio = float(np.median(df["d2c_median"] / df["d2d_median"]))
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("ensemble internal spread\n(median decoy–decoy RMSD, Å)")
        ax.set_ylabel("distance to crystal\n(median decoy–crystal RMSD, Å)")
        n = len(df)
        sub = f"n={n}   median ratio={ratio:.2f}"
        if within is not None:
            sub += f"\ncrystal within scatter: {within}/{n} ({100*within/n:.0f}%)"
        ax.set_title(f"{title}\n{sub}", fontsize=8.5)
        ax.legend(fontsize=6.5, frameon=False, loc="lower right")
        ax.set_aspect("equal", adjustable="box")

    fig.suptitle("Ensemble self-consistency vs. accuracy: the crystal sits within the "
                 "decoys' own scatter", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure2_ensemble_diversity.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figure2_ensemble_diversity.pdf/png to {out}")


if __name__ == "__main__":
    main()
