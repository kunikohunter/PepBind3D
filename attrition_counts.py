"""Per-stage record attrition for the IEDB pHLA curation pipeline.

Reproduces the filtering in IEDBTestPipeline.process_peplist against the raw
IEDB bulk download, logging record counts after each stage. Produces the
Supplementary Table S1 funnel (Rocco comment 5).

This mirrors the pipeline logic; it does not re-run modeling. The final row
should match the released dataset (49,488 measurements / 49,268 pairs). Any
gap indicates a filter applied outside this script.

Usage:
    python attrition_counts.py <mhc_ligand_full.csv> <final_metadata.csv> [out.csv]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


VALID_RESPONSES_RAW = [
    "dissociation constant KD (~EC50)",
    "dissociation constant KD",
    "dissociation constant (~IC50)",
    "dissociation constant (KD)",
    "half maximal inhibitory concentration (IC50)",
]
# After standardization these collapse to two classes:
KD_STD   = "dissociation constant (KD)"
IC50_STD = "half maximal inhibitory concentration (IC50)"


def flatten_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Join the 2-row IEDB header into single strings, matching the pipeline."""
    df = df.copy()
    df.columns = [
        c[1] if (isinstance(c, tuple) and str(c[0]).startswith("Unnamed"))
        else (" ".join(str(x) for x in c).strip() if isinstance(c, tuple) else str(c))
        for c in df.columns
    ]
    return df


def find_col(df: pd.DataFrame, *needles: str) -> str:
    for c in df.columns:
        cl = c.lower()
        if all(n.lower() in cl for n in needles):
            return c
    raise KeyError(f"No column matching {needles}. Available: {list(df.columns)[:40]}")


def main() -> None:
    raw_fn   = Path(sys.argv[1])
    final_fn = Path(sys.argv[2])
    out_fn   = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(
        "/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/"
        "supplementary_table_S1_attrition.csv")

    print(f"Reading {raw_fn} ...")
    df = pd.read_csv(raw_fn, header=[0, 1], low_memory=False)
    df = flatten_columns(df)

    allele_col = find_col(df, "MHC Restriction", "Name")
    resp_col   = find_col(df, "Response measured")
    units_col  = find_col(df, "Units")
    qty_col    = find_col(df, "Quantitative measurement")
    epi_col    = find_col(df, "Epitope", "Name")

    # --- Stage counts (each stage applied sequentially) ---
    n_raw = len(df)

    ab = df[df[allele_col].astype(str).str.startswith(("HLA-A", "HLA-B"))]
    n_ab = len(ab)

    vu = ab.dropna(subset=[qty_col, units_col])
    n_vu = len(vu)

    resp_std = vu[resp_col].replace({
        "dissociation constant KD (~EC50)": KD_STD,
        "dissociation constant KD":         KD_STD,
        "dissociation constant (~IC50)":    KD_STD,
    })
    vr = vu[resp_std.isin([KD_STD, IC50_STD])]
    n_vr = len(vr)

    final = pd.read_csv(final_fn, low_memory=False)
    n_final_rows  = len(final)
    n_final_pairs = (final[["allele", "peptide"]].drop_duplicates().shape[0]
                     if {"allele", "peptide"}.issubset(final.columns) else float("nan"))

    # --- Build the 4-column funnel table (matches manuscript Supplementary Table S2) ---
    rows = [
        ("0. Raw IEDB MHC ligand records",
         "IEDB MHC ligand bulk download",                       n_raw,   None),
        ("1. HLA-A / HLA-B allele restriction",
         "MHC allele name begins with HLA-A or HLA-B",          n_ab,    n_raw - n_ab),
        ("2. Quantitative value + assay units present",
         "Both a quantitative measurement value and assay units", n_vu,  n_ab - n_vu),
        ("3. Retained assay response (KD or IC50)",
         "Assay response is KD or IC50 (variant KD labels normalized)", n_vr, n_vu - n_vr),
        ("4-6. Deduplication, flagging, and residue filters (net)",
         "Duplicate resolution, flagged-record removal, non-canonical (+) exclusion",
         n_final_rows, n_vr - n_final_rows),
    ]

    out = pd.DataFrame(
        rows,
        columns=["Stage", "Filter applied", "Records remaining", "Records removed"],
    )

    out.to_csv(out_fn, index=False)

    print("\nAttrition funnel:")
    print(out.to_string(index=False))
    print(f"\nFinal released rows (measurements): {n_final_rows:,}")
    print(f"Final released unique pairs:        {n_final_pairs:,}")
    print(f"Wrote {out_fn}")

    if n_final_rows != 49488:
        print(f"\nWARNING: final rows {n_final_rows:,} != 49,488 expected.")
    else:
        print("\nFinal row count matches 49,488.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    main()
