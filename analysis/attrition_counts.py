"""Per-stage record attrition for the IEDB pHLA curation pipeline.

Reproduces the filtering in IEDBTestPipeline.process_peplist against the raw
IEDB bulk download, logging record counts after each stage. Produces the
Supplementary Table S1 funnel (Rocco comment 5).

This mirrors the pipeline logic; it does not re-run modeling. The final row
should match the released dataset (118,985 measurements / 112,561 pairs). Any
gap indicates a filter applied outside this script.

Usage:
    python attrition_counts.py <mhc_ligand_full.csv> <final_metadata.csv> [out.csv]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT  # noqa: E402


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
    raw_fn     = Path(sys.argv[1])
    curated_fn = Path(sys.argv[2])
    final_fn   = Path(sys.argv[3])
    out_fn   = Path(sys.argv[4]) if len(sys.argv) > 4 else Path(
        DATA_ROOT / "IEDB_validation" / "supplementary_table_S1_attrition.csv")

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

    # All three classical class I loci. This filter was HLA-A/B only, which was
    # right for the 37-allele submission but wrong once HLA-C was added: the
    # last row of the table reads its totals from the three-locus release, so
    # the stage differences were being taken across incompatible pools.
    ab = df[df[allele_col].astype(str).str.startswith(("HLA-A", "HLA-B", "HLA-C"))]
    n_ab = len(ab)

    vu = ab.dropna(subset=[qty_col, units_col])
    n_vu = len(vu)

    # "dissociation constant (~IC50)" is the typo carried from the pipeline
    # (IEDBTestPipeline.py:136, missing "KD"). That string never occurs in the
    # raw export, so the third variant was never normalized and its ~42,660
    # records fell out here while the release kept them -- which made stage 3
    # smaller than the final row and the funnel report a negative removal.
    resp_std = vu[resp_col].replace({
        "dissociation constant KD (~EC50)": KD_STD,
        "dissociation constant KD":         KD_STD,
        "dissociation constant KD (~IC50)": KD_STD,
    })
    vr = vu[resp_std.isin([KD_STD, IC50_STD])]
    n_vr = len(vr)

    # The released file is curation AND structure generation: pairs that were
    # curated but never produced an ensemble are simply absent from it, so
    # ending the funnel on the release silently folds that loss into "curation".
    # `curated_fn` is the curation endpoint, and the two are reported as
    # separate stages.
    curated = pd.read_csv(curated_fn, low_memory=False)
    n_cur_rows  = len(curated)
    n_cur_pairs = curated[["allele", "peptide"]].drop_duplicates().shape[0]

    final = pd.read_csv(final_fn, low_memory=False)
    n_final_rows  = len(final)
    n_final_pairs = (final[["allele", "peptide"]].drop_duplicates().shape[0]
                     if {"allele", "peptide"}.issubset(final.columns) else float("nan"))

    # Why each unreleased pair was dropped, so the last two stages can be named
    # by cause rather than lumped as "structure generation".
    orphans_fn = Path(final_fn).parent.parent / "release_v2_final" / "orphan_pairs_no_structure.csv"
    if not orphans_fn.exists():
        orphans_fn = Path(sys.argv[5]) if len(sys.argv) > 5 else orphans_fn
    orph = pd.read_csv(orphans_fn)
    rel_pair_set = set(map(tuple, final[["allele", "peptide"]].drop_duplicates().values))
    unreleased = orph[[(a, p) not in rel_pair_set
                       for a, p in zip(orph["allele"], orph["peptide"])]]
    n_noncanon = int((unreleased["exclusion_reason"] == "non_canonical_or_PTM").sum())
    n_no_template = int(len(unreleased) - n_noncanon)

    # --- Build the 4-column funnel table (matches manuscript Supplementary Table S2) ---
    rows = [
        ("0. Raw IEDB MHC ligand records",
         "IEDB MHC ligand bulk download",                       n_raw,   None),
        ("1. HLA-A / HLA-B / HLA-C allele restriction",
         "MHC allele name begins with HLA-A, HLA-B or HLA-C",          n_ab,    n_raw - n_ab),
        ("2. Quantitative value + assay units present",
         "Both a quantitative measurement value and assay units", n_vu,  n_ab - n_vu),
        ("3. Retained assay response (KD or IC50)",
         "Assay response is KD or IC50 (variant KD labels normalized)", n_vr, n_vu - n_vr),
        ("4. Deduplication and removal of curation-flagged records",
         "Duplicate resolution and removal of records flagged for manual review",
         n_cur_rows, n_vr - n_cur_rows),
        # The sequence filter runs at structure generation, not during curation:
        # the curation table still holds 233 "+" peptides and sequences up to 30
        # residues. Both remaining stages are named by cause, taken from the
        # exclusion_reason column of orphan_pairs_no_structure.csv.
        ("5. Exclusion of non-canonical or modified peptides",
         "Peptides carrying the IEDB \"+\" notation for modified residues",
         n_cur_rows - n_noncanon, n_noncanon),
        ("6. Exclusion of peptides with no length-matched template",
         "No template of the peptide's length for that allele in the MHC "
         "template database",
         n_final_rows, n_no_template),
    ]

    out = pd.DataFrame(
        rows,
        columns=["Stage", "Filter applied", "Records remaining", "Records removed"],
    )

    out.to_csv(out_fn, index=False)

    cur_pairs = set(map(tuple, curated[["allele", "peptide"]].drop_duplicates().values))
    rel_pairs = set(map(tuple, final[["allele", "peptide"]].drop_duplicates().values))
    orphans = cur_pairs - rel_pairs
    stray = rel_pairs - cur_pairs
    assert not stray, (
        f"{len(stray):,} released pairs are absent from the curation table, so "
        "these are not two stages of one funnel")
    print(f"\ncuration endpoint : {n_cur_rows:,} rows / {n_cur_pairs:,} pairs")
    print(f"released          : {n_final_rows:,} rows / {n_final_pairs:,} pairs")
    print(f"structureless     : {len(orphans):,} pairs")

    print("\nAttrition funnel:")
    print(out.to_string(index=False))
    print(f"\nFinal released rows (measurements): {n_final_rows:,}")
    print(f"Final released unique pairs:        {n_final_pairs:,}")
    print(f"Wrote {out_fn}")

    # Check the funnel is monotonic instead of against a typed-in total: a
    # hardcoded expectation just goes stale, whereas a stage that removes a
    # negative number of records means a filter here does not match the
    # pipeline and the table would be nonsense.
    removed = out["Records removed"].dropna()
    if (removed < 0).any():
        bad = out.loc[out["Records removed"] < 0, "Stage"].tolist()
        print(f"\nWARNING: these stages removed a negative number of records, so "
              f"a filter here does not match the pipeline: {bad}")
    else:
        print(f"\nFunnel is monotonic; final row {n_final_rows:,} measurements / "
              f"{n_final_pairs:,} pairs.")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(
            "usage: attrition_counts.py <raw_export> <curation_endpoint> "
            "<released_metadata> [out_csv]")
    main()
