"""
Re-run the Validation 1 crystal-structure match on the merged release.

The 52 pairs quoted in the submitted manuscript were matched against the local
MHC template database from v1's 37-allele metadata (01_structural_validation.ipynb,
cell 6). The merged release has 95 alleles and 112,378 pairs, including HLA-C,
so the match has to be redone -- the count is not carried forward.

This reproduces the notebook's matching logic exactly, as a script so the number
is re-runnable (CLAUDE.md: never invent a number):

  * inner merge of unflagged metadata rows against database.info on
    (allele, peptide), where database.info's MHC_Allele -> allele and
    Epitope_Description -> peptide;
  * where several PDB entries exist for one pair, keep the highest-resolution
    one (sort by resolution, drop duplicates);
  * report per-locus and per-source-version breakdowns so it is clear which
    matches are new.

It does NOT compute RMSD. RMSD needs decoy PDBs, and the v2 decoys exist only
inside silent files on this host, so extraction is a separate (and much more
expensive) step. This script's job is to establish how many pairs qualify.

Usage:
    python3 crystal_match.py --out-dir <dir>
    python3 crystal_match.py --out-dir <dir> --metadata <path>   # e.g. v1, to compare
    python3 crystal_match.py --self-test
"""
import argparse
from pathlib import Path

import pandas as pd

BASE = Path("<HOME>/main_project/data/IEDB_data_clean")
# the restored-column merged metadata; add_release_columns.py output, which has
# `flagged` back, so the notebook's unflagged filter can be applied
METADATA = BASE / "release_v2_final" / "metadata_full.csv"
V1_METADATA = BASE / "huggingface" / "metadata.csv"
MHC_DB = Path("<HOME>/Data/MHC_database/database.info")

DB_RENAME = {
    "MHC_Allele": "allele",
    "Epitope_Description": "peptide",
    "PDB_ID": "matched_pdb_id",
    "MHC_PDB_Chain1": "mhc_chain_id",
    "Antigen_PDB_Chain(s)": "peptide_chain_id",
    "Resolution_(Angstrom)": "resolution_angstrom",
}


def match(md, db):
    """Inner-join metadata pairs to template-database entries on (allele,
    peptide), keeping the highest-resolution PDB per pair."""
    if "flagged" in md.columns:
        md = md[md["flagged"] == False]  # noqa: E712
    keep = [c for c in ["allele", "allele_compact", "peptide", "source_version"]
            if c in md.columns]
    pairs = md[keep].drop_duplicates(subset=["allele", "peptide"])

    db = db.rename(columns=DB_RENAME)
    db = db[[c for c in DB_RENAME.values() if c in db.columns]]

    m = pairs.merge(db, on=["allele", "peptide"], how="inner")
    # several crystal structures can exist for one peptide-allele pair; the
    # notebook keeps the best-resolution one, so the RMSD reference is unique
    m = (m.sort_values("resolution_angstrom")
          .drop_duplicates(subset=["allele", "peptide"], keep="first")
          .reset_index(drop=True))
    return m


def self_test():
    """Known answer: three metadata pairs and a template DB where one pair has
    two structures (1.0 and 2.5 A), one pair has one, one pair has none, and a
    flagged row must be excluded even though it would otherwise match."""
    md = pd.DataFrame({
        "allele":         ["A*02:01", "B*07:02", "A*01:01", "C*07:01"],
        "allele_compact": ["A0201",   "B0702",   "A0101",   "C0701"],
        "peptide":        ["GILGFVFTL", "RPPIFIRRL", "NOPEPTIDE", "FLAGGEDONE"],
        "source_version": ["v1", "v2", "v1", "v2"],
        "flagged":        [False, False, False, True],
    })
    db = pd.DataFrame({
        "MHC_Allele":            ["A*02:01", "A*02:01", "B*07:02", "C*07:01"],
        "Epitope_Description":   ["GILGFVFTL", "GILGFVFTL", "RPPIFIRRL", "FLAGGEDONE"],
        "PDB_ID":                ["LOWRES", "BESTRES", "ONLYONE", "SHOULDNOTAPPEAR"],
        "MHC_PDB_Chain1":        ["A", "A", "A", "A"],
        "Antigen_PDB_Chain(s)":  ["C", "C", "C", "C"],
        "Resolution_(Angstrom)": [2.5, 1.0, 1.9, 1.5],
    })
    m = match(md, db)
    assert len(m) == 2, f"expected 2 matches, got {len(m)}:\n{m}"
    got = dict(zip(m["peptide"], m["matched_pdb_id"]))
    assert got["GILGFVFTL"] == "BESTRES", f"resolution tie-break failed: {got}"
    assert got["RPPIFIRRL"] == "ONLYONE", got
    assert "FLAGGEDONE" not in got, "a flagged row was matched; the filter is broken"
    assert "NOPEPTIDE" not in got
    print("Self-test PASSED: highest-resolution structure wins per pair, flagged "
          "rows are excluded, and pairs with no template-DB entry are dropped.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir")
    ap.add_argument("--metadata", default=str(METADATA))
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not args.out_dir:
        raise SystemExit("--out-dir is required")
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    self_test()
    print()

    md = pd.read_csv(args.metadata, low_memory=False)
    db = pd.read_csv(MHC_DB)
    print(f"metadata: {args.metadata}")
    n_unflagged = int((md["flagged"] == False).sum()) if "flagged" in md else len(md)  # noqa: E712
    n_pairs = md.drop_duplicates(subset=["allele", "peptide"]).shape[0]
    print(f"  {len(md):,} rows, {n_unflagged:,} unflagged, {n_pairs:,} unique pairs, "
          f"{md['allele'].nunique()} alleles")
    print(f"template database: {len(db):,} structures "
          f"({db['MHC_Allele'].str[0].value_counts().to_dict()})")

    m = match(md, db)
    m.to_csv(out / "crystal_matched_pairs.csv", index=False)

    print(f"\nPairs with an experimental match: {len(m)}")
    print(f"  resolution: median {m['resolution_angstrom'].median():.2f} A, "
          f"range {m['resolution_angstrom'].min():.2f}-{m['resolution_angstrom'].max():.2f} A")
    print(f"  by locus:   {m['allele'].str[0].value_counts().to_dict()}")
    if "source_version" in m:
        print(f"  by version: {m['source_version'].value_counts().to_dict()}")
    print(f"  distinct alleles: {m['allele'].nunique()}, distinct PDB IDs: "
          f"{m['matched_pdb_id'].nunique()}")
    print(f"\nwrote {out / 'crystal_matched_pairs.csv'}")


if __name__ == "__main__":
    main()
