"""
Restore the release columns that the v1+v2 metadata reconciliation dropped.

release_v2_final/metadata_curated.csv came out of the merge with only the
curation columns (16 of them). The released v1 metadata.csv carries ten more that Data
Records documents and that users' code reads:

    flagged  has_structures  num_pdbs  pdb_dir
    I_sc_best  I_sc_mean  reweighted_sc_best  reweighted_sc_mean
    total_score_best  total_score_mean
    rosetta_best_score  rosetta_mean_score

This adds them back, deriving each from a file on disk rather than copying v1
values forward:

  * the six score columns come from the per-pair score summaries produced by
    parse_scorefiles.py (scores_out/ for v1, scores_out_v2/ for v2), so they
    describe the silent files actually in release_v2_final/structures/ -- which
    matters, because the v1 silents were regenerated after the two conversion
    bugs were fixed.
  * num_pdbs / has_structures come from the same summaries' n_decoys, i.e. the
    decoys really present, not an assumed 25. (Nine v1 pairs have 22-24.)
  * pdb_dir is constructed, matching v1's "structures/{allele}/{peptide}" form.
  * rosetta_best_score / rosetta_mean_score are legacy duplicates of
    total_score_best / total_score_mean -- verified against released v1, where
    rosetta_best_score == total_score_best for every row. Carried for backward
    compatibility with code written against v1.
  * flagged is carried from released v1 by measurement identity. Only 15 rows
    are True and all 15 survive the merge.

pep_sc is written only with --with-pep-sc. parse_scorefiles.py computes it and
Table S6 uses it, but Data Records says "three scores", so shipping a fourth is
a documentation decision, not a default.

KNOWN LIMITATION, reported on every run: the v2 curation run's flag
determination was not preserved anywhere on disk (no v2-era file carries a
flagged column), so v2 rows are written False. If the v2 pipeline flagged
records, that is not represented here.

Usage:
    python3 add_release_columns.py --out <path>          # writes a new file
    python3 add_release_columns.py --self-test
"""
import argparse
from pathlib import Path

import pandas as pd

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT, MHC_DB_ROOT  # noqa: E402

BASE = DATA_ROOT
# The curation output, before score and provenance columns are attached. Named
# for what it is: the release file is metadata.csv, written by this script.
MERGED = BASE / "release_v2_final" / "metadata_curated.csv"
V1_RELEASED = BASE / "huggingface" / "metadata.csv"
# ORDER MATTERS. load_scores() dedups on (allele_dir, peptide) with keep="last",
# so a later file overrides an earlier one for the same pair. The 1,162
# re-docked pairs must therefore come LAST: their entries in scores_out/ are the
# defective originals, whose I_sc was computed for the wrong chain (median
# -11.60 REU against -68.57 after the re-dock) or for a truncated peptide.
SCORE_SUMMARIES = [
    BASE / "IEDB_validation" / "scores_out" / "score_summary.csv",
    BASE / "IEDB_validation" / "scores_out_v2" / "score_summary.csv",
    BASE / "IEDB_validation" / "scores_out_1162" / "score_summary.csv",
    BASE / "IEDB_validation" / "scores_out_183" / "score_summary.csv",
]

MHC_DB = MHC_DB_ROOT / "database.info"

SCORE_METRICS = ("I_sc", "reweighted_sc", "total_score")
PEP_SC = "pep_sc"
# v1 column order, so the restored file is a drop-in for code that reads v1.
V1_COLUMN_ORDER = [
    "allele_iedb", "allele", "allele_compact", "peptide", "peptide_length",
    "measurement_type", "measurement_value", "measurement_units",
    "assay_method", "assay_response", "pubmed_id", "parent_protein",
    "protein_accession", "source_organism", "source_version",
    "flagged", "self_templated", "has_structures", "num_pdbs",
    "I_sc_best", "I_sc_mean", "reweighted_sc_best", "reweighted_sc_mean",
    "total_score_best", "total_score_mean",
    "rosetta_best_score", "rosetta_mean_score", "pdb_dir",
]

# Measurement identity used to carry `flagged` across the merge. allele_iedb is
# the filesystem-safe spelling in both files, so it joins directly.
FLAG_KEY = ["allele_iedb", "peptide", "measurement_type", "measurement_value"]


def load_scores(with_pep_sc):
    metrics = list(SCORE_METRICS) + ([PEP_SC] if with_pep_sc else [])
    cols = ["allele_dir", "peptide", "n_decoys"] + [
        f"{m}_{a}" for m in metrics for a in ("best", "mean")]
    frames = []
    for p in SCORE_SUMMARIES:
        s = pd.read_csv(p)
        missing = [c for c in cols if c not in s.columns]
        if missing:
            raise SystemExit(f"{p} is missing columns {missing}; rerun parse_scorefiles.py")
        frames.append(s[cols])
    scores = pd.concat(frames, ignore_index=True)
    n_before = len(scores)
    scores = scores.drop_duplicates(subset=["allele_dir", "peptide"], keep="last")
    if len(scores) != n_before:
        print(f"  note: {n_before - len(scores):,} duplicate (allele, peptide) score rows "
              f"collapsed (last wins)")
    return scores, metrics


def mark_self_templated(md, db_peptides):
    """Flag pairs that were threaded onto their own crystal structure.

    `get_peptide_template` ranks same-length peptides from every allele of the
    same gene by BLOSUM62 similarity to the target, so an identical peptide wins;
    self-exclusion runs only under --ignore_epitope_match, which the production
    runs did not pass. Membership of the target peptide in database.info is
    exactly the self-exclusion criterion (HLA_db.py matches on
    Epitope_Description == query_sequence, across all alleles), so it identifies
    precisely the affected pairs.

    These structures are the MOST accurate in the release -- they were built from
    real coordinates. The column exists so users can exclude them when measuring
    structure-prediction accuracy, where including them inflates the result."""
    return md["peptide"].astype(str).isin(db_peptides)


def add_columns(md, scores, flags, metrics, db_peptides=None):
    """Attach the release columns. Returns (dataframe, stats dict)."""
    # assay_pdb_id is dropped, not carried: it is empty for every row in both
    # the earlier upload (0 of 49,488) and here (0 of 118,985). Data Records
    # described it as "the experimental PDB ID if reported by IEDB", so shipping
    # it would document a column that contains nothing. The crystal structures
    # used for validation are matched through the local template database
    # instead, and are listed in the validation table.
    md = md.drop(columns=["assay_pdb_id"], errors="ignore")

    # pubmed_id arrives as a float because the column holds NaNs, so it
    # serializes as "22508927.0". Cast to a nullable integer so the release
    # carries "22508927" and a reader can use it as an identifier directly.
    if "pubmed_id" in md.columns:
        md["pubmed_id"] = pd.to_numeric(md["pubmed_id"], errors="coerce").astype("Int64")
    md = md.merge(scores, how="left",
                  left_on=["allele_compact", "peptide"],
                  right_on=["allele_dir", "peptide"]).drop(columns=["allele_dir"])

    md["has_structures"] = md["n_decoys"].notna()
    md["num_pdbs"] = md["n_decoys"].astype("Int64")
    md = md.drop(columns=["n_decoys"])
    md["pdb_dir"] = "structures/" + md["allele_compact"].astype(str) + "/" + md["peptide"].astype(str)

    # legacy aliases; asserted equal to total_score in released v1
    md["rosetta_best_score"] = md["total_score_best"]
    md["rosetta_mean_score"] = md["total_score_mean"]

    if db_peptides is not None:
        md["self_templated"] = mark_self_templated(md, db_peptides)

    md = md.merge(flags, how="left", on=FLAG_KEY)
    n_flagged_matched = int(md["flagged"].notna().sum())
    # the merge leaves an object column of True/NaN; fill then cast explicitly
    # rather than letting fillna silently downcast (deprecated in pandas 2.2)
    md["flagged"] = md["flagged"].notna() & (md["flagged"] == True)  # noqa: E712

    order = [c for c in V1_COLUMN_ORDER if c in md.columns]
    extra = [c for c in md.columns if c not in order]
    md = md[order + extra]

    stats = {
        "rows": len(md),
        "no_structures": int((~md["has_structures"]).sum()),
        "flagged_true": int(md["flagged"].sum()),
        "flagged_carried": n_flagged_matched,
        "missing_score": int(md[f"{metrics[0]}_best"].isna().sum()),
        "num_pdbs_lt_25": int((md["num_pdbs"] < 25).sum()),
        "self_templated": int(md["self_templated"].sum()) if "self_templated" in md else None,
    }
    return md, stats


def self_test():
    """Known-answer test: two pairs, one with scores and 24 decoys, one with no
    structures at all, and one flagged measurement. Every derived column has an
    analytically known value."""
    md = pd.DataFrame({
        "allele_iedb": ["HLA-A_02_01", "HLA-A_02_01", "HLA-B_07_02"],
        "allele_compact": ["A0201", "A0201", "B0702"],
        "peptide": ["AAAAAAAAA", "AAAAAAAAA", "CCCCCCCCC"],
        "measurement_type": ["IC50", "Kd", "IC50"],
        "measurement_value": [10.0, 20.0, 30.0],
        "source_version": ["v1", "v1", "v2"],
    })
    scores = pd.DataFrame({
        "allele_dir": ["A0201"], "peptide": ["AAAAAAAAA"], "n_decoys": [24],
        "I_sc_best": [-70.0], "I_sc_mean": [-65.0],
        "reweighted_sc_best": [-500.0], "reweighted_sc_mean": [-490.0],
        "total_score_best": [-580.0], "total_score_mean": [-570.0],
    })
    # only the Kd measurement of the first pair is flagged
    flags = pd.DataFrame({
        "allele_iedb": ["HLA-A_02_01"], "peptide": ["AAAAAAAAA"],
        "measurement_type": ["Kd"], "measurement_value": [20.0], "flagged": [True],
    })
    # AAAAAAAAA has a crystal in the template DB, CCCCCCCCC does not
    out, stats = add_columns(md, scores, flags, list(SCORE_METRICS),
                             db_peptides={"AAAAAAAAA", "UNRELATEDPEP"})
    assert list(out["self_templated"]) == [True, True, False], list(out["self_templated"])
    assert stats["self_templated"] == 2, stats

    # omitting db_peptides must leave the column off entirely, not write False
    out_no, stats_no = add_columns(md, scores, flags, list(SCORE_METRICS))
    assert "self_templated" not in out_no.columns
    assert stats_no["self_templated"] is None


    assert list(out["has_structures"]) == [True, True, False], list(out["has_structures"])
    # Int64 keeps the missing entry as pd.NA, not None -- compare the present
    # values and the missing mask separately.
    assert list(out["num_pdbs"][:2]) == [24, 24], list(out["num_pdbs"])
    assert list(out["num_pdbs"].isna()) == [False, False, True], list(out["num_pdbs"])
    assert list(out["pdb_dir"]) == ["structures/A0201/AAAAAAAAA",
                                    "structures/A0201/AAAAAAAAA",
                                    "structures/B0702/CCCCCCCCC"]
    # the flag must land on exactly the one measurement it belongs to, not on
    # the other measurement of the same pair
    assert list(out["flagged"]) == [False, True, False], list(out["flagged"])
    assert out["rosetta_best_score"].iloc[0] == -580.0
    assert out["rosetta_mean_score"].iloc[0] == -570.0
    assert pd.isna(out["I_sc_best"].iloc[2]), "a pair with no scores must stay NaN, not 0"
    assert stats["no_structures"] == 1 and stats["flagged_true"] == 1
    assert stats["num_pdbs_lt_25"] == 2, stats

    # Column order must start with the v1 order. The `or True` that used to
    # terminate this line made the whole assert vacuous, so it never guarded
    # anything: reordering the released columns would have gone unnoticed.
    expected = [c for c in V1_COLUMN_ORDER if c in out.columns][:5]
    assert list(out.columns)[:5] == expected, (
        f"released column order changed: {list(out.columns)[:5]} != {expected}")
    print("Self-test PASSED: derived columns match known values, the flag lands on "
          "the correct measurement of a two-measurement pair, and a pair without "
          "scores stays NaN rather than zero.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", help="output CSV. Never overwrites the input.")
    ap.add_argument("--with-pep-sc", action="store_true",
                    help="also write pep_sc_best/pep_sc_mean (documentation decision)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not args.out:
        raise SystemExit("--out is required")
    out = Path(args.out)
    if out.resolve() == MERGED.resolve():
        raise SystemExit("refusing to overwrite the input metadata; choose another --out")

    self_test()
    print()

    md = pd.read_csv(MERGED, low_memory=False)
    print(f"merged metadata: {len(md):,} rows, {len(md.columns)} columns")

    scores, metrics = load_scores(args.with_pep_sc)
    print(f"score summaries: {len(scores):,} pairs")

    v1 = pd.read_csv(V1_RELEASED, low_memory=False)
    # sanity-check the legacy alias claim before relying on it
    same = (v1["rosetta_best_score"].round(4) == v1["total_score_best"].round(4)).all()
    if not same:
        raise SystemExit("rosetta_best_score != total_score_best in released v1; "
                         "the legacy-alias assumption is wrong, fix add_columns()")
    flags = v1.loc[v1["flagged"] == True, FLAG_KEY + ["flagged"]].drop_duplicates(  # noqa: E712
        subset=FLAG_KEY)
    print(f"flagged rows in released v1: {len(flags)}")

    db_peptides = set(pd.read_csv(MHC_DB)["Epitope_Description"].astype(str))
    print(f"template database: {len(db_peptides):,} distinct peptides")

    result, stats = add_columns(md, scores, flags, metrics, db_peptides)
    result.to_csv(out, index=False)

    print(f"\nwrote {out}")
    print(f"  rows                  {stats['rows']:,}")
    print(f"  columns               {len(result.columns)}")
    print(f"  pairs missing a score {stats['missing_score']:,}")
    print(f"  rows w/o structures   {stats['no_structures']:,}")
    print(f"  num_pdbs < 25         {stats['num_pdbs_lt_25']:,}")
    if stats.get("self_templated") is not None:
        n = stats["self_templated"]
        pairs = result.drop_duplicates(subset=["allele_compact", "peptide"])
        print(f"  self_templated rows   {n:,} "
              f"({int(pairs['self_templated'].sum()):,} distinct pairs, "
              f"{100*pairs['self_templated'].mean():.3f}% of pairs)")
    print(f"  flagged == True       {stats['flagged_true']} "
          f"(carried {stats['flagged_carried']} from released v1)")

    print("\nLIMITATION: the v2 curation run's flag determination is not preserved on "
          "disk, so v2 rows are written flagged=False. If the v2 pipeline flagged "
          "records, they are not marked here.")
    if "assay_pdb_id" in result.columns and result["assay_pdb_id"].notna().sum() == 0:
        print("NOTE: assay_pdb_id is empty for every row (it is also empty in released "
              "v1). Either populate it from the raw IEDB file or drop the column and "
              "its sentence in Data Records.")


if __name__ == "__main__":
    main()
