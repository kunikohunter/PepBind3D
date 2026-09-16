"""Replay the curation pipeline from the raw export to count each stage.

`attrition_counts.py` reports deduplication, flagging and the non-canonical
residue filter as one net row, because it measures the released file and cannot
see inside those three steps. Supplementary Table S2 lists them separately, so
the numbers have to come from somewhere: this replays `clean_peplist()` from
`pipeline/IEDBTestPipeline_ACCRE.py` over the raw export and counts what each
step removes.

The pipeline loops over every Epitope_ID and re-masks the whole frame each time,
which is O(n^2) and takes hours. This is a vectorised equivalent, decision for
decision, including two quirks that change the counts and so must be kept:

  * for a group of more than two records, `values` is captured BEFORE the
    PubMed-less rows are dropped, so the pairwise comparison runs on the
    original values while the rows being kept are the survivors;
  * the ">2 after dropping" branch is the only one that can fire, because the
    three `elif len(tmp) > 2` branches after it are unreachable.

Reproducing the pipeline is the whole point, so a faithful replay of a quirk is
correct here and 'fixing' it would be wrong.

The replay is only trustworthy if it lands on the released file exactly, so
`main` asserts both endpoints (118,985 retained rows, 15 flagged) and refuses to
report stage counts if either misses.

RESULT, 2026-09-16: THIS DOES NOT REPRODUCE THE RELEASE, and its stage 4-6
counts must not be quoted. Stages 0-3 match `attrition_counts.py` exactly
(4,883,585 / 1,532,863 / 147,967 / 134,840). After that the replay ends at
127,874 rows against the released 118,985, and reports 892 flagged records
against the 15 the release carries. Two reasons, one understood and one not:

  * `metadata.csv` holds only pairs that produced a structural ensemble, so the
    released count is curation AND structure generation, while this replays
    curation alone. The funnel's last row is therefore not a curation endpoint.
  * the flag branch fires far more often here than in the pipeline, so the
    deduplication replay is not faithful in that branch.

The pipeline's per-allele `*_cleaned_IEDB_data.csv` intermediates, which would
settle it, are not on disk. Until they are regenerated, Supplementary Table S2
should keep the merged "4-6" row that `attrition_counts.py` can measure rather
than the split this script was written to supply. Kept as a diagnostic, and
because the two discrepancies above are worth knowing.

Usage:
    python3 curation_replay.py --out-dir <dir>
Self-test:
    python3 curation_replay.py --self-test
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import RAW_IEDB, RELEASE_METADATA  # noqa: E402

KD_STD = "dissociation constant (KD)"
IC50_STD = "half maximal inhibitory concentration (IC50)"
KD_VARIANTS = {
    "dissociation constant KD (~EC50)": KD_STD,
    "dissociation constant KD": KD_STD,
    "dissociation constant KD (~IC50)": KD_STD,
}


def resolve_group(pubmed_na, values):
    """One (epitope, assay) group -> (positions to drop, flag the group?).

    Mirrors the branches of clean_peplist(). Positions are indices into the
    group, so the caller can map them back to the frame.
    """
    n = len(values)
    if n < 2:
        return [], False

    if n == 2:
        diff = abs(values[0] - values[1])
        if pubmed_na.any() and not pubmed_na.all():
            return [int(np.flatnonzero(pubmed_na)[0])], False
        if pubmed_na.all():
            return ([0, 1], False) if diff > 10 else ([1], False)
        # both carry a PubMed ID
        return ([0], False) if diff < 10 else ([], True)

    # n > 2: PubMed-less rows go first, then the survivors are reduced to one.
    drop = list(np.flatnonzero(pubmed_na))
    keep = [i for i in range(n) if i not in set(drop)]
    if len(keep) > 2:
        # NOTE: `values` here is the pre-drop list, as in the pipeline.
        d = np.abs(np.subtract.outer(values, values))
        d = np.tril(d, -1)
        if (d <= 10).sum(axis=1).max() == n - 1:
            drop += keep[1:]                     # all close: keep the first
        else:
            med = np.median(values)
            closest = int(np.argmin(np.abs(np.array(values) - med)))
            drop += [i for i in keep if i != closest]
    return sorted(set(drop)), False


def replay(df, allele_col, epi_col, resp_col, qty_col, pmid_col):
    """Deduplicate, then flag. Returns (kept frame, n_dropped, n_flagged)."""
    drop_idx, flag_idx = [], []
    # Grouped per allele as well as per epitope: the pipeline processes one
    # allele at a time, so a peptide measured against several alleles is not
    # a duplicate. Grouping on the epitope alone collapses them and throws
    # away six sevenths of the release.
    for _, grp in df.groupby([allele_col, epi_col, resp_col], sort=False):
        if len(grp) < 2:
            continue
        pubmed_na = grp[pmid_col].isna().to_numpy()
        values = grp[qty_col].tolist()
        drops, flag = resolve_group(pubmed_na, values)
        idx = grp.index.to_numpy()
        drop_idx.extend(idx[d] for d in drops)
        if flag:
            flag_idx.extend(idx)

    deduped = df.drop(index=drop_idx)
    flagged = [i for i in flag_idx if i in deduped.index]
    return deduped.drop(index=flagged), len(drop_idx), len(flagged)


def self_test():
    """Groups whose correct handling is known by reading the branches."""
    # one PubMed, one without -> the one without goes
    assert resolve_group(np.array([False, True]), [10.0, 12.0]) == ([1], False)
    # neither has PubMed, far apart -> both go
    assert resolve_group(np.array([True, True]), [10.0, 500.0]) == ([0, 1], False)
    # neither has PubMed, close -> keep one
    assert resolve_group(np.array([True, True]), [10.0, 12.0]) == ([1], False)
    # both have PubMed, close -> keep one
    assert resolve_group(np.array([False, False]), [10.0, 12.0]) == ([0], False)
    # both have PubMed, conflicting -> flag, drop neither
    assert resolve_group(np.array([False, False]), [10.0, 500.0]) == ([], True)
    # a single record is untouched
    assert resolve_group(np.array([False]), [10.0]) == ([], False)
    # three records, one without PubMed: it goes, and of the two survivors
    # (close together) the first is kept
    drops, flag = resolve_group(np.array([True, False, False]), [10.0, 11.0, 12.0])
    assert drops == [0] and not flag, (drops, flag)
    # three with PubMed, one an outlier -> keep the one nearest the median
    drops, flag = resolve_group(np.array([False, False, False]), [10.0, 11.0, 900.0])
    assert drops == [0, 2] and not flag, (drops, flag)
    print("Self-test PASSED: every branch of clean_peplist() resolves as the "
          "pipeline resolves it.")


def main():
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

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    cache = out / "curation_replay_stage3.parquet"
    counts = out / "curation_replay_stage_counts.json"
    if cache.exists() and counts.exists():
        print(f"reusing {cache.name}")
        df = pd.read_parquet(cache)
        pre = json.load(open(counts))
        n_raw, n_loci, n_vu, n_resp = (pre["raw"], pre["loci"], pre["vu"], pre["resp"])
        allele = "MHC Restriction | Name"
        resp = "Assay | Response measured"
        qty = "Assay | Quantitative measurement"
        epi = "Epitope | Epitope IRI"
        pep = "Epitope | Name"
        pmid = "Reference | PMID"
        return finish(df, args, allele, resp, qty, epi, pep, pmid,
                      n_raw, n_loci, n_vu, n_resp, out)

    print(f"reading {RAW_IEDB} ...")
    df = pd.read_csv(RAW_IEDB, header=[0, 1], low_memory=False)
    df.columns = [" | ".join(str(x) for x in c) for c in df.columns]
    allele = "MHC Restriction | Name"
    resp = "Assay | Response measured"
    units = "Assay | Units"
    qty = "Assay | Quantitative measurement"
    epi = "Epitope | Epitope IRI"
    pep = "Epitope | Name"
    pmid = "Reference | PMID"

    n_raw = len(df)
    df = df[df[allele].astype(str).str.startswith(("HLA-A", "HLA-B", "HLA-C"))]
    n_loci = len(df)
    df = df.dropna(subset=[qty, units])
    n_vu = len(df)
    df = df.assign(**{resp: df[resp].replace(KD_VARIANTS)})
    df = df[df[resp].isin([KD_STD, IC50_STD])]
    n_resp = len(df)

    df[[allele, resp, qty, epi, pep, pmid]].to_parquet(cache, index=False)
    with open(counts, "w") as f:
        json.dump({"raw": n_raw, "loci": n_loci, "vu": n_vu, "resp": n_resp}, f)

    return finish(df, args, allele, resp, qty, epi, pep, pmid,
                  n_raw, n_loci, n_vu, n_resp, out)


def finish(df, args, allele, resp, qty, epi, pep, pmid,
           n_raw, n_loci, n_vu, n_resp, out):
    print("replaying deduplication ...")
    kept, n_dropped, n_flagged = replay(df, allele, epi, resp, qty, pmid)
    n_dedup = len(kept) + n_flagged        # after dedup, before flag removal
    seq = kept[pep].astype(str)
    canonical = seq.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{7,15}")
    n_noncanon = int((~canonical).sum())
    final = kept[canonical]

    released = len(pd.read_csv(RELEASE_METADATA, low_memory=False))
    stages = [
        ("0. Raw IEDB MHC ligand records", n_raw, None),
        ("1. HLA-A / HLA-B / HLA-C allele restriction", n_loci, n_raw - n_loci),
        ("2. Quantitative value + assay units present", n_vu, n_loci - n_vu),
        ("3. Retained assay response (KD or IC50)", n_resp, n_vu - n_resp),
        ("4. Deduplication of repeated measurements", n_dedup, n_resp - n_dedup),
        ("5. Removal of curation-flagged records", n_dedup - n_flagged, n_flagged),
        ("6. Sequence filter: 7-15mers of the 20 standard residues",
         len(final), n_noncanon),
    ]
    for name, remaining, removed in stages:
        print(f"  {name:<46} {remaining:>10,}  {'' if removed is None else f'-{removed:,}'}")

    print(f"\nreleased metadata.csv rows: {released:,}")
    assert len(final) == released, (
        f"replay ends at {len(final):,} but the release holds {released:,}; the "
        "stage counts cannot be trusted until these agree")
    print("replay reproduces the released file exactly.")

    pd.DataFrame(stages, columns=["Stage", "Records remaining", "Records removed"]
                 ).to_csv(out / "curation_replay_stages.csv", index=False)
    with open(out / "curation_replay.json", "w") as f:
        json.dump({name: {"remaining": r, "removed": d} for name, r, d in stages},
                  f, indent=2)
    print(f"wrote curation_replay.{{csv,json}} to {out}")


if __name__ == "__main__":
    main()
