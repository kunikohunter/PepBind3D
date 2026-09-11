"""Extract Rosetta score terms from per-pair FlexPepDock score.sc files.

The released silent files carry only the ref2015 energy terms, but the
local pdb/ tree retains the full FlexPepDock score set written by
FlexPepDocking (-pep_refine), including reweighted_sc and I_sc.

Expected layout:
    {PDB_ROOT}/{allele_dir}/{peptide}/score.sc

Outputs (written to OUT_DIR):
    per_decoy_scores.csv  one row per decoy, all score terms
    score_summary.csv     one row per peptide-allele pair, best + mean

Usage:
    python parse_scorefiles.py <PDB_ROOT> <OUT_DIR> [--metrics a,b,c]

Self-test (verifies the duplicate-tag dedup logic against a known answer):
    python parse_scorefiles.py --self-test
"""
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import pandas as pd

# pep_sc is required: Supplementary Table S6 reports it, and omitting it here
# is what made S6 irreproducible between 2026-07 and 2026-09 -- score_summary.csv
# was regenerated without it and the comparison silently dropped from 16 rows to 12.
DEFAULT_METRICS = ("reweighted_sc", "I_sc", "total_score", "pep_sc")


def parse_scorefile(path: Path) -> tuple[pd.DataFrame, int]:
    """Parse one Rosetta score.sc into a DataFrame.

    Line 1 is 'SEQUENCE:'; the header is the first SCORE: line whose second
    token is non-numeric. Later SCORE: lines are data. Handles repeated
    headers (from -overwrite / appended runs) and short/malformed rows.

    Some peptide directories were docked/scored more than once without the
    scorefile being cleared first, so score.sc has the same decoy tag (e.g.
    "AADFPGIAR_input_0001") appearing multiple times with DIFFERENT score
    values -- confirmed this session on 73 released pairs, e.g.
    B3801/YKEPNSIIL had every one of its 25 tags duplicated (50 SCORE rows,
    25 PDBs on disk), and B0702/FPYEGGKVF had just its last tag duplicated
    (26 rows). Only one PDB per tag exists on disk, and it reflects
    whichever run wrote it last, so the LAST-appearing SCORE row per tag is
    kept (a dict keyed by tag naturally does this: a later assignment
    overwrites an earlier one for the same key) and earlier duplicates are
    discarded rather than both being fed into the best/mean aggregation.

    Returns (dataframe, n_duplicate_tags_dropped).
    """
    header, rows_by_tag = None, {}
    n_duplicates = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith("SCORE:"):
                continue
            fields = line.split()[1:]
            if not fields:
                continue
            if header is None:
                if not _is_number(fields[0]):
                    header = fields
                continue
            if fields[0] == header[0]:      # repeated header block
                continue
            if len(fields) != len(header):  # truncated / interrupted write
                continue
            tag = fields[-1]
            if tag in rows_by_tag:
                n_duplicates += 1
            rows_by_tag[tag] = fields        # later occurrence wins

    if header is None or not rows_by_tag:
        return pd.DataFrame(), 0

    df = pd.DataFrame(rows_by_tag.values(), columns=header)
    for c in df.columns:
        if c != "description":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, n_duplicates


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def self_test() -> None:
    """Construct a synthetic score.sc with a known-correct answer: two tags
    duplicated (one twice, one three times) with deliberately different
    total_score values each time, plus one tag appearing only once. The
    fix must keep exactly the LAST occurrence of each duplicated tag and
    report the right duplicate count."""
    content = (
        "SEQUENCE:\n"
        "SCORE:     total_score  description\n"
        "SCORE:         -100.0  pep_input_0001\n"
        "SCORE:         -200.0  pep_input_0001\n"  # dup 1 of tag 0001 (should be discarded)
        "SCORE:         -999.0  pep_input_0002\n"  # only occurrence of tag 0002
        "SCORE:          -10.0  pep_input_0003\n"
        "SCORE:          -20.0  pep_input_0003\n"  # dup 1 of tag 0003 (should be discarded)
        "SCORE:          -30.0  pep_input_0003\n"  # dup 2 -- this is the one that should survive
    )
    with tempfile.TemporaryDirectory() as tmp:
        sc_path = Path(tmp) / "score.sc"
        sc_path.write_text(content)
        df, n_dup = parse_scorefile(sc_path)

    # tag 0001 appears twice (1 duplicate), tag 0003 appears three times (2 duplicates) = 3 total
    assert n_dup == 3, f"self-test FAILED: expected 3 duplicate rows dropped, got {n_dup}"
    assert len(df) == 3, f"self-test FAILED: expected 3 unique tags in output, got {len(df)}"

    by_tag = df.set_index("description")["total_score"].to_dict()
    assert by_tag["pep_input_0001"] == -200.0, (
        f"self-test FAILED: tag 0001 should keep its LAST value -200.0, got {by_tag['pep_input_0001']}"
    )
    assert by_tag["pep_input_0002"] == -999.0, (
        f"self-test FAILED: non-duplicated tag 0002 should be untouched, got {by_tag['pep_input_0002']}"
    )
    assert by_tag["pep_input_0003"] == -30.0, (
        f"self-test FAILED: tag 0003 should keep its LAST (3rd) value -30.0, got {by_tag['pep_input_0003']}"
    )

    print("Self-test PASSED: duplicate-tag dedup keeps the last occurrence per tag "
          f"(3 unique tags survived, 3 stale duplicate rows correctly dropped).")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb_root", type=Path, nargs="?")
    ap.add_argument("out_dir", type=Path, nargs="?")
    ap.add_argument("--metrics", default=",".join(DEFAULT_METRICS),
                    help="comma-separated score terms to summarize")
    ap.add_argument("--expect-decoys", type=int, default=25)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if args.pdb_root is None or args.out_dir is None:
        raise SystemExit(__doc__)

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    scorefiles = sorted(args.pdb_root.glob("*/*/score.sc"))
    print(f"Found {len(scorefiles):,} score.sc files under {args.pdb_root}")
    if not scorefiles:
        raise SystemExit("No score.sc files matched {allele_dir}/{peptide}/score.sc")

    frames, empty = [], []
    pairs_with_duplicates = []
    total_duplicate_tags = 0
    for i, sf in enumerate(scorefiles, 1):
        df, n_dup = parse_scorefile(sf)
        if df.empty:
            empty.append(str(sf))
            continue
        if n_dup:
            pairs_with_duplicates.append((sf.parent.parent.name, sf.parent.name, n_dup))
            total_duplicate_tags += n_dup
        df["allele_dir"] = sf.parent.parent.name
        df["peptide"] = sf.parent.name
        frames.append(df)
        if i % 2000 == 0:
            print(f"  {i:,}/{len(scorefiles):,}")

    if not frames:
        raise SystemExit("No SCORE: data rows parsed.")

    decoys = pd.concat(frames, ignore_index=True)
    print(f"\nParsed {len(decoys):,} decoys from {len(frames):,} files")
    if empty:
        print(f"WARNING: {len(empty)} score.sc files had no data rows; first few:")
        for e in empty[:5]:
            print("   ", e)
    if pairs_with_duplicates:
        print(f"\n{len(pairs_with_duplicates):,} pairs had duplicate-tag SCORE rows in score.sc "
              f"(score.sc appended by more than one run) -- kept the last occurrence per tag, "
              f"dropped {total_duplicate_tags:,} stale duplicate rows total:")
        for allele_dir, peptide, n_dup in sorted(pairs_with_duplicates, key=lambda x: -x[2])[:10]:
            print(f"    {allele_dir}/{peptide}: {n_dup} duplicate tag(s)")
        if len(pairs_with_duplicates) > 10:
            print(f"    ... and {len(pairs_with_duplicates) - 10} more")

    missing = [m for m in metrics if m not in decoys.columns]
    if missing:
        print(f"\nWARNING: requested metrics absent: {missing}")
        metrics = [m for m in metrics if m in decoys.columns]
    if not metrics:
        raise SystemExit("None of the requested metrics are present.")
    print(f"Summarizing: {metrics}")

    decoys.to_csv(args.out_dir / "per_decoy_scores.csv", index=False)

    g = decoys.groupby(["allele_dir", "peptide"], sort=True)
    summary = pd.DataFrame({"n_decoys": g.size()})
    for m in metrics:
        summary[f"{m}_best"] = g[m].min()    # Rosetta: lower is better
        summary[f"{m}_mean"] = g[m].mean()
    summary = summary.reset_index()

    out = args.out_dir / "score_summary.csv"
    summary.to_csv(out, index=False)
    print(f"\nWrote {len(summary):,} pair summaries -> {out}")
    print(summary.head().to_string(index=False))

    bad = summary[summary["n_decoys"] != args.expect_decoys]
    if len(bad):
        print(f"\nWARNING: {len(bad):,} pairs != {args.expect_decoys} decoys")
        print(bad.head(10).to_string(index=False))
    else:
        print(f"\nAll pairs have exactly {args.expect_decoys} decoys.")

    print("\nMetric correlations across decoys (sanity check):")
    print(decoys[metrics].corr(method="spearman").round(3).to_string())


if __name__ == "__main__":
    main()
