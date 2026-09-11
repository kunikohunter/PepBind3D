#!/data/p_csb_meiler/huntek1/envs/ensemble/bin/python
"""Merge per-target score_target.py outputs into arm-level summary/decoys CSVs.

Usage:
    merge_results.py <arm>

Reads:
    cofold_arm/<arm>/results/tmp/*.summary.csv
    cofold_arm/<arm>/results/tmp/*.decoys.csv

Writes:
    cofold_arm/<arm>/results/<arm>_summary.csv
    cofold_arm/<arm>/results/<arm>_decoys.csv
"""
import csv
import sys
from pathlib import Path

COFOLD_ROOT = Path("/data/p_csb_meiler/huntek1/benchmark/cofold_arm")
ARMS = ("af2_reuse", "boltz2", "boltz1", "chai1", "protenix")


def merge_csvs(paths, out_path):
    header = None
    rows = []
    for p in sorted(paths):
        with open(p, newline="") as f:
            reader = csv.reader(f)
            file_header = next(reader, None)
            if file_header is None:
                continue
            if header is None:
                header = file_header
            elif file_header != header:
                print(f"WARNING: header mismatch in {p}, skipping file", file=sys.stderr)
                continue
            for row in reader:
                rows.append(row)
    if header is None:
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    return len(rows)


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <arm>", file=sys.stderr)
        sys.exit(1)
    arm = sys.argv[1]
    if arm not in ARMS:
        print(f"ERROR: unknown arm {arm!r}, expected one of {ARMS}", file=sys.stderr)
        sys.exit(1)

    tmp_dir = COFOLD_ROOT / arm / "results" / "tmp"
    summary_paths = sorted(tmp_dir.glob("*.summary.csv"))
    decoy_paths = sorted(tmp_dir.glob("*.decoys.csv"))

    n_summary_rows = merge_csvs(summary_paths, COFOLD_ROOT / arm / "results" / f"{arm}_summary.csv")
    n_decoy_rows = merge_csvs(decoy_paths, COFOLD_ROOT / arm / "results" / f"{arm}_decoys.csv")

    print(f"Merged {arm}: {len(summary_paths)} targets, {n_summary_rows} summary rows, "
          f"{n_decoy_rows} decoy rows from {len(decoy_paths)} target files.")


if __name__ == "__main__":
    main()
