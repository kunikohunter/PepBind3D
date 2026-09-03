#!/usr/bin/env python3
"""
Staging job: run after the authoritative threading chain (13758514 -> ... -> 13758524)
completes. Collects every <peptide>_input.pdb produced under
/data/p_csb_meiler/huntek1/hlac_generator/output/IEDB_data/<allele>/output/...,
verifies each stamps REMARK 220 VERSION 2024.09+release.06b3cf8, and stages it at
/data/p_csb_meiler/huntek1/hlac_production/pdb/<COMPACT>/<peptide>/<peptide>_input.pdb

Hard-fails (nonzero exit) listing every peptide that is missing, version-mismatched,
or otherwise bad, so a downstream --dependency=afterok chain halts rather than docking
against an incomplete or wrong-build staging tree.
"""
import csv, os, sys, shutil, glob

GEN_ROOT = "/data/p_csb_meiler/huntek1/hlac_generator/output/IEDB_data"
STAGE_ROOT = "/data/p_csb_meiler/huntek1/hlac_production/pdb"
TARGETS_CSV = "/home/huntek1/main_project/scripts/hla_c/hlac_threading_targets.csv"
EXPECTED_VERSION = "2024.09+release.06b3cf8"

def compact(allele):
    return allele.replace("*", "").replace(":", "")

def find_input_pdb(allele, peptide):
    # Scope strictly to this allele's own output tree -- 386 of 1861 peptides in the
    # target list are shared across multiple alleles by sequence, so an unscoped
    # filesystem-wide search could silently stage the wrong allele's structure.
    pattern = os.path.join(GEN_ROOT, allele, "output", "*", peptide, f"{peptide}_input.pdb")
    hits = glob.glob(pattern)
    if not hits:
        return None
    return hits[0]

def check_version(pdb_path):
    try:
        with open(pdb_path, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK 220  VERSION"):
                    return EXPECTED_VERSION in line, line.strip()
    except OSError as e:
        return False, f"<read error: {e}>"
    return False, "<no REMARK 220 VERSION line found>"

def main():
    with open(TARGETS_CSV, newline="") as f:
        content = f.read()
    rows = list(csv.DictReader(content.replace("\r\n", "\n").replace("\r", "\n").splitlines()))

    total = len(rows)
    missing = []
    version_fail = []
    staged = 0

    for row in rows:
        allele = row["allele"].strip()
        peptide = row["peptide"].strip()
        comp = compact(allele)

        src = find_input_pdb(allele, peptide)
        if src is None:
            missing.append((allele, peptide))
            continue

        ok, remark = check_version(src)
        if not ok:
            version_fail.append((allele, peptide, src, remark))
            continue

        dest_dir = os.path.join(STAGE_ROOT, comp, peptide)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, f"{peptide}_input.pdb")
        shutil.copy2(src, dest)
        staged += 1

    print(f"Targets in CSV: {total}")
    print(f"Staged OK: {staged}")
    print(f"Missing (no threaded PDB found): {len(missing)}")
    print(f"Version mismatch / unreadable REMARK 220: {len(version_fail)}")

    if missing:
        print("\n--- MISSING PEPTIDES ---")
        for a, p in missing:
            print(f"MISSING\t{a}\t{p}")

    if version_fail:
        print("\n--- VERSION MISMATCHES ---")
        for a, p, src, remark in version_fail:
            print(f"VERSION_FAIL\t{a}\t{p}\t{src}\t{remark}")

    if missing or version_fail:
        print(f"\nFATAL: {len(missing)} missing + {len(version_fail)} version-mismatched "
              f"of {total} targets. Refusing to declare staging complete.", file=sys.stderr)
        sys.exit(1)

    if staged != total:
        print(f"FATAL: staged count {staged} != target count {total}", file=sys.stderr)
        sys.exit(1)

    if staged != 1861:
        print(f"FATAL: expected 1861 staged peptides, got {staged}", file=sys.stderr)
        sys.exit(1)

    print(f"\n=== OK: all {staged} peptides staged and REMARK 220 VERSION-verified "
          f"({EXPECTED_VERSION}) ===")

if __name__ == "__main__":
    main()
