"""Generate per-allele FASTA files for a set of validation pairs.

Reads a pair list (any CSV with `allele` and `peptide` columns), groups by
allele, and writes one FASTA per allele plus a manifest. The FASTAs feed
`run_regeneration.sh`, which re-docks the pairs with `--ignore_epitope_match`
so the ensembles are free of self-template leakage.

Usage:
    python build_validation_fastas.py                       # original 52 pairs
    python build_validation_fastas.py --pairs-csv <csv> --out-base <dir> \
                                      --suffix validation_v2

Note the default `--pairs-csv` was corrected: it pointed at
`01_structural/rmsd_per_pair.csv`, a directory that no longer exists (the
current output lives in `01_structural_regen/`).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT  # noqa: E402


# -----------------------------------------------------------------------------
# Defaults (overridable on the command line)
# -----------------------------------------------------------------------------

RMSD_CSV    = DATA_ROOT / "IEDB_validation" / "01_structural_regen" / "rmsd_per_pair.csv"
OUTPUT_BASE = DATA_ROOT / "IEDB_validation" / "regeneration"

COL_ALLELE  = 'allele'
COL_PEPTIDE = 'peptide'


def allele_to_dir(allele: str) -> str:
    """Filesystem-safe allele identifier. 'HLA-A*02:01' -> 'A0201'."""
    s = allele
    if s.startswith('HLA-'):
        s = s[4:]
    return s.replace('*', '').replace(':', '')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairs-csv', default=str(RMSD_CSV),
                    help='CSV with `allele` and `peptide` columns')
    ap.add_argument('--out-base', default=str(OUTPUT_BASE))
    ap.add_argument('--suffix', default='validation',
                    help='FASTA name suffix: {allele_dir}_{suffix}.fasta')
    args = ap.parse_args()

    pairs_csv = Path(args.pairs_csv)
    output_base = Path(args.out_base)
    FASTA_DIR = output_base / 'fastas'
    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(pairs_csv)
    missing = {COL_ALLELE, COL_PEPTIDE} - set(df.columns)
    if missing:
        raise SystemExit(f'{pairs_csv} is missing required column(s): {sorted(missing)}')
    print(f'Loaded {len(df)} validation pairs from {pairs_csv.name}')

    n_alleles = df[COL_ALLELE].nunique()
    print(f'Spans {n_alleles} unique allele(s)')

    summary = []
    for allele, group in df.groupby(COL_ALLELE):
        peptides = sorted(set(group[COL_PEPTIDE]))
        out_path = FASTA_DIR / f'{allele_to_dir(allele)}_{args.suffix}.fasta'
        with open(out_path, 'w') as f:
            for pep in peptides:
                f.write(f'>{pep}\n{pep}\n')
        summary.append({
            'allele':        allele,
            'allele_dir':    allele_to_dir(allele),
            'n_peptides':    len(peptides),
            'fasta_path':    str(out_path),
        })
        print(f'  {allele:12s}  {len(peptides):3d} peptides  ->  {out_path.name}')

    summary_df = pd.DataFrame(summary)
    summary_path = output_base / f'allele_manifest_{args.suffix}.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f'\nManifest: {summary_path}')
    print(f'Total peptides: {summary_df["n_peptides"].sum()}')


if __name__ == '__main__':
    main()
