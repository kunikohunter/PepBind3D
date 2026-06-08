"""Generate per-allele FASTA files for the 52 validation pairs.

Reads the validation pair list from the structural validation output
(`rmsd_per_pair.csv`), groups by allele, and writes one FASTA per allele
into the regeneration output directory.

Usage:
    python build_validation_fastas.py

Paths are hardcoded at the top; edit them if your layout differs.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

RMSD_CSV    = Path('/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/01_structural/rmsd_per_pair.csv')
OUTPUT_BASE = Path('/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration')
FASTA_DIR   = OUTPUT_BASE / 'fastas'

COL_ALLELE  = 'allele'
COL_PEPTIDE = 'peptide'


def allele_to_dir(allele: str) -> str:
    """Filesystem-safe allele identifier. 'HLA-A*02:01' -> 'A0201'."""
    s = allele
    if s.startswith('HLA-'):
        s = s[4:]
    return s.replace('*', '').replace(':', '')


def main():
    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(RMSD_CSV)
    print(f'Loaded {len(df)} validation pairs from {RMSD_CSV.name}')

    n_alleles = df[COL_ALLELE].nunique()
    print(f'Spans {n_alleles} unique allele(s)')

    summary = []
    for allele, group in df.groupby(COL_ALLELE):
        peptides = sorted(set(group[COL_PEPTIDE]))
        out_path = FASTA_DIR / f'{allele_to_dir(allele)}_validation.fasta'
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
    summary_path = OUTPUT_BASE / 'allele_manifest.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f'\nManifest: {summary_path}')
    print(f'Total peptides: {summary_df["n_peptides"].sum()}')


if __name__ == '__main__':
    main()
