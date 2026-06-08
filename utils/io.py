"""IO utilities: reading score.sc files, metadata.csv, locating decoys."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, List

import pandas as pd


def read_score_file(score_path: Path) -> pd.DataFrame:
    """Read a Rosetta score.sc file into a DataFrame.

    Score files have the format:
        SEQUENCE: <peptide>
        SCORE: total_score reweighted_sc fa_atr ... description
        SCORE: <numeric values...>          ... <decoy_tag>
        ...

    The first SCORE line is the header; subsequent SCORE lines are data.
    Lines without a SCORE prefix are ignored.
    """
    score_path = Path(score_path)

    header: Optional[List[str]] = None
    rows: List[List[str]] = []

    with open(score_path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("SCORE:"):
                continue
            tokens = line.split()[1:]  # drop the SCORE: prefix
            if header is None:
                header = tokens
            else:
                # If a row has unexpected length (rare, but possible if a
                # value field contained whitespace), pad/truncate to header.
                if len(tokens) != len(header):
                    if len(tokens) < len(header):
                        tokens = tokens + [""] * (len(header) - len(tokens))
                    else:
                        tokens = tokens[: len(header)]
                rows.append(tokens)

    if header is None:
        raise ValueError(f"No SCORE: header found in {score_path}")

    df = pd.DataFrame(rows, columns=header)

    # Coerce numeric columns; leave 'description' as string.
    for col in df.columns:
        if col == "description":
            continue
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def get_decoy_paths(
    pdb_root: Path,
    allele_dir: str,
    peptide: str,
) -> Tuple[Path, List[Path], pd.DataFrame]:
    """Return (peptide_dir, decoy_paths, score_df) for a peptide–allele pair.

    Decoy paths are returned in the order they appear in score.sc, i.e.
    parallel to score_df rows. The 'description' column of score_df gives
    the PDB filename stem.

    Parameters
    ----------
    pdb_root : Path
        Root of the local PDB store, e.g.
        /home/huntek1/main_project/data/IEDB_data_clean/pdb/
    allele_dir : str
        Filesystem-safe allele directory name, e.g. 'A0101'.
    peptide : str
        Peptide sequence, used as the subdirectory name.
    """
    peptide_dir = Path(pdb_root) / allele_dir / peptide
    if not peptide_dir.is_dir():
        raise FileNotFoundError(f"Peptide directory not found: {peptide_dir}")

    score_path = peptide_dir / "score.sc"
    if not score_path.is_file():
        raise FileNotFoundError(f"score.sc not found in {peptide_dir}")

    score_df = read_score_file(score_path)

    decoy_paths: List[Path] = []
    for desc in score_df["description"]:
        pdb_path = peptide_dir / f"{desc}.pdb"
        if not pdb_path.is_file():
            raise FileNotFoundError(f"Decoy PDB not found: {pdb_path}")
        decoy_paths.append(pdb_path)

    return peptide_dir, decoy_paths, score_df


def allele_to_dir_name(allele: str) -> str:
    """Convert allele strings to filesystem-safe form.

    Examples
    --------
    >>> allele_to_dir_name('HLA-A*02:01')
    'A0201'
    >>> allele_to_dir_name('A*02:01')
    'A0201'
    >>> allele_to_dir_name('A0201')
    'A0201'
    """
    s = allele
    if s.startswith("HLA-"):
        s = s[4:]
    return s.replace("*", "").replace(":", "")


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """Load the master metadata.csv."""
    return pd.read_csv(metadata_path)
