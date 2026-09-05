"""Structure handling: load PDBs, identify chains, align, compute RMSD."""

from __future__ import annotations

import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure


# Standard amino acid three-letter to one-letter mapping.
_THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

_PARSER = PDBParser(QUIET=True)


# ---------------------------------------------------------------------------
# Basic loading and chain inspection
# ---------------------------------------------------------------------------


def load_structure(pdb_path: Path) -> Structure:
    """Parse a PDB file."""
    return _PARSER.get_structure(Path(pdb_path).stem, str(pdb_path))


def chain_sequence(chain: Chain) -> str:
    """Return one-letter sequence of standard residues in a chain."""
    seq: List[str] = []
    for residue in chain:
        if residue.id[0] != " ":  # skip HETATMs and modified residues
            continue
        resname = residue.get_resname()
        if resname in _THREE_TO_ONE:
            seq.append(_THREE_TO_ONE[resname])
    return "".join(seq)


def get_protein_chains(structure: Structure, min_len: int = 5) -> Dict[str, Chain]:
    """Return all protein chains in model 0 with at least min_len residues."""
    chains: Dict[str, Chain] = {}
    for chain in structure[0]:
        if len(chain_sequence(chain)) >= min_len:
            chains[chain.id] = chain
    return chains


# ---------------------------------------------------------------------------
# Chain identification (peptide + MHC) in experimental structures
# ---------------------------------------------------------------------------


def identify_peptide_chain(
    structure: Structure,
    expected_sequence: Optional[str] = None,
    length_range: Tuple[int, int] = (7, 15),
) -> Chain:
    """Identify the peptide chain in a pHLA structure.

    First tries to match the expected sequence exactly. Falls back to the
    chain whose length falls inside length_range. Raises ValueError if no
    chain matches.
    """
    chains = get_protein_chains(structure)

    if expected_sequence:
        matches = [
            c for c in chains.values() if chain_sequence(c) == expected_sequence
        ]
        if matches:
            return matches[0]

    candidates = [
        c
        for c in chains.values()
        if length_range[0] <= len(chain_sequence(c)) <= length_range[1]
    ]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        # Multiple peptide-length chains (e.g. two copies in the asymmetric
        # unit). Take the first by chain ID for reproducibility.
        return sorted(candidates, key=lambda c: c.id)[0]
    raise ValueError(
        f"No peptide chain found in {structure.id} "
        f"(expected sequence: {expected_sequence})"
    )


def identify_mhc_chain(
    structure: Structure,
    exclude_chain_id: Optional[str] = None,
    length_range: Tuple[int, int] = (170, 290),
) -> Chain:
    """Identify the MHC heavy chain.

    Length range 170–290 catches both the trimmed binding cleft (~180
    residues, α₁/α₂ only) and the full heavy chain (~270 residues, α₁/α₂/α₃).
    Excludes the peptide chain.
    """
    chains = get_protein_chains(structure)
    candidates = [
        c
        for c in chains.values()
        if c.id != exclude_chain_id
        and length_range[0] <= len(chain_sequence(c)) <= length_range[1]
    ]
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        # Multiple matches (e.g. two HLAs in the asymmetric unit, or the
        # second is β2m mis-sized). Return the longest (more likely full chain).
        return max(candidates, key=lambda c: len(chain_sequence(c)))
    raise ValueError(f"No MHC heavy chain found in {structure.id}")


# ---------------------------------------------------------------------------
# Alignment + RMSD
# ---------------------------------------------------------------------------


def get_ca_atoms(chain: Chain, max_residues: Optional[int] = None) -> List:
    """List of Cα atoms for standard residues in a chain, in residue order."""
    atoms = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if residue.get_resname() not in _THREE_TO_ONE:
            continue
        if "CA" not in residue:
            continue
        atoms.append(residue["CA"])
        if max_residues is not None and len(atoms) >= max_residues:
            break
    return atoms


def superpose_on_mhc(
    modeled_mhc: Chain,
    experimental_mhc: Chain,
    n_residues: int = 180,  # kept for backward compatibility; ignored
) -> Tuple[Superimposer, float]:
    """Superpose modeled MHC onto experimental MHC using sequence-based pairing.

    Pairs Cα atoms by identical residue match in a global pairwise sequence
    alignment, then performs least-squares superposition. Robust to different
    residue numbering, missing residues, and chain-length differences.

    Returns the Superimposer (ready to .apply() to atoms) and the alignment RMSD.
    """
    from Bio import pairwise2

    mod_residues = [r for r in modeled_mhc
                    if r.id[0] == ' ' and r.get_resname() in _THREE_TO_ONE
                    and 'CA' in r]
    exp_residues = [r for r in experimental_mhc
                    if r.id[0] == ' ' and r.get_resname() in _THREE_TO_ONE
                    and 'CA' in r]

    mod_seq = ''.join(_THREE_TO_ONE[r.get_resname()] for r in mod_residues)
    exp_seq = ''.join(_THREE_TO_ONE[r.get_resname()] for r in exp_residues)

    # Global alignment, light gap penalties (we expect high identity)
    alignments = pairwise2.align.globalxs(mod_seq, exp_seq, -2, -0.5)
    if not alignments:
        raise ValueError("Sequence alignment failed")
    aligned_mod, aligned_exp, *_ = alignments[0]

    # Walk the alignment, pairing residues only where both sides are aligned
    # to the same amino acid (i.e. no gap and identical)
    mod_atoms, exp_atoms = [], []
    mi = ei = 0
    for mc, ec in zip(aligned_mod, aligned_exp):
        if mc != '-' and ec != '-':
            if mc == ec:
                mod_atoms.append(mod_residues[mi]['CA'])
                exp_atoms.append(exp_residues[ei]['CA'])
            mi += 1
            ei += 1
        elif mc != '-':
            mi += 1
        elif ec != '-':
            ei += 1

    if len(mod_atoms) < 50:
        raise ValueError(
            f"Too few paired Cα atoms after sequence alignment: {len(mod_atoms)}"
        )

    sup = Superimposer()
    sup.set_atoms(exp_atoms, mod_atoms)  # (fixed, moving)
    return sup, float(sup.rms)


def peptide_backbone_atoms(chain: Chain) -> List:
    """Return backbone atoms (N, Cα, C, O) for a peptide chain in residue order."""
    atoms = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if residue.get_resname() not in _THREE_TO_ONE:
            continue
        for atom_name in ("N", "CA", "C", "O"):
            if atom_name in residue:
                atoms.append(residue[atom_name])
    return atoms


def _align_peptide_residues(mod_chain: Chain, exp_chain: Chain) -> List[Tuple[object, object]]:
    """Pair modeled/experimental peptide residues by global sequence alignment.

    Returns a list of (modeled_residue, experimental_residue) for every
    alignment column where both sides are present (no gap) and identical.

    This replaces pairing peptide residues by ordinal resolved-residue
    position (the previous approach in compute_peptide_rmsd, effectively a
    positional zip of two independently-built atom lists). Ordinal pairing
    is correct only when both chains resolve the same set of residues; if
    the experimental peptide is missing an *internal* residue (e.g. an
    unresolved Trp side chain), every residue after the gap shifts down by
    one position on the experimental side only, so every downstream residue
    is silently compared against the wrong partner while the atom count
    still looks like a clean match. Global sequence alignment (same
    parameters as superpose_on_mhc: Bio.pairwise2.align.globalxs, gap-open
    -2, gap-extend -0.5) places a gap at the unresolved residue instead, so
    every residue downstream of it still pairs against its true structural
    counterpart.
    """
    from Bio import pairwise2

    mod_residues = [r for r in mod_chain if r.id[0] == " " and r.get_resname() in _THREE_TO_ONE]
    exp_residues = [r for r in exp_chain if r.id[0] == " " and r.get_resname() in _THREE_TO_ONE]

    mod_seq = "".join(_THREE_TO_ONE[r.get_resname()] for r in mod_residues)
    exp_seq = "".join(_THREE_TO_ONE[r.get_resname()] for r in exp_residues)

    alignments = pairwise2.align.globalxs(mod_seq, exp_seq, -2, -0.5)
    if not alignments:
        raise ValueError("Peptide sequence alignment failed")
    aligned_mod, aligned_exp, *_ = alignments[0]

    pairs: List[Tuple[object, object]] = []
    mi = ei = 0
    for mc, ec in zip(aligned_mod, aligned_exp):
        if mc != "-" and ec != "-":
            if mc == ec:
                pairs.append((mod_residues[mi], exp_residues[ei]))
            mi += 1
            ei += 1
        elif mc != "-":
            mi += 1
        elif ec != "-":
            ei += 1
    return pairs


def _peptide_backbone_atom_count(chain: Chain) -> int:
    """Total number of standard backbone (N, Cα, C, O) atoms present in a peptide chain."""
    total = 0
    for residue in chain:
        if residue.id[0] != " " or residue.get_resname() not in _THREE_TO_ONE:
            continue
        total += sum(1 for atom_name in ("N", "CA", "C", "O") if atom_name in residue)
    return total


def _peptide_atom_pairs(mod_chain: Chain, exp_chain: Chain) -> Tuple[List[Tuple[object, object]], int, int]:
    """Backbone atom pairs (modeled_atom, experimental_atom) via aligned residues.

    Residue correspondence comes from `_align_peptide_residues` (sequence
    alignment), not ordinal position. Within each aligned residue pair, an
    atom is included only if present on both sides (handles atom-level, as
    opposed to residue-level, disorder). Also returns the total backbone
    atom counts on each side (modeled, experimental), for a coverage check.
    """
    residue_pairs = _align_peptide_residues(mod_chain, exp_chain)
    paired: List[Tuple[object, object]] = []
    for mod_res, exp_res in residue_pairs:
        for atom_name in ("N", "CA", "C", "O"):
            if atom_name in mod_res and atom_name in exp_res:
                paired.append((mod_res[atom_name], exp_res[atom_name]))
    n_mod_atoms = _peptide_backbone_atom_count(mod_chain)
    n_exp_atoms = _peptide_backbone_atom_count(exp_chain)
    return paired, n_mod_atoms, n_exp_atoms


# Minimum fraction of the modeled peptide's full backbone atom complement
# (4 atoms/residue) that must survive alignment and be present on both
# sides before a peptide RMSD is trusted. Analogous in spirit to
# superpose_on_mhc's ">= 50 paired Cα atoms" guard: a peptide RMSD computed
# from only a handful of overlapping atoms (e.g. after a bad alignment, or
# a mostly-unresolved reference peptide) is not a meaningful number and
# must fail loudly rather than being silently reported.
_MIN_PEPTIDE_ATOM_COVERAGE = 0.5


def compute_peptide_rmsd(
    modeled_pdb: Path,
    experimental_pdb: Path,
    expected_peptide_seq: str,
    n_align_residues: int = 180,
) -> Dict:
    """Peptide backbone RMSD between modeled and experimental structures.

    Workflow:
      1. Load both PDBs.
      2. Modeled: chain A = MHC, chain B = peptide (project convention).
      3. Experimental: identify peptide chain by sequence, MHC by length.
      4. Superpose on MHC binding cleft (first n_align_residues Cα atoms).
      5. Apply transform to all modeled atoms.
      6. Compute RMSD on peptide backbone atoms (N, Cα, C, O), pairing
         residues by global sequence alignment (not ordinal position) so
         an unresolved internal experimental residue does not shift the
         pairing of every downstream residue.

    Returns dict: rmsd, mhc_alignment_rmsd, n_atoms, n_atoms_dropped,
    modeled_peptide_chain, experimental_peptide_chain.
    """
    mod_struct = load_structure(modeled_pdb)
    exp_struct = load_structure(experimental_pdb)

    mod_mhc = mod_struct[0]["A"]
    mod_pep = mod_struct[0]["B"]

    exp_pep = identify_peptide_chain(exp_struct, expected_peptide_seq)
    exp_mhc = identify_mhc_chain(exp_struct, exclude_chain_id=exp_pep.id)

    sup, mhc_rmsd = superpose_on_mhc(mod_mhc, exp_mhc, n_align_residues)
    sup.apply(list(mod_struct.get_atoms()))

    atom_pairs, n_mod_atoms, n_exp_atoms = _peptide_atom_pairs(mod_pep, exp_pep)
    n = len(atom_pairs)
    n_dropped = (n_mod_atoms - n) + (n_exp_atoms - n)
    if n == 0:
        raise ValueError("No peptide backbone atoms to compare")

    min_required = _MIN_PEPTIDE_ATOM_COVERAGE * 4 * len(expected_peptide_seq)
    if n < min_required:
        raise ValueError(
            f"Only {n} paired peptide backbone atoms (< {min_required:.1f} "
            f"required, {_MIN_PEPTIDE_ATOM_COVERAGE:.0%} of "
            f"{4 * len(expected_peptide_seq)} full backbone atoms) -- "
            f"peptide correspondence is not trustworthy"
        )

    diffs = np.array(
        [mod_atom.coord - exp_atom.coord for mod_atom, exp_atom in atom_pairs]
    )
    rmsd = float(np.sqrt(np.mean(np.sum(diffs ** 2, axis=1))))

    return {
        "rmsd": rmsd,
        "mhc_alignment_rmsd": mhc_rmsd,
        "n_atoms": n,
        "n_atoms_dropped": n_dropped,
        "modeled_peptide_chain": mod_pep.id,
        "experimental_peptide_chain": exp_pep.id,
    }


# ---------------------------------------------------------------------------
# Experimental PDB retrieval
# ---------------------------------------------------------------------------


def get_experimental_pdb(
    pdb_id: str,
    template_dir: Path,
    cache_dir: Path,
) -> Path:
    """Locate or download an experimental PDB structure.

    Lookup order: local template database, then cache, then RCSB download.
    """
    pdb_id = pdb_id.upper()

    template_path = Path(template_dir) / f"{pdb_id}.pdb"
    if template_path.is_file():
        return template_path

    cache_path = Path(cache_dir) / f"{pdb_id}.pdb"
    if cache_path.is_file():
        return cache_path

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, cache_path)
    return cache_path
