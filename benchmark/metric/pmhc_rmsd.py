"""RMSD scoring for peptide-MHC structure predictions.

Ports the RMSD metric used for the published PepBind3D validation (Figure 2),
whose reference implementation lives on Tungsten at:

    /home/huntek1/main_project/scripts/IEDB_validation/utils/structure.py
    /home/huntek1/main_project/scripts/IEDB_validation/01_structural_validation.ipynb

Environment
-----------
Use the shared "ensemble" env, which has numpy/scipy/pandas and biopython
(confirmed: numpy 1.26.4, biopython 1.87). `bin/activate` does not reliably
work from a non-interactive shell on this cluster, so prefer calling the
interpreter directly:

    /data/p_csb_meiler/huntek1/envs/ensemble/bin/python pmhc_rmsd.py

or, if you do want an activated shell:

    source /data/p_csb_meiler/huntek1/envs/ensemble/bin/activate
    python pmhc_rmsd.py

Provenance
----------
`superpose_on_mhc` and `peptide_backbone_atoms` below are ported verbatim
(algorithm-for-algorithm) from `utils/structure.py::superpose_on_mhc` and
`::peptide_backbone_atoms`. `find_peptide_chain_in_modeled` and
`find_mhc_chain_in_modeled` are ported verbatim from the same-named cells in
`01_structural_validation.ipynb`. See the docstring of each function for the
exact correspondence and any deviation.

What was added on top (fixes, not metric changes; see task write-up):
  1. Explicit (residue_index, atom_name) atom pairing for the peptide
     backbone RMSD, computed only over atoms present in both structures,
     with counts of atoms used/dropped reported. The original silently
     zipped two lists that could be different lengths.
  2. Content-based chain identification is used everywhere (never chain-ID
     based) to survive AlphaFold2's inconsistent A/B vs B/C chain labeling
     across seeds.
  3. Multi-copy support for reference structures with >1 complex in the
     asymmetric unit, with a configurable `copy_policy` ('first' | 'best' |
     'mean'), default 'best'.
  4. Both mmCIF and PDB input are accepted for both prediction and
     reference files.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure

# Standard amino acid three-letter to one-letter mapping.
# Verbatim from utils/structure.py.
_THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

_CIF_PARSER = MMCIFParser(QUIET=True)
_PDB_PARSER = PDBParser(QUIET=True)


# ---------------------------------------------------------------------------
# Loading (DEVIATION: original only read PDB via PDBParser; here we dispatch
# on extension so mmCIF (Boltz/Chai-1/Protenix) and PDB (AF2/Rosetta) both
# work. Purely an I/O extension -- no effect on geometry or RMSD values.)
# ---------------------------------------------------------------------------


def load_structure(path: Path) -> Structure:
    """Parse a PDB or mmCIF file, dispatching on file extension."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        return _CIF_PARSER.get_structure(path.stem, str(path))
    if suffix in (".pdb", ".ent"):
        return _PDB_PARSER.get_structure(path.stem, str(path))
    raise ValueError(f"Unrecognized structure file extension: {path}")


def chain_sequence(chain: Chain) -> str:
    """One-letter sequence of standard residues in a chain. Verbatim port."""
    seq: List[str] = []
    for residue in chain:
        if residue.id[0] != " ":  # skip HETATMs and modified residues
            continue
        resname = residue.get_resname()
        if resname in _THREE_TO_ONE:
            seq.append(_THREE_TO_ONE[resname])
    return "".join(seq)


# ---------------------------------------------------------------------------
# Chain identification in MODELED (prediction) structures.
# Ported verbatim from 01_structural_validation.ipynb.
# ---------------------------------------------------------------------------


def find_peptide_chain_in_modeled(struct: Structure, expected_peptide: str) -> Chain:
    """Find the modeled chain matching the expected peptide sequence.

    Verbatim port of the notebook cell of the same name. This is exactly
    the mechanism needed for fix #2 (chain-label normalisation): it never
    assumes the peptide is chain "B", it matches by sequence content. The
    "fallback to chain B" branch is retained only for parity with the
    original; score_prediction() emits a warning if that fallback fires.
    """
    for chain in struct[0]:
        if chain_sequence(chain) == expected_peptide:
            return chain
    # Fallback to chain B for backward compatibility
    return struct[0]["B"]


def find_mhc_chain_in_modeled(struct: Structure, peptide_chain_id: str) -> Chain:
    """Find the MHC chain (longest non-peptide protein chain).

    Verbatim port of the notebook cell of the same name.
    """
    candidates = [
        c for c in struct[0]
        if c.id != peptide_chain_id and len(chain_sequence(c)) >= 100
    ]
    if not candidates:
        raise ValueError(f"No MHC chain found in {struct.id}")
    return max(candidates, key=lambda c: len(chain_sequence(c)))


# ---------------------------------------------------------------------------
# Chain identification (+ multi-copy grouping) in REFERENCE (experimental)
# structures. NEW: needed for fix #3 (multi-copy support). The original
# utils/structure.py::identify_peptide_chain / identify_mhc_chain handled at
# most one copy (picking "the first by chain ID" when several peptide-length
# chains existed). Here we instead enumerate every copy and pair each
# peptide chain with its spatially nearest MHC-length chain, so a 4-copy
# asymmetric unit (e.g. PDB 8ENH) yields 4 independent (peptide, MHC) pairs
# instead of one arbitrarily chosen pair.
# ---------------------------------------------------------------------------


def _chain_centroid(chain: Chain) -> np.ndarray:
    coords = np.array([atom.coord for residue in chain for atom in residue])
    return coords.mean(axis=0)


def _seq_identity_fraction(candidate_seq: str, reference_seq: str) -> float:
    """Fraction of `reference_seq` matched by the best global alignment.

    Used to score reference MHC-chain candidates against the modeled MHC
    chain's sequence (fix for the b2m/TCR chain-selection bug -- see
    find_copies_in_reference). A cheap, symmetric-enough proxy for identity:
    the number of matched positions in the best ungapped-scoring global
    alignment (pairwise2.align.globalxx), normalised by the reference
    length. Two very differently-sized chains (e.g. b2m vs. a heavy chain)
    score far lower than same-family chains even under a light-alignment
    (mismatch/gap-free) scoring scheme.
    """
    from Bio import pairwise2

    if not candidate_seq or not reference_seq:
        return 0.0
    score = pairwise2.align.globalxx(candidate_seq, reference_seq, score_only=True)
    return float(score) / len(reference_seq)


def find_copies_in_reference(
    structure: Structure,
    expected_peptide: str,
    peptide_length_range: Tuple[int, int] = (8, 15),
    mhc_min_len: int = 100,
    mhc_reference_seq: Optional[str] = None,
) -> List[Tuple[Chain, Chain]]:
    """Enumerate (peptide_chain, mhc_chain) pairs, one per complex copy.

    Peptide chains are identified by exact sequence match to
    `expected_peptide` (content-based, never by chain ID -- see fix #2);
    if none match exactly, falls back to any chain in
    `peptide_length_range`.

    Each peptide chain is paired with an MHC-length (>= mhc_min_len aa)
    chain, without reusing an MHC chain across copies. Selection is a
    two-stage rule, applied in this order:

      1. FILTER: when `mhc_reference_seq` is given (the modeled MHC chain's
         sequence), every remaining candidate is scored by sequence-identity
         to that reference (see `_seq_identity_fraction`) and candidates
         scoring well below the best-scoring one are dropped. This is what
         excludes beta-2-microglobulin and TCR chains -- they are a
         different protein family from the class I heavy chain and score
         far lower -- while every genuine heavy-chain copy in a multi-copy
         asymmetric unit (near-identical sequence to the reference) survives
         the filter.
      2. TIEBREAK: among the surviving (filtered) candidates, the one with
         the nearest centroid to the peptide chain is chosen. This is what
         keeps peptide/heavy-chain pairing correct *within the right copy*
         of a multi-copy ASU -- proximity alone (the original rule) could
         instead lose to a nearer b2m/TCR chain, and identity alone (an
         earlier version of this fix) could pick a heavy chain from the
         *wrong* copy since every copy's heavy chain scores equally well on
         identity.

    When no `mhc_reference_seq` is supplied, falls back to the original
    proximity-only rule (no filter stage).
    """
    chains = list(structure[0])

    pep_candidates = [c for c in chains if chain_sequence(c) == expected_peptide]
    if not pep_candidates:
        pep_candidates = [
            c for c in chains
            if peptide_length_range[0] <= len(chain_sequence(c)) <= peptide_length_range[1]
        ]
    if not pep_candidates:
        raise ValueError(
            f"No peptide chain found in {structure.id} (expected sequence: {expected_peptide})"
        )

    mhc_candidates = [c for c in chains if len(chain_sequence(c)) >= mhc_min_len]
    if not mhc_candidates:
        raise ValueError(f"No MHC-length (>= {mhc_min_len} aa) chain found in {structure.id}")

    # Stage 1 (filter): restrict to chains that are plausibly the same
    # protein as the modeled MHC chain, by sequence identity. A generous
    # relative margin (candidates within 20 percentage points of the best
    # identity score) keeps every true heavy-chain copy in a multi-copy ASU
    # (which should all score close to the top, near-identical sequence)
    # while still rejecting b2m/TCR chains, which score far lower.
    if mhc_reference_seq:
        identity = {c.id: _seq_identity_fraction(chain_sequence(c), mhc_reference_seq)
                    for c in mhc_candidates}
        best_identity = max(identity.values())
        filtered = [c for c in mhc_candidates if identity[c.id] >= best_identity - 0.2]
        if filtered:
            mhc_candidates = filtered

    mhc_available = list(mhc_candidates)
    pairs: List[Tuple[Chain, Chain]] = []
    for pep in pep_candidates:
        # Stage 2 (tiebreak): nearest centroid to the peptide, among the
        # (already content-filtered) surviving candidates.
        pep_centroid = _chain_centroid(pep)
        mhc_available.sort(key=lambda m: np.linalg.norm(pep_centroid - _chain_centroid(m)))
        best_mhc = mhc_available.pop(0)
        pairs.append((pep, best_mhc))

    return pairs


# ---------------------------------------------------------------------------
# Superposition. Verbatim algorithmic port of
# utils/structure.py::superpose_on_mhc.
# ---------------------------------------------------------------------------


def superpose_on_mhc(
    modeled_mhc: Chain,
    experimental_mhc: Chain,
    n_residues: int = 180,  # kept for backward compatibility; ignored
) -> Tuple[Superimposer, float, int]:
    """Superpose modeled MHC onto experimental MHC using sequence-based pairing.

    Verbatim port of utils/structure.py::superpose_on_mhc: global pairwise
    sequence alignment (Bio.pairwise2.align.globalxs, gap-open -2,
    gap-extend -0.5) between the modeled and experimental MHC sequences;
    Ca atoms are paired only where both sides are aligned (no gap) AND the
    residues are identical; requires >= 50 paired Ca atoms; superimposes
    with experimental FIXED and modeled MOVING
    (`Superimposer.set_atoms(exp_atoms, mod_atoms)`).

    DEVIATION (bookkeeping only, not an algorithm change): also returns the
    number of paired Ca atoms as a 3rd tuple element, for reporting
    (n_paired_ca). The original returned a 2-tuple.
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
            f"Too few paired Ca atoms after sequence alignment: {len(mod_atoms)}"
        )

    sup = Superimposer()
    sup.set_atoms(exp_atoms, mod_atoms)  # (fixed, moving)
    return sup, float(sup.rms), len(mod_atoms)


def peptide_backbone_atoms(chain: Chain) -> List:
    """Backbone atoms (N, Ca, C, O) for a peptide chain, in residue order.

    Verbatim port of utils/structure.py::peptide_backbone_atoms.
    """
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


# ---------------------------------------------------------------------------
# NEW: fix #1, explicit atom pairing guard.
#
# The original built two backbone-atom lists independently (one per chain)
# and zipped them positionally:
#     n = min(len(mod_atoms), len(exp_atoms))
#     diffs = [mod_atoms[i].coord - exp_atoms[i].coord for i in range(n)]
# If the modeled and experimental peptide are missing *different* atoms
# (e.g. exp is missing residue 3's O, modeled is missing residue 1's N),
# positional zip silently pairs the wrong atoms together. We instead key
# every backbone atom by (residue_index_in_chain, atom_name) and intersect.
# ---------------------------------------------------------------------------


def _peptide_atom_map(chain: Chain) -> Dict[Tuple[int, str], object]:
    """Map (residue_index, atom_name) -> Atom for standard residues, chain order."""
    atom_map: Dict[Tuple[int, str], object] = {}
    idx = 0
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if residue.get_resname() not in _THREE_TO_ONE:
            continue
        for atom_name in ("N", "CA", "C", "O"):
            if atom_name in residue:
                atom_map[(idx, atom_name)] = residue[atom_name]
        idx += 1
    return atom_map


# ---------------------------------------------------------------------------
# Top-level scoring API.
# ---------------------------------------------------------------------------


def score_prediction(
    pred_path: Path,
    ref_path: Path,
    expected_peptide: str,
    copy_policy: str = "best",
) -> Dict:
    """Score a peptide-MHC prediction against a reference structure.

    Parameters
    ----------
    pred_path : predicted structure file (mmCIF or PDB).
    ref_path : experimental/reference structure file (mmCIF or PDB).
    expected_peptide : expected one-letter peptide sequence, used for
        content-based chain identification in both structures.
    copy_policy : how to combine results across multiple complex copies in
        the reference asymmetric unit. One of:
          'first' -- use only the first (peptide, MHC) copy found.
          'best'  -- use the copy with the lowest peptide_backbone_rmsd (default).
          'mean'  -- average the numeric results across all copies.

    Returns
    -------
    dict with (at least): peptide_backbone_rmsd, n_atom_pairs,
    n_atoms_dropped, n_paired_ca, superposition_rmsd, copy_used,
    n_copies_available, all_copy_rmsds, warnings.
    """
    if copy_policy not in ("first", "best", "mean"):
        raise ValueError(f"Unknown copy_policy: {copy_policy!r}")

    warnings: List[str] = []

    mod_struct = load_structure(pred_path)
    ref_struct = load_structure(ref_path)

    mod_pep = find_peptide_chain_in_modeled(mod_struct, expected_peptide)
    if chain_sequence(mod_pep) != expected_peptide:
        warnings.append(
            f"Modeled peptide chain (id={mod_pep.id!r}) sequence does not exactly "
            f"match expected_peptide; used chain-B fallback."
        )
    mod_mhc = find_mhc_chain_in_modeled(mod_struct, mod_pep.id)
    mod_atom_map = _peptide_atom_map(mod_pep)

    copies = find_copies_in_reference(
        ref_struct, expected_peptide, mhc_reference_seq=chain_sequence(mod_mhc)
    )
    n_copies_available = len(copies)
    if copy_policy == "first":
        copies = copies[:1]

    per_copy: List[Dict] = []
    for exp_pep, exp_mhc in copies:
        try:
            sup, mhc_rmsd, n_paired_ca = superpose_on_mhc(mod_mhc, exp_mhc)
        except ValueError as e:
            warnings.append(f"Copy (pep={exp_pep.id}, mhc={exp_mhc.id}): {e}")
            continue

        rot, tran = sup.rotran
        transformed = {k: (np.dot(a.coord, rot) + tran) for k, a in mod_atom_map.items()}
        exp_atom_map = _peptide_atom_map(exp_pep)

        common_keys = sorted(set(transformed) & set(exp_atom_map))
        n_dropped = len(set(transformed) ^ set(exp_atom_map))
        if not common_keys:
            warnings.append(
                f"Copy (pep={exp_pep.id}, mhc={exp_mhc.id}): no overlapping "
                f"peptide backbone atoms between prediction and reference"
            )
            continue

        diffs = np.array([transformed[k] - exp_atom_map[k].coord for k in common_keys])
        rmsd = float(np.sqrt(np.mean(np.sum(diffs ** 2, axis=1))))

        per_copy.append({
            "peptide_backbone_rmsd": rmsd,
            "n_atom_pairs": len(common_keys),
            "n_atoms_dropped": n_dropped,
            "n_paired_ca": n_paired_ca,
            "superposition_rmsd": mhc_rmsd,
            "copy_used": f"pep={exp_pep.id},mhc={exp_mhc.id}",
        })

    if not per_copy:
        raise ValueError(
            f"No usable complex copy found for {pred_path} vs {ref_path} "
            f"(peptide={expected_peptide}); warnings: {warnings}"
        )

    all_copy_rmsds = [r["peptide_backbone_rmsd"] for r in per_copy]

    if copy_policy in ("first", "best"):
        chosen = min(per_copy, key=lambda r: r["peptide_backbone_rmsd"]) if copy_policy == "best" else per_copy[0]
        result = dict(chosen)
    else:  # mean
        keys = ["peptide_backbone_rmsd", "n_atom_pairs", "n_atoms_dropped",
                "n_paired_ca", "superposition_rmsd"]
        result = {k: float(np.mean([r[k] for r in per_copy])) for k in keys}
        result["copy_used"] = "mean of: " + "; ".join(r["copy_used"] for r in per_copy)

    result["n_copies_available"] = n_copies_available
    result["all_copy_rmsds"] = all_copy_rmsds
    result["warnings"] = warnings
    return result


# ---------------------------------------------------------------------------
# Self-test.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import shutil
    import tempfile

    from Bio.PDB import PDBIO

    REFS_DIR = Path("/data/p_csb_meiler/huntek1/benchmark/refs/cif")
    SINGLE_COPY_CIF = REFS_DIR / "22EE.cif"       # 1 copy: A=heavy(274) B=b2m(100) C=pep(9)
    SINGLE_COPY_PEPTIDE = "SLWQLLQAA"
    MULTI_COPY_CIF = REFS_DIR / "28IL.cif"         # 2 copies, TCR-pMHC, pep = NYNYLYRLF
    MULTI_COPY_PEPTIDE = "NYNYLYRLF"

    tmpdir = Path(tempfile.mkdtemp(prefix="pmhc_rmsd_selftest_"))
    io = PDBIO()

    def write_pdb(structure, path):
        io.set_structure(structure)
        io.save(str(path))

    def fmt(d):
        keys = ["peptide_backbone_rmsd", "n_atom_pairs", "n_atoms_dropped",
                 "n_paired_ca", "superposition_rmsd", "copy_used",
                 "n_copies_available", "all_copy_rmsds", "warnings"]
        return "\n".join(f"    {k}: {d.get(k)}" for k in keys)

    print("=" * 70)
    print("TEST (a): score a structure against itself -> RMSD ~ 0.000")
    print("=" * 70)
    s_self = load_structure(SINGLE_COPY_CIF)
    self_pdb = tmpdir / "22EE_self.pdb"
    write_pdb(s_self, self_pdb)
    res_a = score_prediction(self_pdb, SINGLE_COPY_CIF, SINGLE_COPY_PEPTIDE, copy_policy="best")
    print(fmt(res_a))
    assert res_a["peptide_backbone_rmsd"] < 1e-3, "TEST (a) FAILED"
    assert res_a["superposition_rmsd"] < 1e-3, "TEST (a) FAILED (superposition)"
    print("PASSED\n")

    print("=" * 70)
    print("TEST (b): known rigid-body transform (translate 5A, rotate 30deg)")
    print("           -> peptide RMSD after superposition still ~ 0.000")
    print("=" * 70)
    s_tf = load_structure(SINGLE_COPY_CIF)
    theta = np.radians(30.0)
    R = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta),  np.cos(theta), 0.0],
        [0.0,            0.0,           1.0],
    ], dtype=float)
    t = np.array([5.0, 0.0, 0.0], dtype=float)
    for atom in s_tf.get_atoms():
        atom.coord = np.dot(atom.coord, R.T) + t
    tf_pdb = tmpdir / "22EE_rigid_transformed.pdb"
    write_pdb(s_tf, tf_pdb)
    res_b = score_prediction(tf_pdb, SINGLE_COPY_CIF, SINGLE_COPY_PEPTIDE, copy_policy="best")
    print(fmt(res_b))
    assert res_b["peptide_backbone_rmsd"] < 1e-2, "TEST (b) FAILED"
    print("PASSED\n")

    print("=" * 70)
    print("TEST (c): perturb peptide coords by known amount -> RMSD matches")
    print("           analytical expectation")
    print("=" * 70)
    s_pert = load_structure(SINGLE_COPY_CIF)
    shift = np.array([1.0, 2.0, 2.0], dtype=float)  # magnitude = 3.0
    expected_shift_magnitude = float(np.linalg.norm(shift))
    for chain in s_pert[0]:
        if chain_sequence(chain) == SINGLE_COPY_PEPTIDE:
            for residue in chain:
                for atom in residue:
                    atom.coord = atom.coord + shift
    pert_pdb = tmpdir / "22EE_peptide_perturbed.pdb"
    write_pdb(s_pert, pert_pdb)
    res_c = score_prediction(pert_pdb, SINGLE_COPY_CIF, SINGLE_COPY_PEPTIDE, copy_policy="best")
    print(fmt(res_c))
    print(f"    analytically expected RMSD: {expected_shift_magnitude:.6f}")
    assert abs(res_c["peptide_backbone_rmsd"] - expected_shift_magnitude) < 1e-3, "TEST (c) FAILED"
    assert res_c["superposition_rmsd"] < 1e-3, "TEST (c) FAILED (MHC should be unperturbed)"
    print("PASSED\n")

    print("=" * 70)
    print("TEST (d, bonus): multi-copy reference (28IL, 2 copies) + copy_policy")
    print("=" * 70)
    s_multi = load_structure(MULTI_COPY_CIF)
    multi_pdb = tmpdir / "28IL_self.pdb"
    write_pdb(s_multi, multi_pdb)
    for policy in ("first", "best", "mean"):
        res_d = score_prediction(multi_pdb, MULTI_COPY_CIF, MULTI_COPY_PEPTIDE, copy_policy=policy)
        print(f"  copy_policy={policy}:")
        print(fmt(res_d))
    assert res_d["n_copies_available"] == 2, "TEST (d) FAILED: expected 2 copies"
    print("PASSED\n")

    print("=" * 70)
    print("TEST (e, bonus): mmCIF vs PDB input tolerance")
    print("=" * 70)
    res_e_cif = score_prediction(SINGLE_COPY_CIF, SINGLE_COPY_CIF, SINGLE_COPY_PEPTIDE)
    res_e_pdb = score_prediction(self_pdb, SINGLE_COPY_CIF, SINGLE_COPY_PEPTIDE)
    print("  cif->cif:", res_e_cif["peptide_backbone_rmsd"])
    print("  pdb->cif:", res_e_pdb["peptide_backbone_rmsd"])
    print("PASSED\n")

    shutil.rmtree(tmpdir, ignore_errors=True)
    print("ALL SELF-TESTS PASSED")
