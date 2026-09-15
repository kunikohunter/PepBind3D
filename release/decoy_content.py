"""
Decoy content validation: does a decoy PDB actually hold the receptor and the
peptide it claims?

Shared by release/convert_to_silent.py (as a gate before a silent file is
written) and analysis/screen_decoy_content.py (retrospectively, over what is
already on disk). One implementation deliberately, so the release gate and the
audit cannot drift apart.

Pure stdlib on purpose: the converter runs in a PyRosetta environment and is
fanned out over 112k pairs, so this must not drag in biopython or matplotlib.
That rules out importing `utils`, whose __init__ eagerly imports both.

The defects this catches come from two independent bugs in the generation
pipeline, both diagnosed on the ACCRE side and both confirmed here from the
residue sequences:

  * the receptor trim builds DeleteRegionMover with end=f"{count}{chain}",
    using a residue COUNT as a residue NUMBER. When the MHC chain has numbering
    gaps, count < max residue number, so every residue numbered above the count
    SURVIVES the trim and becomes its own chain. The remnant is therefore
    untrimmed RECEPTOR (e.g. B5801/TRTSPNIPK carried chain B =
    "HVQHEGLPKPLTLRWEP", alpha2/alpha3 junction sequence), not template peptide
    overhang. The remnant's size equals the number of numbering gaps, which
    predicts every observed case.
  * Rosetta drops occupancy=0 atoms by default, so a template with an
    unresolved peptide middle loads short and the query threads onto fewer
    positions, truncating it.

FlexPepDock then refines whichever chain it takes to be the peptide, which for
an extra_chain pair is the remnant -- so the real peptide never moves and the
25 "decoys" share one pose.

Two invariants follow from the trim being applied to the HEAVY CHAIN only
(DeleteRegionMover start=receptor_rescount+1 end=template_rescount on mhcChain),
and both are useful as diagnostics rather than just description:

  * the remnant SIZE equals the heavy chain's numbering-gap count (8ZV9 11,
    3REW 4, 1JF1 2, 1QR1 0 -- and 1QR1 accordingly never produced an
    extra_chain, only truncation);
  * the remnant SEQUENCE is always a C-terminal heavy-chain fragment. A remnant
    that looked peptide-like would therefore indicate a DIFFERENT bug, not this
    one.

The two classes are mechanically disjoint -- gapped heavy-chain numbering causes
extra_chain, occupancy-0 peptide residues cause truncation -- so a pair can be
attributed to one cause or the other after the fact, and no pair needed both
fixes for the same reason.

  extra_chain       3+ chains. If the remnant is too small to touch the
                    receptor every interface term is exactly 0.000, which is
                    score-visible; if it is large enough, I_sc looks perfectly
                    normal but describes the WRONG chain, which no score-based
                    check can see.
  truncated_peptide 2 chains, but chain B holds only a prefix of the curated
                    peptide -- the structure does not contain the sequence it
                    claims.

See the threading step for where this defect originates, and below for why
chain composition
rather than the score terms is the reliable test.
"""
import collections

RECEPTOR_CHAIN = "A"
PEPTIDE_CHAIN = "B"
# the modeled receptor is the alpha1/alpha2 cleft; 181 residues in practice
RECEPTOR_RANGE = (170, 200)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

OK = "ok"


def chain_sequences(path):
    """{chain_id: one-letter sequence} from CA records, in file order.

    Parsed by fixed column positions rather than whitespace splitting: PDB is a
    fixed-column format and adjacent fields can abut in wide files.

    Residues are deduplicated on (chain, resSeq+iCode) so a residue with
    alternate conformations counts once. Rosetta output does not need this --
    measured over 400 released decoys, all 76,122 CA records carry a blank
    altLoc and there are no duplicate keys -- but the same function is useful
    against experimental structures, which do carry altLocs, and there
    double-counting would silently inflate a chain and turn a truncated peptide
    into an apparent match. Matches the dedup in ACCRE's independent screen.
    """
    seqs = collections.OrderedDict()
    seen = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[12:16] != " CA ":
                continue
            key = (line[21], line[22:27])
            if key in seen:
                continue
            seen.add(key)
            seqs.setdefault(line[21], []).append(
                THREE_TO_ONE.get(line[17:20].strip().upper(), "X"))
    return {k: "".join(v) for k, v in seqs.items()}


def classify(seqs, peptide):
    """Return (verdict, detail). verdict == OK means the decoy holds exactly
    the receptor on chain A and the complete curated peptide on chain B."""
    if not seqs:
        return "unreadable", "no CA records"
    chains = list(seqs)
    pep_seq = seqs.get(PEPTIDE_CHAIN, "")
    n_rec = len(seqs.get(RECEPTOR_CHAIN, ""))

    if len(chains) > 2:
        extra = {c: len(seqs[c]) for c in chains
                 if c not in (RECEPTOR_CHAIN, PEPTIDE_CHAIN)}
        holder = [c for c in chains if seqs[c] == peptide]
        return ("extra_chain",
                f"chains={ {c: len(s) for c, s in seqs.items()} }; extra={extra}; "
                f"peptide_on={holder or 'none'}")
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        return "wrong_chain_ids", f"chains={chains}"
    if not (RECEPTOR_RANGE[0] <= n_rec <= RECEPTOR_RANGE[1]):
        return "receptor_size", f"chain A has {n_rec} residues"
    if pep_seq == peptide:
        return OK, ""
    if len(pep_seq) < len(peptide) and peptide.startswith(pep_seq):
        return ("truncated_peptide",
                f"chain B = {pep_seq} ({len(pep_seq)}/{len(peptide)}), "
                f"missing {peptide[len(pep_seq):]}")
    if len(pep_seq) != len(peptide):
        return "peptide_length", f"chain B = {pep_seq} ({len(pep_seq)} vs {len(peptide)})"
    return "peptide_mismatch", f"chain B = {pep_seq}, expected {peptide}"


def check_decoy(path, peptide):
    """Convenience wrapper: (verdict, detail) for one decoy PDB."""
    try:
        seqs = chain_sequences(path)
    except OSError as e:
        return "unreadable", str(e)
    return classify(seqs, peptide)


def self_test():
    """Known-answer test on synthetic PDB text, one case per defect class."""
    import tempfile
    from pathlib import Path

    def pdb(chain_res, altloc_dup=False):
        out, i = [], 1
        for ch, residues in chain_res:
            for r in residues:
                out.append(f"ATOM  {i:5d}  CA  {r} {ch}{i:4d}"
                           f"      0.000   0.000   0.000  1.00  0.00           C")
                if altloc_dup:
                    # same residue number, alternate conformation 'B'
                    out.append(f"ATOM  {i:5d}  CA B{r} {ch}{i:4d}"
                               f"      0.000   0.000   0.000  1.00  0.00           C")
                i += 1
        return "\n".join(out) + "\n"

    rec = ["ALA"] * 181
    cases = {
        "ok":       ([("A", rec), ("B", ["VAL", "VAL", "ALA", "ASN"])], "VVAN", OK),
        "trunc":    ([("A", rec), ("B", ["VAL", "VAL"])], "VVAN", "truncated_peptide"),
        "extra":    ([("A", rec), ("B", ["TRP", "GLU"]),
                      ("C", ["VAL", "VAL", "ALA", "ASN"])], "VVAN", "extra_chain"),
        "mismatch": ([("A", rec), ("B", ["GLY"] * 4)], "VVAN", "peptide_mismatch"),
        "recsize":  ([("A", ["ALA"] * 20), ("B", ["VAL", "VAL", "ALA", "ASN"])],
                     "VVAN", "receptor_size"),
    }
    with tempfile.TemporaryDirectory() as td:
        for name, (chain_res, peptide, expected) in cases.items():
            p = Path(td) / f"{name}.pdb"
            p.write_text(pdb(chain_res))
            got, detail = check_decoy(p, peptide)
            assert got == expected, f"{name}: got {got} ({detail}), want {expected}"
        # the extra-chain case must say where the real peptide actually is
        assert "peptide_on=['C']" in check_decoy(Path(td) / "extra.pdb", "VVAN")[1]
        # an empty file must not pass as ok
        e = Path(td) / "empty.pdb"; e.write_text("")
        assert check_decoy(e, "VVAN")[0] == "unreadable"

        # altLoc: every residue written twice as alternate conformations must
        # still read as its true length. Without dedup the truncated case below
        # would count 4 residues and pass as ok -- a FALSE NEGATIVE, the
        # dangerous direction.
        a = Path(td) / "altloc_ok.pdb"
        a.write_text(pdb([("A", rec), ("B", ["VAL", "VAL", "ALA", "ASN"])],
                         altloc_dup=True))
        assert check_decoy(a, "VVAN") == (OK, ""), check_decoy(a, "VVAN")
        b = Path(td) / "altloc_trunc.pdb"
        b.write_text(pdb([("A", rec), ("B", ["VAL", "VAL"])], altloc_dup=True))
        assert check_decoy(b, "VVAN")[0] == "truncated_peptide", check_decoy(b, "VVAN")
    print("decoy_content self-test PASSED: healthy decoys pass; extra chains, "
          "truncated peptides, wrong residues, bad receptor size and empty files "
          "are each classified correctly, and the real peptide's chain is named.")


if __name__ == "__main__":
    self_test()
