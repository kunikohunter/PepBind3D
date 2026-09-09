#!/usr/bin/env python3
"""Rebuild docking inputs for targets whose peptide was split into two chains.

DEFECT. Threading places the target peptide onto a peptide template. When the
template is LONGER than the target peptide, the surplus template residues were
left in the model and assigned their own chain, giving three chains instead of
two:

    healthy   A: 1-181 receptor   B: 182-190 peptide
    broken    A: 1-181 receptor   B: 182-183 remnant   C: 184-193 peptide

FlexPepDock then refines the chain it takes to be the peptide -- the remnant --
and never moves the real peptide. All 25 decoys therefore share one identical
peptide backbone RMSD (spread 0.0000 A, against a median of 1.014 A across
unaffected targets), so the "ensemble" is the threaded starting pose repeated
25 times. Detection by all-zero interface terms finds only the subset where the
remnant is too small to form an interface (9O5S, 9SKO, 9XMR); 8RJH and 8ZV9
have normal-looking interface terms for the WRONG chain. Chain count is the
reliable test, decoy spread the confirmation.

FIX. Drop the remnant chain, relabel the peptide chain to B, and renumber
residues consecutively from 1. The peptide coordinates are taken unchanged from
decoy 1: since the peptide never moved, every decoy carries the identical
threaded pose, so this reconstructs the intended input rather than inventing
one. Receptor backbone is likewise unmoved under -pep_refine.

Hydrogens are stripped; Rosetta rebuilds them.
"""
import collections
import glob
import os
import sys

BENCH = "/data/p_csb_meiler/huntek1/benchmark"
ARM = f"{BENCH}/rosetta_arm/arm"
OUT = f"{BENCH}/rosetta_arm/refix"


def chain_ca_counts(path):
    counts = collections.OrderedDict()
    for line in open(path):
        if line.startswith("ATOM") and line[12:16] == " CA ":
            counts.setdefault(line[21], 0)
            counts[line[21]] += 1
    return counts


def rebuild(src, peptide_seq, dst):
    """Write receptor (chain A) + peptide (chain B), renumbered from 1."""
    counts = chain_ca_counts(src)
    if len(counts) < 3:
        return None
    # The peptide chain is the one whose residue count matches the target
    # peptide. Never guess by position: in 8RJH/8ZV9 the remnant (11) is LONGER
    # than the peptide (9), so "last chain" and "longest chain" disagree.
    pep_chain = [c for c, n in counts.items() if n == len(peptide_seq) and c != "A"]
    if len(pep_chain) != 1:
        return f"ambiguous peptide chain: {counts}, peptide len {len(peptide_seq)}"
    pep_chain = pep_chain[0]

    out_lines = []
    serial = 0
    for new_chain, keep in (("A", "A"), ("B", pep_chain)):
        resnum = 0
        last_res = None
        for line in open(src):
            if not line.startswith("ATOM") or line[21] != keep:
                continue
            if line[76:78].strip() == "H" or line[12:16].strip().startswith(("H", "1H", "2H", "3H")):
                continue
            res_id = line[22:27]
            if res_id != last_res:
                resnum += 1
                last_res = res_id
            serial += 1
            out_lines.append(
                line[:6] + f"{serial:5d}" + line[11:21] + new_chain
                + f"{resnum:4d}" + line[26:]
            )
        out_lines.append("TER\n")
    out_lines.append("END\n")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        fh.writelines(out_lines)
    return None


def main():
    targets = sys.argv[1:]
    if not targets:
        sys.exit("usage: fix_split_peptide_inputs.py <pdb_id> [...]")
    import csv
    jl = {r["pdb_id"]: r["peptide_seq"] for r in
          csv.DictReader(open(f"{BENCH}/cofold_arm/joblist.tsv"), delimiter="\t")}
    for t in targets:
        pep = jl[t]
        src = sorted(glob.glob(f"{ARM}/{t}/docking/{pep}_input_0001.pdb"))
        if not src:
            print(f"{t}: no decoy 1 found")
            continue
        dst = f"{OUT}/{t}/{pep}_input.pdb"
        err = rebuild(src[0], pep, dst)
        if err:
            print(f"{t}: {err}")
            continue
        counts = chain_ca_counts(dst)
        ok = list(counts.values()) == [181, len(pep)]
        print(f"{t}: wrote {dst}  chains={dict(counts)}  {'OK' if ok else 'CHECK'}")


if __name__ == "__main__":
    main()
