#!/usr/bin/env python3
"""M7: how much does the receptor actually vary, and would varying it help?

Reviewer 2 read the figures correctly: each allele has ONE AlphaFold2 receptor
model, and FlexPepDock -pep_refine samples the peptide against it, so decoys
within a pair differ almost entirely in the peptide. Two measurements:

  PART 1 -- receptor spread WITHIN a docking run.
      CA RMSD of the MHC chain between decoy 1 and decoys 2..25, per target.
      Establishes the magnitude of the limitation rather than asserting it.

  PART 2 -- is there real groove variation to be had?
      For every allele with two or more crystal references in the benchmark,
      the pairwise CA RMSD between those crystals over the groove-forming
      region, plus the deviation of the threaded AF2 receptor from each. If
      crystals of one allele differ by less than the docking error, then
      sampling the receptor cannot help and the rigid-receptor choice costs
      little. If they differ by more, it is a real limitation.

Groove region is the alpha1/alpha2 platform, which is what the threaded receptor
covers; the alpha3 domain and beta-2-microglobulin are excluded because the
peptide never contacts them.

RESIDUE PAIRING. Cross-structure comparisons go through
pmhc_rmsd.superpose_on_mhc, which pairs residues by SEQUENCE ALIGNMENT. Do not
pair the first n CA atoms by index: the threaded receptor is 181 residues while
the crystal heavy chains are 263-276 with missing loops, so index i is not the
same residue in both. A first attempt at this script did exactly that and
returned a median AF2-vs-crystal RMSD of 3.841 A with an IQR of 3.815-3.872 --
a 0.057 A spread over 167 independent comparisons, which is the signature of a
systematic offset rather than a measurement. Same failure class as the ordinal
peptide-atom pairing already fixed in pmhc_rmsd.py.
"""
import csv, glob, os, sys, itertools, collections
import numpy as np
sys.path.insert(0, "/data/p_csb_meiler/huntek1/benchmark/metric")
import pmhc_rmsd as M

B = "/data/p_csb_meiler/huntek1/benchmark"
REFS = f"{B}/refs/cif"
JOBLIST = f"{B}/cofold_arm/joblist.tsv"
RECEPTOR_DIR = "/home/huntek1/Data/MHC_database/default_receptor"
GROOVE_MAX_RES = 180


def ca_coords(chain, limit=GROOVE_MAX_RES):
    out = []
    for i, r in enumerate([x for x in chain if x.id[0] == " "]):
        if i >= limit:
            break
        if "CA" in r:
            out.append(np.array(r["CA"].coord, dtype=float))
    return np.array(out)


def kabsch_rmsd(P, Q):
    n = min(len(P), len(Q))
    P, Q = P[:n], Q[:n]
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, W = np.linalg.svd(Pc.T @ Qc)
    if np.linalg.det(V @ W) < 0:
        V[:, -1] = -V[:, -1]
    U = V @ W
    return float(np.sqrt(np.mean(np.sum((Pc @ U - Qc) ** 2, axis=1)))), n


def part1_within_run():
    print("=== PART 1: receptor spread across the 25 decoys of one docking run ===")
    rows = []
    dirs = sorted(glob.glob(f"{B}/rosetta_arm/arm/*/docking"))
    for d in dirs:
        target = d.split("/")[-2]
        files = sorted(glob.glob(f"{d}/*_input_[0-9][0-9][0-9][0-9].pdb"))
        if len(files) < 5:
            continue
        try:
            ref = M.load_structure(files[0])
            rc = [c for c in ref[0] if len(M.chain_sequence(c)) >= 100]
            if not rc:
                continue
            base = ca_coords(rc[0])
            vals = []
            for f in files[1:]:
                s = M.load_structure(f)
                sc = [c for c in s[0] if len(M.chain_sequence(c)) >= 100]
                if not sc:
                    continue
                r, _ = kabsch_rmsd(ca_coords(sc[0]), base)
                vals.append(r)
            if vals:
                rows.append((target, float(np.median(vals)), float(np.max(vals))))
        except Exception as e:
            print(f"  skip {target}: {type(e).__name__}: {e}")
    if rows:
        med = np.array([r[1] for r in rows])
        mx = np.array([r[2] for r in rows])
        print(f"  targets: {len(rows)}")
        print(f"  receptor CA RMSD vs decoy 1, median across decoys: "
              f"median {np.median(med):.4f} A, 95th pct {np.percentile(med,95):.4f} A")
        print(f"  worst single decoy per target: median {np.median(mx):.4f} A, "
              f"max {mx.max():.4f} A")
    return rows


def part2_between_crystals():
    print("\n=== PART 2: groove variation between crystal structures of one allele ===")
    jl = list(csv.DictReader(open(JOBLIST), delimiter="\t"))
    by_allele = collections.defaultdict(list)
    for r in jl:
        by_allele[r["allele"]].append((r["pdb_id"], r["peptide_seq"]))

    xtal_rows, af2_rows = [], []
    for allele, entries in sorted(by_allele.items()):
        if len(entries) < 2:
            continue
        coords = {}
        chains = {}
        for pdb_id, pep in entries:
            path = f"{REFS}/{pdb_id}.cif"
            if not os.path.exists(path):
                continue
            try:
                exp = M.load_structure(path)
                # pick the reference MHC the same way the metric does
                probe = sorted(glob.glob(
                    f"{B}/cofold_arm/boltz2/out/{pdb_id}/boltz_results_pmhc_input/"
                    f"predictions/pmhc_input/*.cif"))
                if not probe:
                    continue
                ps = M.load_structure(probe[0])
                mp = M.find_peptide_chain_in_modeled(ps, pep)
                mm = M.find_mhc_chain_in_modeled(ps, mp.id)
                cps = M.find_copies_in_reference(
                    exp, pep, mhc_reference_seq=M.chain_sequence(mm))
                if not cps:
                    continue
                coords[pdb_id] = ca_coords(cps[0][1])
                chains[pdb_id] = cps[0][1]
            except Exception:
                continue
        if len(coords) < 2:
            continue
        ids = sorted(coords)
        for a, b in itertools.combinations(ids, 2):
            try:
                _, r, n = M.superpose_on_mhc(chains[a], chains[b])
            except ValueError:
                continue
            xtal_rows.append((allele, a, b, r, n))
        rec = f"{RECEPTOR_DIR}/{allele}.pdb"
        if os.path.exists(rec):
            try:
                rs = M.load_structure(rec)
                rch = [c for c in rs[0] if len(M.chain_sequence(c)) >= 100]
                if rch:
                    for pid in ids:
                        try:
                            _, r, n = M.superpose_on_mhc(rch[0], chains[pid])
                        except ValueError:
                            continue
                        af2_rows.append((allele, pid, r, n))
            except Exception:
                pass

    if xtal_rows:
        v = np.array([r[3] for r in xtal_rows])
        print(f"  crystal-vs-crystal, same allele: {len(xtal_rows)} pairs across "
              f"{len(set(r[0] for r in xtal_rows))} alleles")
        print(f"    median {np.median(v):.3f} A   IQR {np.percentile(v,25):.3f}-"
              f"{np.percentile(v,75):.3f}   max {v.max():.3f}")
    if af2_rows:
        v = np.array([r[2] for r in af2_rows])
        print(f"  threaded AF2 receptor vs crystal: {len(af2_rows)} comparisons")
        print(f"    median {np.median(v):.3f} A   IQR {np.percentile(v,25):.3f}-"
              f"{np.percentile(v,75):.3f}   max {v.max():.3f}")

    with open(f"{B}/receptor/receptor_variation.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "allele", "a", "b", "ca_rmsd", "n_ca"])
        for allele, a, b, r, n in xtal_rows:
            w.writerow(["crystal_vs_crystal", allele, a, b, f"{r:.4f}", n])
        for allele, pid, r, n in af2_rows:
            w.writerow(["af2_receptor_vs_crystal", allele, "af2", pid, f"{r:.4f}", n])
    return xtal_rows, af2_rows


if __name__ == "__main__":
    part1_within_run()
    part2_between_crystals()
