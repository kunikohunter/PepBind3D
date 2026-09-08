"""Per-residue conformational spread across each ensemble.

The accuracy analysis showed error rising from the anchors toward the peptide
centre in every method, and the reference crystals are themselves less ordered
there (mean B-factor 66.8 at anchors, 98.8 centrally). That leaves a question
the error analysis cannot answer: are the 25 samples failing to converge on the
bulge, or are they representing a region that genuinely varies?

For each target this superposes every sample onto the reference MHC, then
measures the spread of each peptide residue's CA position across the ensemble
(RMSF about the ensemble mean). If spread tracks distance from the nearer
terminus, and tracks the reference B-factor, the ensemble is recovering real
disorder rather than merely disagreeing with itself.
"""
import csv, os, re, sys, glob, collections
import numpy as np
sys.path.insert(0, "/data/p_csb_meiler/huntek1/benchmark/metric")
from scipy.stats import spearmanr
import pmhc_rmsd as M

B = "/data/p_csb_meiler/huntek1/benchmark"
REFS = f"{B}/refs/cif"
ARM = sys.argv[1]

SUMMARY = {
    "Rosetta":    f"{B}/rosetta_arm/results/arm_summary.csv",
    "AlphaFold2": f"{B}/cofold_arm/af2_reuse/results/af2_reuse_summary.csv",
    "Boltz-2":    f"{B}/cofold_arm/boltz2/results/boltz2_summary.csv",
    "Boltz-1":    f"{B}/cofold_arm/boltz1/results/boltz1_summary.csv",
    "Chai-1":     f"{B}/cofold_arm/chai1/results/chai1_summary.csv",
    "Protenix":   f"{B}/cofold_arm/protenix/results/protenix_summary.csv",
}[ARM]


def sample_files(target):
    if ARM == "Rosetta":
        return sorted(glob.glob(f"{B}/rosetta_arm/arm/{target}/docking/*_input_[0-9][0-9][0-9][0-9].pdb"))
    if ARM == "Chai-1":
        return sorted(glob.glob(f"{B}/cofold_arm/chai1/out/{target}/structures/**/*.cif", recursive=True))
    if ARM == "Protenix":
        # Only the prediction CIFs. A bare recursive glob would also pick up
        # the MSA directory and the input JSON echo.
        return sorted(glob.glob(
            f"{B}/cofold_arm/protenix/out/{target}/pmhc_input/seed_101/predictions/pmhc_input_sample_*.cif"))
    d = {"AlphaFold2": "af2_reuse", "Boltz-2": "boltz2", "Boltz-1": "boltz1"}[ARM]
    hits = sorted(glob.glob(f"{B}/cofold_arm/{d}/out/{target}/**/*.cif", recursive=True))
    if not hits:
        hits = sorted(f for f in glob.glob(f"{B}/cofold_arm/{d}/out/{target}/**/*.pdb", recursive=True)
                      if "unrelaxed" in os.path.basename(f))
    return hits


rows = []
skipped = []
n_ok = 0
for r in csv.DictReader(open(SUMMARY)):
    target, pep = r["pdb_id"], r["peptide_seq"]
    ref_path = f"{REFS}/{target}.cif"
    files = sample_files(target)
    if not os.path.exists(ref_path):
        skipped.append((target, "no reference cif")); continue
    if len(files) < 5:
        skipped.append((target, f"only {len(files)} samples")); continue
    try:
        exp = M.load_structure(ref_path)
        # reference peptide/MHC copy, chosen the same way the metric does
        probe = M.load_structure(files[0])
        probe_pep = M.find_peptide_chain_in_modeled(probe, pep)
        probe_mhc = M.find_mhc_chain_in_modeled(probe, probe_pep.id)
        copies = M.find_copies_in_reference(exp, pep, mhc_reference_seq=M.chain_sequence(probe_mhc))  # 3rd positional arg is peptide_length_range, NOT mhc_reference_seq
        if not copies:
            skipped.append((target, "no peptide/MHC copy in reference")); continue
        # Try each copy, not just the first, mirroring score_prediction() in the
        # metric (pmhc_rmsd.py ~line 519), which loops over copies and falls
        # through on ValueError. This is a consistency measure only: with
        # mhc_reference_seq passed correctly above, copies[0] is already the
        # right copy for every target in this benchmark. It is NOT what fixed
        # the nine previously-dropped targets -- that was the keyword argument.
        # Kept so the two code paths cannot disagree on which copy they score.
        exp_pep = exp_mhc = None
        first_err = None
        probe_mp = M.find_peptide_chain_in_modeled(probe, pep)
        for cand_pep, cand_mhc in copies:
            try:
                M.superpose_on_mhc(M.find_mhc_chain_in_modeled(probe, probe_mp.id), cand_mhc)
                exp_pep, exp_mhc = cand_pep, cand_mhc
                break
            except ValueError as e:
                first_err = first_err or e
        if exp_mhc is None:
            skipped.append((target, f"no copy superposed ({first_err})")); continue
        exp_res = [x for x in exp_pep if x.id[0] == " "]
        bf = {}
        for i, x in enumerate(exp_res):
            v = [a.get_bfactor() for a in x if a.element != "H" and a.get_bfactor() is not None]
            if v:
                bf[i] = float(np.mean(v))
        coords = collections.defaultdict(list)
        for f in files:
            s = M.load_structure(f)
            mp = M.find_peptide_chain_in_modeled(s, pep)
            mm = M.find_mhc_chain_in_modeled(s, mp.id)
            sup, _, _ = M.superpose_on_mhc(mm, exp_mhc)
            sup.apply([a for a in s.get_atoms()])
            for i, (mres, _eres) in enumerate(M._align_peptide_residues(mp, exp_pep)):
                if "CA" in mres:
                    coords[i].append(np.array(mres["CA"].coord, dtype=float))
        if not coords:
            skipped.append((target, "no CA coords recovered")); continue
        L = len(coords)
        for i, pts in sorted(coords.items()):
            if len(pts) < 5:
                continue
            P = np.vstack(pts)
            rmsf = float(np.sqrt(np.mean(np.sum((P - P.mean(axis=0)) ** 2, axis=1))))
            # Median PAIRWISE CA deviation over all n(n-1)/2 sample pairs.
            # Reported alongside RMSF because the sampling-parity claim is
            # stated in pairwise terms. Deliberately NOT the minimum pairwise
            # deviation (an extreme order statistic, near zero whenever any
            # two of 25 samples converge) and NOT max-minus-min (conflates
            # diversity with error). For an isotropic Gaussian cloud the mean
            # pairwise distance is sqrt(2)*RMSF, so the two are monotonically
            # related; both are reported so neither has to be defended alone.
            iu = np.triu_indices(len(P), k=1)
            pair_d = np.linalg.norm(P[iu[0]] - P[iu[1]], axis=1)
            rows.append(dict(arm=ARM, pdb_id=target, pos=i, L=L, n_samples=len(pts),
                             dist_from_terminus=min(i, L - 1 - i),
                             rmsf=rmsf,
                             median_pairwise=float(np.median(pair_d)),
                             mean_pairwise=float(np.mean(pair_d)),
                             ref_bfactor=bf.get(i)))
        n_ok += 1
        if n_ok % 25 == 0:
            print(f"  {n_ok} targets", flush=True)
    except Exception as e:
        # Was a bare "except Exception: continue". That silently dropped 9
        # targets from EVERY arm (9ASF 9ASG 9DL1 9FE1 9HKQ 9MIN 9NFC 9NNF
        # 9O5S) with no warning and no count, so the run reported 163/158
        # targets and looked complete. Never swallow silently here.
        skipped.append((target, f"{type(e).__name__}: {e}"))
        continue

out = f"{B}/positional/spread_{ARM.replace('-', '')}.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["arm", "pdb_id", "pos", "L", "n_samples",
                                       "dist_from_terminus", "rmsf",
                                       "median_pairwise", "mean_pairwise",
                                       "ref_bfactor"])
    w.writeheader(); w.writerows(rows)

print(f"\n{ARM}: {n_ok} targets, {len(rows)} residues -> {out}")
if skipped:
    print(f"  SKIPPED {len(skipped)} targets:")
    for t, why in skipped:
        print(f"    {t}: {why}")
if rows:
    d = [r["dist_from_terminus"] for r in rows]
    v = [r["rmsf"] for r in rows]
    rho, p = spearmanr(d, v)
    anc = [r["rmsf"] for r in rows if r["dist_from_terminus"] <= 1]
    cen = [r["rmsf"] for r in rows if r["dist_from_terminus"] >= 3]
    print(f"  rho(dist, ensemble spread) = {rho:.3f}  p = {p:.2e}")
    print(f"  mean RMSF  anchor {np.mean(anc):.3f} A   central {np.mean(cen):.3f} A")
    mp = [r["median_pairwise"] for r in rows]
    prho, pp = spearmanr(d, mp)
    print(f"  rho(dist, median pairwise CA) = {prho:.3f}  p = {pp:.2e}")
    print(f"  median pairwise CA  anchor "
          f"{np.mean([r['median_pairwise'] for r in rows if r['dist_from_terminus'] <= 1]):.3f} A"
          f"   central "
          f"{np.mean([r['median_pairwise'] for r in rows if r['dist_from_terminus'] >= 3]):.3f} A")
    bb = [(r["rmsf"], float(r["ref_bfactor"])) for r in rows if r.get("ref_bfactor") is not None]
    if len(bb) > 30:
        brho, bp = spearmanr([x[0] for x in bb], [x[1] for x in bb])
        print(f"  rho(ensemble spread, reference B-factor) = {brho:.3f}  p = {bp:.2e}   n={len(bb)}")
