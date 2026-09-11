"""Per-residue peptide error along the sequence, by method.

A single peptide-backbone RMSD can hide a systematic pattern: anchor residues
placed correctly while the central bulge is wrong. That distinction matters
because the interface, not the whole peptide, is what I_sc scores -- so a model
can look accurate globally while disagreeing with its own score label at the
same coordinates.

For each target's confidence-selected model this computes, after superposing on
the MHC exactly as the headline metric does:
  - per-residue backbone deviation
  - each residue's distance from the nearer terminus, min(i, L-1-i)
  - Spearman rho between the two, per method
  - anchor (P2 and P-omega) versus central-residue error
  - which peptide residues contact the MHC (interface), and their error
"""
import csv, os, re, sys, glob, collections, json
import numpy as np
sys.path.insert(0, "/data/p_csb_meiler/huntek1/benchmark/metric")
from scipy.stats import spearmanr
import pmhc_rmsd as M

B = "/data/p_csb_meiler/huntek1/benchmark"
REFS = f"{B}/refs/cif"
INTERFACE_CUTOFF = 5.0     # Angstrom, peptide heavy atom to MHC heavy atom

ARMS = {
    "Rosetta":    (f"{B}/rosetta_arm/results/arm_summary.csv",
                   lambda t, d: glob.glob(f"{B}/rosetta_arm/arm/{t}/docking/*_input_{int(d.split('_')[-1]):04d}.pdb")
                   if d and d.split('_')[-1].isdigit() else []),
    "AlphaFold2": (f"{B}/cofold_arm/af2_reuse/results/af2_reuse_summary.csv", None),
    "Boltz-2":    (f"{B}/cofold_arm/boltz2/results/boltz2_summary.csv", None),
    "Boltz-1":    (f"{B}/cofold_arm/boltz1/results/boltz1_summary.csv", None),
    "Chai-1":     (f"{B}/cofold_arm/chai1/results/chai1_summary.csv", None),
    "Protenix":   (f"{B}/cofold_arm/protenix/results/protenix_summary.csv", None),
}
DECOYS = {
    "Rosetta":    f"{B}/rosetta_arm/results/arm_decoys.csv",
    "AlphaFold2": f"{B}/cofold_arm/af2_reuse/results/af2_reuse_decoys.csv",
    "Boltz-2":    f"{B}/cofold_arm/boltz2/results/boltz2_decoys.csv",
    "Boltz-1":    f"{B}/cofold_arm/boltz1/results/boltz1_decoys.csv",
    "Chai-1":     f"{B}/cofold_arm/chai1/results/chai1_decoys.csv",
    "Protenix":   f"{B}/cofold_arm/protenix/results/protenix_decoys.csv",
}
SEARCH = {
    "Rosetta":    [f"{B}/rosetta_arm/arm/{{t}}/docking"],
    "AlphaFold2": [f"{B}/cofold_arm/af2_reuse/out/{{t}}"],
    "Boltz-2":    [f"{B}/cofold_arm/boltz2/out/{{t}}"],
    "Boltz-1":    [f"{B}/cofold_arm/boltz1/out/{{t}}"],
    "Chai-1":     [f"{B}/cofold_arm/chai1/out/{{t}}"],
    "Protenix":   [f"{B}/cofold_arm/protenix/out/{{t}}/pmhc_input/seed_101/predictions"],
}


def selected_file(arm, target, decoy_name):
    """Locate the structure file for a named selected decoy.

    Chai-1 names its samples "trunk_<T>_idx_<I>" in the results table but writes
    them to structures/trunk_<T>/pred.model_idx_<I>.cif, so a substring match
    never fires; it needs an explicit mapping.
    """
    if arm == "Chai-1":
        m = re.match(r"trunk_(\d+)_idx_(\d+)", decoy_name or "")
        if m:
            f = f"{B}/cofold_arm/chai1/out/{target}/structures/trunk_{m.group(1)}/pred.model_idx_{m.group(2)}.cif"
            if os.path.exists(f):
                return f
        hits = sorted(glob.glob(f"{B}/cofold_arm/chai1/out/{target}/structures/**/*.cif", recursive=True))
        return hits[0] if hits else None
    # Exact stem match, NOT a substring test. A substring test collides:
    # decoy "model_1" is contained in pmhc_input_model_10.cif ... _model_19.cif,
    # and likewise "sample_1" in pmhc_input_sample_10.cif, so the wrong
    # structure could be scored for decoys 1-9 in the Boltz and Protenix arms.
    # Every arm's filename is either exactly the decoy name (AlphaFold2) or the
    # decoy name preceded by an underscore (Boltz, Protenix).
    def _is_match(path):
        stem = os.path.splitext(os.path.basename(path))[0]
        return stem == decoy_name or stem.endswith(f"_{decoy_name}")

    if not decoy_name:
        return None
    for root in SEARCH[arm]:
        base = root.format(t=target)
        for ext in ("pdb", "cif"):
            hits = sorted(f for f in glob.glob(f"{base}/**/*.{ext}", recursive=True)
                          if _is_match(f))
            if hits:
                return hits[0]
    return None


def analyse(pred_path, ref_path, peptide):
    """Per-residue deviations after MHC superposition, plus interface flags."""
    mod = M.load_structure(pred_path)
    exp = M.load_structure(ref_path)
    mod_pep = M.find_peptide_chain_in_modeled(mod, peptide)
    mod_mhc = M.find_mhc_chain_in_modeled(mod, mod_pep.id)
    copies = M.find_copies_in_reference(exp, peptide, mhc_reference_seq=M.chain_sequence(mod_mhc))  # 3rd positional arg is peptide_length_range, NOT mhc_reference_seq
    if not copies:
        return None
    exp_pep, exp_mhc = copies[0][0], copies[0][1]
    sup, _, npc = M.superpose_on_mhc(mod_mhc, exp_mhc)
    sup.apply([a for a in mod.get_atoms()])
    pairs, _, _ = M._peptide_atom_pairs(mod_pep, exp_pep)
    per_res = collections.defaultdict(list)
    for (i, _atom), ma, ea in pairs:
        per_res[i].append(float(np.linalg.norm(ma.coord - ea.coord)))
    if not per_res:
        return None
    # interface: experimental peptide residues near the experimental MHC
    mhc_atoms = np.array([a.coord for a in exp_mhc.get_atoms() if a.element != "H"])
    exp_res = [r for r in exp_pep if r.id[0] == " "]
    iface = {}
    for i, r in enumerate(exp_res):
        ra = np.array([a.coord for a in r if a.element != "H"])
        if len(ra) and len(mhc_atoms):
            d = np.linalg.norm(ra[:, None, :] - mhc_atoms[None, :, :], axis=-1).min()
            iface[i] = bool(d <= INTERFACE_CUTOFF)
    # Reference crystal B-factors per peptide residue. If the central bulge is
    # genuinely mobile, the crystal itself is less certain there, and part of
    # the apparent "model error" is reference uncertainty rather than a wrong
    # prediction. This is what distinguishes the two explanations.
    bf = {}
    for i, r in enumerate(exp_res):
        vals = [a.get_bfactor() for a in r if a.element != "H" and a.get_bfactor() is not None]
        if vals:
            bf[i] = float(np.mean(vals))
    L = len(per_res)
    out = []
    for i, devs in sorted(per_res.items()):
        out.append(dict(pos=i, L=L,
                        dist_from_terminus=min(i, L - 1 - i),
                        rmsd=float(np.sqrt(np.mean(np.square(devs)))),
                        interface=iface.get(i),
                        ref_bfactor=bf.get(i)))
    return out


rows = []
for arm, (summ, _) in ARMS.items():
    if not os.path.exists(summ):
        print(f"  skip {arm}: no summary"); continue
    dec = {}
    if os.path.exists(DECOYS[arm]):
        for r in csv.DictReader(open(DECOYS[arm])):
            k = r["pdb_id"]
            v = r.get("peptide_backbone_rmsd")
            if v:
                dec.setdefault(k, []).append((float(v), r.get("decoy", "")))
    n_ok = 0
    for r in csv.DictReader(open(summ)):
        t, pep = r["pdb_id"], r["peptide_seq"]
        ref = f"{REFS}/{t}.cif"
        if not os.path.exists(ref):
            continue
        sel = r.get("selected")
        name = ""
        if t in dec and sel:
            try:
                name = min(dec[t], key=lambda x: abs(x[0] - float(sel)))[1]
            except Exception:
                name = ""
        f = selected_file(arm, t, name)
        if not f:
            continue
        try:
            res = analyse(f, ref, pep)
        except Exception:
            continue
        if not res:
            continue
        n_ok += 1
        for d in res:
            d.update(arm=arm, pdb_id=t)
            rows.append(d)
    print(f"  {arm}: {n_ok} targets analysed", flush=True)

out = f"{B}/positional/per_residue_error.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["arm", "pdb_id", "pos", "L",
                                       "dist_from_terminus", "rmsd", "interface", "ref_bfactor"])
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {out} ({len(rows)} residue observations)\n")

print(f"{'arm':12s} {'n_res':>6s} {'rho':>7s} {'p':>10s} {'anchor':>8s} {'central':>8s} {'iface':>8s} {'non-if':>8s}")
for arm in ARMS:
    sub = [r for r in rows if r["arm"] == arm]
    if len(sub) < 20:
        continue
    rho, p = spearmanr([r["dist_from_terminus"] for r in sub], [r["rmsd"] for r in sub])
    anchor = [r["rmsd"] for r in sub if r["dist_from_terminus"] <= 1]
    central = [r["rmsd"] for r in sub if r["dist_from_terminus"] >= 3]
    ifa = [r["rmsd"] for r in sub if r["interface"]]
    nif = [r["rmsd"] for r in sub if r["interface"] is False]
    f = lambda v: f"{np.mean(v):8.3f}" if v else "       -"
    print(f"{arm:12s} {len(sub):6d} {rho:7.3f} {p:10.2e} {f(anchor)} {f(central)} {f(ifa)} {f(nif)}")

print()
print("Reference-crystal B-factors: is the central bulge genuinely less well determined?")
bsub = [r for r in rows if r.get("ref_bfactor") not in (None, "")]
if bsub:
    seen, uniq = set(), []
    for r in bsub:                     # one observation per (target, residue)
        k = (r["pdb_id"], r["pos"])
        if k not in seen:
            seen.add(k); uniq.append(r)
    brho, bp = spearmanr([r["dist_from_terminus"] for r in uniq],
                         [float(r["ref_bfactor"]) for r in uniq])
    anc = [float(r["ref_bfactor"]) for r in uniq if r["dist_from_terminus"] <= 1]
    cen = [float(r["ref_bfactor"]) for r in uniq if r["dist_from_terminus"] >= 3]
    print(f"  residues: {len(uniq)}   rho(dist, B-factor) = {brho:.3f}  p = {bp:.2e}")
    print(f"  mean B-factor  anchor {np.mean(anc):.1f}   central {np.mean(cen):.1f}")
    print("  A positive rho means the crystal is itself less certain toward the peptide centre,")
    print("  so part of the apparent model error is reference uncertainty, not misprediction.")
else:
    print("  no B-factors recovered from the references")
