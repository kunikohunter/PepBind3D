"""
Phase 1 item 1.2 / Reviewer 2's ensemble-scatter question.
R2 asked, verbatim: "How distinct are the 25 decoys from each other? ... In the
validations, where ensembles are compared to crystal structures, I wonder if
the RMSD values fall into the same range as e.g. pairwise RMSD between the
ensemble decoys."

So for each validation pair we need two peptide-backbone RMSD distributions,
computed with the SAME metric (superpose on MHC, RMSD on peptide N/Cα/C/O):
  - decoy-to-decoy  : all 25*24/2 = 300 pairwise RMSDs  -> the ensemble's own scatter
  - decoy-to-crystal: all 25 decoy-vs-crystal RMSDs      -> how far the ensemble sits from truth

"Crystal inside the scatter" means the ensemble's displacement from the crystal
is comparable to (or smaller than) its internal spread. We report, per pair,
median/mean/max of each distribution and the ratio median_d2c / median_d2d,
then summarize across pairs.

RMSD is computed by reusing utils.structure.superpose_on_mhc and
_peptide_atom_pairs verbatim (the same functions notebook 01
uses) -- no second RMSD implementation. Structures are loaded
once per pair and transformed non-destructively (the superposition rotran is
applied to a numpy copy of the peptide coords, never to the shared structure),
so a decoy reused as mobile in 24 comparisons is never corrupted.

Usage:
    $ python3 ensemble_diversity.py --out-dir <dir>
Self-test:
    $ python3 ensemble_diversity.py --self-test
"""
import argparse
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.structure import (  # noqa: E402

    load_structure,
    superpose_on_mhc,
    identify_peptide_chain,
    identify_mhc_chain,
    chain_sequence,
    _peptide_atom_pairs,
)

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT, MHC_DB_ROOT  # noqa: E402


RMSD_PER_PAIR = DATA_ROOT / "IEDB_validation/01_structural_regen/rmsd_per_pair.csv"
# Both re-docked trees. rmsd_per_pair.csv carries all 76 crystal-matched pairs;
# reading only the first tree leaves the 24 with no decoys, so they fail the
# "< 2 decoys" check, are skipped to stderr, and the run silently reports the
# old 52-pair result.
REGEN_ROOTS = [DATA_ROOT / "IEDB_validation/regeneration/pdb",
               DATA_ROOT / "IEDB_validation/regeneration_v2/pdb"]
REGEN_ROOT = REGEN_ROOTS[0]   # self-test substrate only


def resolve_regen_root(allele_dir, peptide):
    """Return the re-docked tree holding this pair."""
    for root in REGEN_ROOTS:
        if (root / allele_dir / peptide).is_dir():
            return root
    raise FileNotFoundError(
        f"{allele_dir}/{peptide} is in neither re-docked tree")
TEMPLATE_DIR = MHC_DB_ROOT / "templates"
CACHE_DIR = MHC_DB_ROOT / "pdb_cache"


def allele_to_dir(a: str) -> str:
    s = a[4:] if a.startswith("HLA-") else a
    return s.replace("*", "").replace(":", "")


def _apply_rotran(coords, rotran):
    """Reproduce Bio.PDB.Superimposer.apply on a numpy coord array without
    mutating any atom: moved = coord @ rot + tran (BioPython convention)."""
    rot, tran = rotran
    return np.dot(coords, rot) + tran


def get_pep_mhc_chains(struct, pep_seq):
    """Return (peptide_chain, mhc_chain) for a modeled decoy. Tries the strict
    A=MHC/B=peptide convention first (regeneration set), falls back to
    sequence/length-based identification (robust to varying chain
    IDs, which are not guaranteed to be A/B)."""
    model = struct[0]
    chain_ids = [c.id for c in model]
    if "A" in chain_ids and "B" in chain_ids and chain_sequence(model["B"]) == pep_seq:
        return model["B"], model["A"]  # B really is the peptide
    pep = identify_peptide_chain(struct, pep_seq)
    mhc = identify_mhc_chain(struct, exclude_chain_id=pep.id)
    return pep, mhc


def peptide_rmsd_between(mob_pep, mob_mhc, ref_pep, ref_mhc):
    """Peptide-backbone RMSD after superposing the mobile MHC onto the
    reference MHC. Uses superpose_on_mhc (reused) for the transform and
    _peptide_atom_pairs (reused) for residue-aligned N/Cα/C/O pairing.
    Non-destructive: mobile atoms are never moved in place."""
    sup, mhc_rmsd = superpose_on_mhc(mob_mhc, ref_mhc)
    pairs = _peptide_atom_pairs(mob_pep, ref_pep)[0]
    if not pairs:
        return None, mhc_rmsd
    mob_coords = np.array([mob_atom.coord for mob_atom, _ in pairs])
    ref_coords = np.array([ref_atom.coord for _, ref_atom in pairs])
    moved = _apply_rotran(mob_coords, sup.rotran)
    rmsd = float(np.sqrt(np.mean(np.sum((moved - ref_coords) ** 2, axis=1))))
    return rmsd, mhc_rmsd


def analyze_pair(allele, peptide, matched_pdb_id):
    """Return (decoy_to_decoy_rmsds[list], decoy_to_crystal_rmsds[list]) for
    one validation pair, or (None, None, reason) if it can't be computed."""
    allele_dir = allele_to_dir(allele)
    try:
        pep_dir = resolve_regen_root(allele_dir, peptide) / allele_dir / peptide
    except FileNotFoundError as e:
        return None, None, str(e)
    decoy_paths = sorted(pep_dir.glob(f"{peptide}_input_[0-9][0-9][0-9][0-9].pdb"))
    if len(decoy_paths) < 2:
        return None, None, f"only {len(decoy_paths)} decoys"

    # Load each decoy once; chain A = MHC, chain B = peptide (project convention).
    decoys = []
    for p in decoy_paths:
        s = load_structure(p)
        decoys.append((s[0]["B"], s[0]["A"]))  # (peptide, mhc)

    # decoy-to-decoy: all unique pairs, non-destructive.
    d2d = []
    for i in range(len(decoys)):
        for j in range(i + 1, len(decoys)):
            r, _ = peptide_rmsd_between(decoys[i][0], decoys[i][1],
                                        decoys[j][0], decoys[j][1])
            if r is not None:
                d2d.append(r)

    # decoy-to-crystal.
    d2c = []
    from utils.structure import get_experimental_pdb
    try:
        crystal_path = get_experimental_pdb(matched_pdb_id, TEMPLATE_DIR, CACHE_DIR)
        cstruct = load_structure(crystal_path)
        cpep = identify_peptide_chain(cstruct, peptide)
        cmhc = identify_mhc_chain(cstruct, exclude_chain_id=cpep.id)
        for pep_ch, mhc_ch in decoys:
            r, _ = peptide_rmsd_between(pep_ch, mhc_ch, cpep, cmhc)
            if r is not None:
                d2c.append(r)
    except Exception as e:  # noqa: BLE001
        return d2d, None, f"crystal RMSD failed: {e}"

    return d2d, d2c, None




def self_test():
    """Known-answer test: translate a real decoy's peptide by a fixed vector
    and confirm the measured peptide RMSD equals the translation magnitude.
    Because superposition is on the (untranslated) MHC, a rigid shift of only
    the peptide is NOT removed and must appear in full as the RMSD."""
    # find any validation decoy to use as substrate
    sample = next((p for p in REGEN_ROOT.glob("*/*/*_input_0001.pdb")), None)
    if sample is None:
        raise SystemExit("self-test: no decoy PDB found under REGEN_ROOT")

    s = load_structure(sample)
    mhc, pep = s[0]["A"], s[0]["B"]

    # identical structure vs itself -> RMSD must be 0
    s2 = deepcopy(s)
    r0, _ = peptide_rmsd_between(pep, mhc, s2[0]["B"], s2[0]["A"])
    assert r0 < 1e-6, f"self-test FAILED: identical structures gave RMSD {r0}, expected 0"

    # translate the peptide of a copy by exactly (1.0, 0, 0) Å; MHC untouched.
    s3 = deepcopy(s)
    shift = np.array([1.0, 0.0, 0.0])
    for atom in s3[0]["B"].get_atoms():
        atom.coord = atom.coord + shift
    r1, _ = peptide_rmsd_between(s3[0]["B"], s3[0]["A"], pep, mhc)
    assert abs(r1 - 1.0) < 1e-4, (
        f"self-test FAILED: peptide translated by 1.0 Å gave RMSD {r1}, expected 1.0 "
        f"(rot/tran convention or atom pairing is wrong)"
    )

    print(f"Self-test PASSED: identical->{r0:.2e} Å (≈0); "
          f"peptide shifted 1.0 Å->{r1:.4f} Å (≈1.0), so superposition + "
          f"non-destructive transform + atom pairing are correct.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=str, default=None)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--limit", type=int, default=None, help="debug: only first N pairs")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if args.out_dir is None:
        raise SystemExit("--out-dir required")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)


    pairs = pd.read_csv(RMSD_PER_PAIR)
    if args.limit:
        pairs = pairs.head(args.limit)

    rows = []
    for _, row in pairs.iterrows():
        allele, peptide = row["allele"], row["peptide"]
        matched = row["matched_pdb_id"]
        d2d, d2c, reason = analyze_pair(allele, peptide, matched)
        if not d2d:
            print(f"  SKIP {allele}/{peptide}: {reason}", file=sys.stderr)
            continue
        d2d = np.array(d2d)
        rec = {
            "allele": allele, "peptide": peptide, "matched_pdb_id": matched,
            "template_identity": row.get("template_identity"),
            "n_d2d": len(d2d),
            "d2d_median": float(np.median(d2d)),
            "d2d_mean": float(d2d.mean()),
            "d2d_max": float(d2d.max()),
        }
        has_crystal = bool(d2c)
        if has_crystal:
            d2c = np.array(d2c)
            rec.update({
                "n_d2c": len(d2c),
                "d2c_median": float(np.median(d2c)),
                "d2c_mean": float(d2c.mean()),
                "d2c_min": float(d2c.min()),
                "ratio_med_d2c_over_d2d": float(np.median(d2c) / np.median(d2d)),
                "crystal_within_d2d_range": bool(np.median(d2c) <= d2d.max()),
            })
        rows.append(rec)
        tag = f"ratio={rec.get('ratio_med_d2c_over_d2d', float('nan')):.2f}" if has_crystal else "(no crystal)"
        print(f"  {allele}/{peptide}: d2d_med={rec['d2d_median']:.2f} "
              f"d2c_med={rec.get('d2c_median', float('nan')):.2f} {tag}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "ensemble_diversity_per_pair.csv", index=False)

    with_crystal = df.dropna(subset=["ratio_med_d2c_over_d2d"]) if "ratio_med_d2c_over_d2d" in df else df.iloc[:0]
    print("\n=== SUMMARY ===")
    print(f"pairs analyzed: {len(df)}  (with crystal RMSD: {len(with_crystal)})")
    print(f"decoy-to-decoy median RMSD across pairs: "
          f"{df['d2d_median'].median():.3f} Å (IQR {df['d2d_median'].quantile(.25):.3f}-{df['d2d_median'].quantile(.75):.3f})")
    if len(with_crystal):
        n_inside = int(with_crystal["crystal_within_d2d_range"].sum())
        print(f"decoy-to-crystal median RMSD across pairs: "
              f"{with_crystal['d2c_median'].median():.3f} Å")
        print(f"median(ratio d2c/d2d) across pairs: {with_crystal['ratio_med_d2c_over_d2d'].median():.2f}")
        print(f"crystal falls within the ensemble's own decoy-decoy range "
              f"in {n_inside}/{len(with_crystal)} pairs")
    print(f"\nWrote {out_dir/'ensemble_diversity_per_pair.csv'}")


if __name__ == "__main__":
    main()
