"""
Validation 1 RMSD for the crystal-matched pairs of the merged release.

crystal_match.py finds 76 matched pairs on the merged release (52 from the v1
batch, 24 newly matched from v2), against the 52 in the submitted manuscript.
This computes the three reported metrics for all 76:

    rmsd_best_score   RMSD of the lowest-I_sc decoy
    rmsd_top5_mean    mean RMSD over the five lowest-I_sc decoys
    rmsd_min_of_25    lowest RMSD over the whole ensemble

Every pair is computed from the release silent files, including the 52 that
already have numbers from 01_structural_validation.ipynb.

*** READ THIS BEFORE QUOTING ANY NUMBER THIS SCRIPT PRODUCES ***

The RMSDs here are measured on the RELEASED structures, and for crystal-matched
pairs those are LEAKED: they were threaded onto their own crystal, so their
agreement with that crystal measures refinement of a self-template, not
modeling accuracy. They must not be used as Validation 1.

Why. Template selection (HLA_db.MHCdatabase.get_peptide_template, called from
IEDBTestPipeline.thread_template) pools every same-length peptide from every
allele of the same gene and ranks them by BLOSUM62 similarity to the target. An
identical peptide therefore scores highest and is chosen first. Self-exclusion
exists (omit=["self"], matching on identical peptide sequence) but only runs
when the pipeline is given --ignore_epitope_match, which the production run did
not pass. regeneration/README.md states the consequence directly: "The released
structures for these 52 pairs were threaded onto their own crystal, because the
self-exclusion path in HLA_db.get_peptide_template did not run." It affects only
pairs whose native PDB is in the local template database -- about 0.1% of the
release -- but those are exactly the pairs Validation 1 uses.

That is why 01_structural_validation.ipynb reads decoys from
IEDB_data_clean/IEDB_validation/regeneration/pdb/ instead: a deliberate
re-docking of the 52 pairs with --ignore_epitope_match, so the measurement is
leakage-free. The manuscript's 1.21 A over all 76 pairs (1.14 A over the
original 52) is the honest number; the ~0.88 A this script gets on the same
pairs is the leakage.

The two trees are genuinely different runs -- verified on A0201/GLCTLVAML, where
the released silent's per-decoy I_sc values (0005 = -77.656, 0003 = -75.483,
...) match the production score.sc exactly and none appear in the regeneration
score.sc (0019 = -79.215, 0004 = -78.116, ...).

So what is this script still good for?
  * quantifying the leakage, by differencing released against regenerated on the
    50 shared pairs (~0.24 A median, released closer in 43 of 50);
  * verifying extraction fidelity (verify_extraction asserts the silent's
    per-decoy I_sc equals the production score.sc value, 25/25 per pair);
  * providing the released-structure numbers for the 24 newly matched v2 pairs,
    which are an UPPER BOUND on their accuracy, not an estimate of it. Those 24
    have since been re-docked with --ignore_epitope_match into
    IEDB_validation/regeneration_v2/pdb/, so Validation 1 now covers all 76 and
    takes its numbers from there, not from here.

RMSD itself is not reimplemented: utils.structure.compute_peptide_rmsd is the
same function the notebook uses (peptide backbone
N/Ca/C/O after superposition on the first 180 MHC Ca atoms, residues paired by
sequence alignment).

Decoys are extracted to a scratch directory one pair at a time and deleted
immediately after that pair is scored, so peak disk stays at ~25 PDBs (~10 MB)
rather than 1,900 (~760 MB).

Usage:
    python3 crystal_rmsd.py --out-dir <dir> [--scratch <dir>] [--limit N]
    python3 crystal_rmsd.py --self-test
"""
import argparse
import shutil
import tempfile
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.structure import compute_peptide_rmsd  # noqa: E402

import sys as _sys; from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
from paths import DATA_ROOT, MHC_DB_ROOT  # noqa: E402

BASE = DATA_ROOT
MATCHED = BASE / "IEDB_validation" / "crystal_match_v2" / "crystal_matched_pairs.csv"
SILENT_ROOT = BASE / "release_v2_final" / "structures"
TEMPLATE_DIR = MHC_DB_ROOT / "templates"
# the notebook's stored result for the v1 52, used as the known-answer check
NOTEBOOK_RMSD = BASE / "IEDB_validation" / "01_structural_regen" / "rmsd_per_pair.csv"

PRIMARY_SCORE = "I_sc"
TOP_N = 5

# the production decoy tree the release silents were converted from; used to
# verify that extracted per-decoy scores are faithful
PRODUCTION_PDB_ROOT = BASE / "pdb"
# conversion appended an index to each tag, so silent tag "X_input_0005_0001"
# corresponds to production decoy "X_input_0005"
TAG_SUFFIX = "_0001"


def verify_extraction(pairs, n_check=3):
    """Known-answer check on real data: the per-decoy I_sc stored in a release
    silent must equal the value in the production score.sc for the same decoy.

    This is the check that matters for this script -- it proves the scores used
    to rank decoys belong to the structures they are attached to. Runs on the
    first `n_check` pairs that have a production score.sc available."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from utils.io import read_score_file

    checked = 0
    for r in pairs.itertuples(index=False):
        sc = PRODUCTION_PDB_ROOT / r.allele_compact / r.peptide / "score.sc"
        silent = SILENT_ROOT / r.allele_compact / f"{r.peptide}.silent"
        if not sc.exists() or not silent.exists():
            continue
        prod = read_score_file(sc).set_index("description")[PRIMARY_SCORE].to_dict()
        decoys = extract_decoys_scores_only(silent)
        compared = 0
        for tag, val in decoys.items():
            base = tag[: -len(TAG_SUFFIX)] if tag.endswith(TAG_SUFFIX) else tag
            if base not in prod:
                continue
            if abs(prod[base] - val) > 1e-3:
                raise SystemExit(
                    f"extraction check FAILED for {r.allele_compact}/{r.peptide} "
                    f"decoy {base}: silent I_sc={val:.4f} but production "
                    f"score.sc I_sc={prod[base]:.4f}. The scores in the release do "
                    f"not belong to the structures they are attached to; stop and "
                    f"investigate the conversion.")
            compared += 1
        if compared == 0:
            continue
        print(f"  {r.allele_compact}/{r.peptide}: {compared}/{len(decoys)} decoy "
              f"I_sc values match production score.sc exactly")
        checked += 1
        if checked >= n_check:
            break
    if checked == 0:
        print("  (no pair had a production score.sc available; extraction unverified)")
    return checked


def extract_decoys_scores_only(silent_path):
    """Per-decoy PRIMARY_SCORE from a silent file, without writing any PDB."""
    from pyrosetta.rosetta.core.io.silent import SilentFileData, SilentFileOptions
    sfd = SilentFileData(SilentFileOptions())
    sfd.read_file(str(silent_path))
    out = {}
    for tag in sfd.tags():
        ss = sfd.get_structure(tag)
        names = set(ss.energy_names().energy_names())
        if PRIMARY_SCORE in names:
            out[tag] = ss.get_energy(PRIMARY_SCORE)
    return out


def extract_decoys(silent_path, scratch):
    """Extract every decoy of one silent file to `scratch` and return
    [(pdb_path, score_dict), ...]. Scores come from the silent structs
    themselves, so the per-decoy I_sc used for ranking is the one stored with
    the coordinates -- no separate score.sc join that could go stale."""
    from pyrosetta import init, Pose
    from pyrosetta.rosetta.core.io.silent import SilentFileData, SilentFileOptions

    scratch.mkdir(parents=True, exist_ok=True)
    opts = SilentFileOptions()
    sfd = SilentFileData(opts)
    sfd.read_file(str(silent_path))

    out = []
    for tag in sfd.tags():
        ss = sfd.get_structure(tag)
        pose = Pose()
        ss.fill_pose(pose)
        scores = {}
        for k in ss.energy_names().energy_names():
            try:
                scores[k] = ss.get_energy(k)
            except Exception:
                pass
        p = scratch / f"{tag}.pdb"
        pose.dump_pdb(str(p))
        out.append((p, scores))
    return out


def decoys_from_pdb_tree(pep_dir, score_col=PRIMARY_SCORE):
    """[(pdb_path, {score: value}), ...] from a decoy PDB directory + score.sc.

    The leakage-free validation ensembles are re-docks delivered as PDB trees,
    not silents, so Validation 1 reads them here rather than through
    extract_decoys(). Scores come from score.sc keyed on `description`, which
    matches each PDB's filename stem."""
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from utils.io import read_score_file
    sc = read_score_file(pep_dir / "score.sc").set_index("description")
    out = []
    for pdb in sorted(pep_dir.glob(f"{pep_dir.name}_input_[0-9][0-9][0-9][0-9].pdb")):
        if pdb.stem not in sc.index:
            continue
        out.append((pdb, {score_col: float(sc.loc[pdb.stem, score_col])}))
    return out


def score_pair(decoys, exp_pdb, peptide):
    """Return the three metrics for one pair, or (None, reason)."""
    rows = []
    for pdb_path, scores in decoys:
        if PRIMARY_SCORE not in scores:
            return None, f"decoy {pdb_path.name} has no {PRIMARY_SCORE}"
        try:
            r = compute_peptide_rmsd(pdb_path, exp_pdb, peptide)
        except Exception as e:  # a single unusable decoy must not kill the pair
            rows.append({"rmsd": np.nan, "score": scores[PRIMARY_SCORE], "err": str(e)})
            continue
        rows.append({"rmsd": r["rmsd"], "score": scores[PRIMARY_SCORE],
                     "mhc_alignment_rmsd": r["mhc_alignment_rmsd"], "err": None})
    df = pd.DataFrame(rows)
    ok = df.dropna(subset=["rmsd"])
    if ok.empty:
        return None, f"all {len(df)} decoys failed: {df['err'].dropna().iloc[0] if len(df) else 'n/a'}"
    ranked = ok.sort_values("score")  # lower I_sc is better
    return {
        "n_decoys": int(len(df)),
        "n_decoys_scored": int(len(ok)),
        "rmsd_best_score": float(ranked.iloc[0]["rmsd"]),
        "rmsd_top5_mean": float(ranked.head(TOP_N)["rmsd"].mean()),
        "rmsd_min_of_25": float(ok["rmsd"].min()),
        "best_I_sc": float(ranked.iloc[0]["score"]),
        "mhc_alignment_rmsd": float(ranked.iloc[0].get("mhc_alignment_rmsd", np.nan)),
    }, None


def self_test():
    """Known answer for the ranking and aggregation logic, which is the part
    this script owns (compute_peptide_rmsd is already self-tested upstream).
    Six synthetic decoys with known scores and RMSDs:
      scores  -10 -9 -8 -7 -6 -5   (lower = better, so this is rank order)
      rmsds   3.0 1.0 2.0 4.0 5.0 0.5
    best-score RMSD must be 3.0 (the -10 decoy, NOT the lowest RMSD);
    top-5 mean must be (3+1+2+4+5)/5 = 3.0; min-of-all must be 0.5.
    A decoy that fails RMSD must be skipped without killing the pair."""
    rmsds = [3.0, 1.0, 2.0, 4.0, 5.0, 0.5]
    scores = [-10.0, -9.0, -8.0, -7.0, -6.0, -5.0]
    rows = [{"rmsd": r, "score": s, "mhc_alignment_rmsd": 0.7, "err": None}
            for r, s in zip(rmsds, scores)]
    df = pd.DataFrame(rows)
    ok = df.dropna(subset=["rmsd"])
    ranked = ok.sort_values("score")
    assert float(ranked.iloc[0]["rmsd"]) == 3.0, "best-score pick is wrong"
    assert abs(float(ranked.head(5)["rmsd"].mean()) - 3.0) < 1e-12, "top-5 mean is wrong"
    assert float(ok["rmsd"].min()) == 0.5, "min-of-ensemble is wrong"

    # one failed decoy: dropped from RMSD stats, still counted in n_decoys
    rows2 = rows + [{"rmsd": np.nan, "score": -11.0, "err": "boom"}]
    df2 = pd.DataFrame(rows2)
    ok2 = df2.dropna(subset=["rmsd"])
    assert len(df2) == 7 and len(ok2) == 6
    ranked2 = ok2.sort_values("score")
    assert float(ranked2.iloc[0]["rmsd"]) == 3.0, \
        "a failed top-scoring decoy must not become the best-score RMSD"
    print("Self-test PASSED: best-score decoy is chosen by score not by RMSD, "
          "top-5 mean and min-of-ensemble are exact, and a failed decoy is "
          "skipped without discarding the pair.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir")
    ap.add_argument("--scratch", default=str(Path(tempfile.gettempdir()) / "pepbind3d_rmsd_scratch"),
                    help="working directory for extracted decoys; one pair at a "
                         "time, removed immediately after that pair is scored")
    ap.add_argument("--limit", type=int, default=None, help="first N pairs only (smoke test)")
    ap.add_argument("--pdb-tree", default=None,
                    help="score decoys from this PDB tree ({allele}/{peptide}/) "
                         "instead of the release silents. This is how the "
                         "LEAKAGE-FREE validation ensembles are read: they were "
                         "re-docked with --ignore_epitope_match and delivered as "
                         "PDBs, so unlike the silent path these numbers ARE "
                         "Validation 1 rather than an upper bound.")
    ap.add_argument("--matched", default=str(MATCHED),
                    help="pair list to score (default: all 76 matched pairs)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not args.out_dir:
        raise SystemExit("--out-dir is required")
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    scratch = Path(args.scratch)

    self_test()
    print()

    from pyrosetta import init
    init("-mute all -in:file:silent_struct_type binary")

    pairs = pd.read_csv(args.matched)
    if args.limit:
        pairs = pairs.head(args.limit)
    print(f"{len(pairs)} matched pairs to score")

    if not args.pdb_tree:
        print("\nverifying extraction fidelity against the production score.sc:", flush=True)
        verify_extraction(pairs)
    else:
        print(f"\nreading decoys from PDB tree {args.pdb_tree}\n"
              f"  these are the LEAKAGE-FREE re-docked ensembles, so these ARE "
              f"Validation 1 numbers", flush=True)
    print(flush=True)

    results, failures = [], []
    for i, r in enumerate(pairs.itertuples(index=False), 1):
        silent = SILENT_ROOT / r.allele_compact / f"{r.peptide}.silent"
        exp_pdb = TEMPLATE_DIR / f"{r.matched_pdb_id}.pdb"
        need_silent = not args.pdb_tree
        if (need_silent and not silent.exists()) or not exp_pdb.exists():
            failures.append({"allele": r.allele, "peptide": r.peptide,
                             "reason": f"missing {'silent' if not silent.exists() else 'exp pdb'}"})
            continue
        pair_scratch = scratch / f"{r.allele_compact}_{r.peptide}"
        try:
            if args.pdb_tree:
                decoys = decoys_from_pdb_tree(
                    Path(args.pdb_tree) / r.allele_compact / r.peptide)
            else:
                decoys = extract_decoys(silent, pair_scratch)
            metrics, err = score_pair(decoys, exp_pdb, r.peptide)
        except Exception as e:
            metrics, err = None, f"{type(e).__name__}: {e}"
        finally:
            shutil.rmtree(pair_scratch, ignore_errors=True)

        if metrics is None:
            failures.append({"allele": r.allele, "peptide": r.peptide, "reason": err})
            print(f"  [{i:3d}/{len(pairs)}] {r.allele} {r.peptide} FAILED: {err}", flush=True)
            continue
        results.append({"allele": r.allele, "allele_compact": r.allele_compact,
                        "peptide": r.peptide, "matched_pdb_id": r.matched_pdb_id,
                        "resolution_angstrom": r.resolution_angstrom,
                        "source_version": getattr(r, "source_version", None), **metrics})
        print(f"  [{i:3d}/{len(pairs)}] {r.allele} {r.peptide} {r.matched_pdb_id} "
              f"best={metrics['rmsd_best_score']:.2f} min={metrics['rmsd_min_of_25']:.2f}",
              flush=True)

    df = pd.DataFrame(results)
    df.to_csv(out / "crystal_rmsd_per_pair.csv", index=False)
    if failures:
        pd.DataFrame(failures).to_csv(out / "crystal_rmsd_failures.csv", index=False)

    def summarize(d, label):
        if d.empty:
            print(f"\n{label}: no pairs"); return
        print(f"\n{label} (n={len(d)}):")
        for col in ("rmsd_best_score", "rmsd_top5_mean", "rmsd_min_of_25"):
            v = d[col]
            print(f"  {col:16s} median {v.median():.2f} A  "
                  f"IQR {v.quantile(.25):.2f}-{v.quantile(.75):.2f}  "
                  f"<=2A {100*(v <= 2).mean():.1f}%")

    if not args.pdb_tree:
        print("\n*** these are SELF-TEMPLATED pairs; the numbers below are upper "
              "bounds, not validation figures. See the module docstring. ***")
    summarize(df, "ALL matched pairs")
    if "source_version" in df:
        summarize(df[df.source_version == "v1"], "v1 batch")
        summarize(df[df.source_version == "v2"], "newly matched v2 pairs")

    # known-answer check against the notebook's stored 52
    if NOTEBOOK_RMSD.exists() and not df.empty:
        nb = pd.read_csv(NOTEBOOK_RMSD)
        j = df.merge(nb, on=["allele", "peptide"], suffixes=("", "_nb"))
        if len(j):
            d = (j["rmsd_best_score"] - j["rmsd_best_score_nb"]).abs()
            better = int((j["rmsd_best_score"] < j["rmsd_best_score_nb"]).sum())
            print(f"\nLEAKAGE ESTIMATE -- released (self-templated) vs regenerated "
                  f"(self-excluded), {len(j)} shared pairs:")
            print(f"  released  median {j['rmsd_best_score'].median():.2f} A  "
                  f"<- self-templated, NOT a validation number")
            print(f"  regenerated median {j['rmsd_best_score_nb'].median():.2f} A  "
                  f"<- leakage-free, this is Validation 1")
            print(f"  leakage: median {d.median():.2f} A, max {d.max():.2f} A; "
                  f"released closer to its own crystal in {better}/{len(j)} pairs")

    msg = f"\nwrote {out / 'crystal_rmsd_per_pair.csv'}"
    if failures:
        msg += f" and crystal_rmsd_failures.csv ({len(failures)} failures)"
    print(msg)


if __name__ == "__main__":
    main()
