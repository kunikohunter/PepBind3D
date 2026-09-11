#!/data/p_csb_meiler/huntek1/envs/ensemble/bin/python
"""Score one target for one co-folding arm against the crystal reference.

Usage:
    score_target.py <arm> <pdb_id>

<arm> in {af2_reuse, boltz2, boltz1, chai1, protenix}

Writes:
    cofold_arm/<arm>/results/tmp/<pdb_id>.summary.csv
    cofold_arm/<arm>/results/tmp/<pdb_id>.decoys.csv

Mirrors the output shape of rosetta_arm/results/arm_summary.csv (see
rosetta_arm/score_arm.py for the reference behavior this reproduces).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

BENCH_ROOT = Path("/data/p_csb_meiler/huntek1/benchmark")
COFOLD_ROOT = BENCH_ROOT / "cofold_arm"
REFS_DIR = BENCH_ROOT / "refs" / "cif"
JOBLIST = COFOLD_ROOT / "joblist.tsv"

sys.path.insert(0, str(BENCH_ROOT / "metric"))
from pmhc_rmsd import score_prediction  # noqa: E402

ARMS = ("af2_reuse", "boltz2", "boltz1", "chai1", "protenix")

SUMMARY_HEADER = [
    "pdb_id", "peptide_seq", "n_scored", "selected", "oracle", "median",
    "rank_of_selected", "percentile_of_selected",
]
DECOY_HEADER = [
    "pdb_id", "peptide_seq", "decoy", "peptide_backbone_rmsd", "confidence",
    "n_atom_pairs", "n_atoms_dropped", "n_paired_ca", "superposition_rmsd",
    "warnings",
]


def load_peptide_seq(pdb_id: str):
    with open(JOBLIST, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["pdb_id"] == pdb_id:
                return row["peptide_seq"]
    return None


def enumerate_samples_af2_reuse(pdb_id: str):
    """Return list of (decoy_name, pred_path, confidence, is_selected).

    Uses the canonical unrelaxed_model_{i}_multimer_v3_pred_{j}.pdb files
    named in MANIFEST_unrelaxed.txt, NOT ranked_*.pdb (rank-reordered, and
    ranked_0 is the lone Amber-relaxed structure -- see MANIFEST_unrelaxed.txt).
    """
    base = COFOLD_ROOT / "af2_reuse" / "out" / pdb_id / "pmhc_input"
    ranking_path = base / "ranking_debug.json"
    manifest_path = base / "MANIFEST_unrelaxed.txt"
    if not ranking_path.exists() or not manifest_path.exists():
        return None
    with open(ranking_path) as f:
        ranking = json.load(f)
    order = ranking["order"]
    iptm_ptm = ranking["iptm+ptm"]

    with open(manifest_path) as f:
        filenames = [
            line.strip() for line in f
            if line.strip() and not line.strip().startswith("#")
        ]

    samples = []
    for fname in filenames:
        pred_path = base / fname
        if not pred_path.exists():
            continue
        # tag = filename minus "unrelaxed_" prefix and ".pdb" suffix
        tag = fname[len("unrelaxed_"):-len(".pdb")]
        conf = iptm_ptm.get(tag)
        samples.append({
            "decoy": fname[:-len(".pdb")],
            "tag": tag,
            "path": pred_path,
            "confidence": conf,
        })

    selected_tag = order[0] if order else None
    selected_decoy = None
    for s in samples:
        if s["tag"] == selected_tag:
            selected_decoy = s["decoy"]
            break

    if samples:
        # sanity check: the AF2 top-ranked (order[0]) sample's confidence
        # should be the max among the unrelaxed samples' confidences.
        confs = [(s["confidence"], s["decoy"]) for s in samples if s["confidence"] is not None]
        if confs:
            max_conf, max_decoy = max(confs, key=lambda x: x[0])
            if max_decoy != selected_decoy:
                print(
                    f"WARNING {pdb_id}: af2_reuse selected decoy {selected_decoy!r} "
                    f"is not the max iptm+ptm sample (max is {max_decoy}={max_conf})",
                    file=sys.stderr,
                )
    return samples, selected_decoy


def enumerate_samples_boltz(arm: str, pdb_id: str):
    base = (
        COFOLD_ROOT / arm / "out" / pdb_id / "boltz_results_pmhc_input"
        / "predictions" / "pmhc_input"
    )
    if not base.exists():
        return None
    samples = []
    for n in range(25):
        pred_path = base / f"pmhc_input_model_{n}.cif"
        conf_path = base / f"confidence_pmhc_input_model_{n}.json"
        if not pred_path.exists():
            continue
        conf = None
        if conf_path.exists():
            try:
                with open(conf_path) as f:
                    conf = json.load(f).get("confidence_score")
            except Exception as e:
                print(f"WARNING {pdb_id}: could not read {conf_path}: {e}", file=sys.stderr)
        samples.append({
            "decoy": f"model_{n}",
            "path": pred_path,
            "confidence": conf,
        })
    if not samples:
        return None
    # selected = max confidence_score, ties -> lowest index
    scored = [s for s in samples if s["confidence"] is not None]
    if not scored:
        return samples, None
    best = max(scored, key=lambda s: (s["confidence"], -int(s["decoy"].split("_")[1])))
    return samples, best["decoy"]


def enumerate_samples_chai1(pdb_id: str):
    base = COFOLD_ROOT / "chai1" / "out" / pdb_id / "structures"
    if not base.exists():
        return None
    samples = []
    for trunk in range(5):
        trunk_dir = base / f"trunk_{trunk}"
        for idx in range(5):
            pred_path = trunk_dir / f"pred.model_idx_{idx}.cif"
            score_path = trunk_dir / f"scores.model_idx_{idx}.npz"
            if not pred_path.exists():
                continue
            conf = None
            if score_path.exists():
                try:
                    npz = np.load(score_path)
                    conf = float(npz["aggregate_score"][0])
                except Exception as e:
                    print(f"WARNING {pdb_id}: could not read {score_path}: {e}", file=sys.stderr)
            samples.append({
                "decoy": f"trunk_{trunk}_idx_{idx}",
                "path": pred_path,
                "confidence": conf,
            })
    if not samples:
        return None
    scored = [s for s in samples if s["confidence"] is not None]
    if not scored:
        return samples, None
    best = max(scored, key=lambda s: s["confidence"])
    return samples, best["decoy"]


def enumerate_samples_protenix(pdb_id: str):
    """Protenix writes one directory per seed; the arm uses a single seed (101)
    with 25 diffusion samples, matching Boltz-2's 25 diffusion samples.

    Confidence is `ranking_score` from the per-sample summary_confidence JSON.
    That is Protenix's own AF3-style ranking quantity, and is the direct
    analogue of Boltz's `confidence_score` and Chai-1's `aggregate_score`:
    each arm is scored on the decoy ITS OWN confidence ranks first, never on
    a quantity borrowed from another model.
    """
    base = COFOLD_ROOT / "protenix" / "out" / pdb_id / "pmhc_input" / "seed_101" / "predictions"
    if not base.exists():
        return None
    samples = []
    for n in range(25):
        pred_path = base / f"pmhc_input_sample_{n}.cif"
        conf_path = base / f"pmhc_input_summary_confidence_sample_{n}.json"
        if not pred_path.exists():
            continue
        conf = None
        if conf_path.exists():
            try:
                with open(conf_path) as f:
                    conf = json.load(f).get("ranking_score")
            except Exception as e:
                print(f"WARNING {pdb_id}: could not read {conf_path}: {e}", file=sys.stderr)
        samples.append({
            "decoy": f"sample_{n}",
            "path": pred_path,
            "confidence": conf,
        })
    if not samples:
        return None
    scored = [s for s in samples if s["confidence"] is not None]
    if not scored:
        return samples, None
    # ties -> lowest sample index, same convention as enumerate_samples_boltz
    best = max(scored, key=lambda s: (s["confidence"], -int(s["decoy"].split("_")[1])))
    return samples, best["decoy"]


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <arm> <pdb_id>", file=sys.stderr)
        sys.exit(1)
    arm, pdb_id = sys.argv[1], sys.argv[2]
    if arm not in ARMS:
        print(f"ERROR: unknown arm {arm!r}, expected one of {ARMS}", file=sys.stderr)
        sys.exit(1)

    peptide_seq = load_peptide_seq(pdb_id)
    if not peptide_seq:
        print(f"SKIP {pdb_id}: no peptide_seq found in joblist.tsv", file=sys.stderr)
        return

    ref_path = REFS_DIR / f"{pdb_id}.cif"
    if not ref_path.exists():
        print(f"SKIP {pdb_id}: no reference")
        return

    if arm == "af2_reuse":
        result = enumerate_samples_af2_reuse(pdb_id)
    elif arm in ("boltz2", "boltz1"):
        result = enumerate_samples_boltz(arm, pdb_id)
    elif arm == "chai1":
        result = enumerate_samples_chai1(pdb_id)
    elif arm == "protenix":
        result = enumerate_samples_protenix(pdb_id)
    else:
        result = None

    if result is None:
        print(f"SKIP {pdb_id}: no samples found for arm {arm}", file=sys.stderr)
        return

    samples, selected_decoy = result
    if not samples:
        print(f"SKIP {pdb_id}: no samples found for arm {arm}", file=sys.stderr)
        return

    # score each sample
    decoy_rows = []
    scored = {}  # decoy -> rmsd
    for s in samples:
        decoy = s["decoy"]
        try:
            res = score_prediction(s["path"], ref_path, peptide_seq, copy_policy="best")
        except Exception as e:
            print(f"WARNING {pdb_id}/{decoy} ({arm}): scoring failed: {e}", file=sys.stderr)
            continue
        rmsd = res.get("peptide_backbone_rmsd")
        decoy_rows.append({
            "pdb_id": pdb_id,
            "peptide_seq": peptide_seq,
            "decoy": decoy,
            "peptide_backbone_rmsd": rmsd,
            "confidence": s["confidence"],
            "n_atom_pairs": res.get("n_atom_pairs"),
            "n_atoms_dropped": res.get("n_atoms_dropped"),
            "n_paired_ca": res.get("n_paired_ca"),
            "superposition_rmsd": res.get("superposition_rmsd"),
            "warnings": ";".join(res.get("warnings") or []),
        })
        if rmsd is not None:
            scored[decoy] = rmsd

    if not scored:
        print(f"SKIP {pdb_id}: no samples successfully scored for arm {arm}", file=sys.stderr)
        return

    n_scored = len(scored)
    rmsds_sorted = sorted(scored.values())
    oracle = rmsds_sorted[0]
    median = rmsds_sorted[n_scored // 2]

    selected_rmsd = scored.get(selected_decoy) if selected_decoy is not None else None
    if selected_rmsd is None:
        rank_of_selected = ""
        percentile_of_selected = ""
    else:
        # 1-indexed rank in ascending order; ties broken by stable order of rmsds_sorted
        rank_of_selected = sorted(scored.values()).index(selected_rmsd) + 1
        if n_scored <= 1:
            percentile_of_selected = float("nan")
        else:
            percentile_of_selected = (rank_of_selected - 1) / (n_scored - 1)

    summary_row = {
        "pdb_id": pdb_id,
        "peptide_seq": peptide_seq,
        "n_scored": n_scored,
        "selected": selected_rmsd if selected_rmsd is not None else "",
        "oracle": oracle,
        "median": median,
        "rank_of_selected": rank_of_selected,
        "percentile_of_selected": percentile_of_selected,
    }

    out_dir = COFOLD_ROOT / arm / "results" / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / f"{pdb_id}.summary.csv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_HEADER)
        writer.writeheader()
        writer.writerow(summary_row)

    decoys_path = out_dir / f"{pdb_id}.decoys.csv"
    with open(decoys_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=DECOY_HEADER)
        writer.writeheader()
        for row in decoy_rows:
            writer.writerow(row)

    print(f"OK {pdb_id} ({arm}): n_scored={n_scored} oracle={oracle:.4f} selected={selected_rmsd}")


if __name__ == "__main__":
    main()
