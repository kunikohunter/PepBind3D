#!/usr/bin/env python3
"""af2_postprocess.py -- post-run cleanup/manifest step for run_alphafold2.sh.

Called once, after run_alphafold.py exits 0, on the AF2 per-target output
directory (the directory containing ranking_debug.json, e.g.
`$OUT_DIR/pmhc_input/`).

Does exactly three things, in this order, and stops (non-zero exit, no
deletion) if any check fails:

  1. Reads ranking_debug.json and confirms it has all 25 prediction keys
     under "iptm+ptm" (5 models x 5 predictions/model =
     model_{1..5}_multimer_v3_pred_{0..4}), and confirms the corresponding
     unrelaxed_model_*.pdb file exists on disk for each of the 25 keys.
  2. Writes two small files that name the canonical scored set explicitly,
     so downstream code cannot accidentally pick up ranked_* or a relaxed
     structure:
       - MANIFEST_unrelaxed.txt : one unrelaxed pdb filename per line,
         exactly the 25 files that constitute the scored ensemble.
       - summary.json : one row per prediction with model_name,
         ranking_confidence (the iptm+ptm value), and path_to_unrelaxed_pdb.
  3. Only after (1) and (2) succeed: deletes result_*.pkl (the ~27 MB/each
     per-residue PAE/plDDT arrays this benchmark does not use). If step 1
     fails, no pkls are touched -- losing them without having captured the
     confidences would be unrecoverable without a re-run.

Usage:
    python3 af2_postprocess.py <af2_target_dir>

Exit codes: 0 on success (manifest+summary written, pkls removed).
            1 on any check failure (nothing deleted).
"""

import json
import sys
from pathlib import Path

MODELS = range(1, 6)          # model_1 .. model_5
PREDS_PER_MODEL = range(0, 5)  # pred_0 .. pred_4


def expected_keys():
    return [f"model_{m}_multimer_v3_pred_{p}" for m in MODELS for p in PREDS_PER_MODEL]


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <af2_target_dir>", file=sys.stderr)
        return 1

    target_dir = Path(sys.argv[1])
    ranking_path = target_dir / "ranking_debug.json"

    if not ranking_path.exists():
        print(f"CHECK FAILED: {ranking_path} does not exist. Not deleting anything.",
              file=sys.stderr)
        return 1

    with open(ranking_path) as f:
        ranking = json.load(f)

    conf = ranking.get("iptm+ptm")
    if not isinstance(conf, dict):
        print(f"CHECK FAILED: {ranking_path} has no 'iptm+ptm' dict. "
              f"Not deleting anything.", file=sys.stderr)
        return 1

    exp_keys = expected_keys()
    missing_keys = [k for k in exp_keys if k not in conf]
    if missing_keys or len(conf) != 25:
        print(f"CHECK FAILED: expected 25 'iptm+ptm' keys "
              f"(model_{{1..5}}_multimer_v3_pred_{{0..4}}), found {len(conf)}. "
              f"Missing: {missing_keys}. Not deleting anything.", file=sys.stderr)
        return 1

    rows = []
    missing_pdbs = []
    for key in exp_keys:
        pdb_path = target_dir / f"unrelaxed_{key}.pdb"
        if not pdb_path.exists():
            missing_pdbs.append(str(pdb_path))
            continue
        rows.append({
            "model_name": key,
            "ranking_confidence": conf[key],
            "path_to_unrelaxed_pdb": str(pdb_path),
        })

    if missing_pdbs:
        print(f"CHECK FAILED: {len(missing_pdbs)} unrelaxed pdb file(s) missing "
              f"on disk despite ranking_debug.json entries: {missing_pdbs}. "
              f"Not deleting anything.", file=sys.stderr)
        return 1

    # All checks passed -- write manifest + summary, then delete result_*.pkl.
    manifest_path = target_dir / "MANIFEST_unrelaxed.txt"
    with open(manifest_path, "w") as f:
        f.write(
            "# Canonical scored set for this target: exactly these 25\n"
            "# unrelaxed_model_{1..5}_multimer_v3_pred_{0..4}.pdb files.\n"
            "# Do NOT use ranked_*.pdb (rank-reordered, and ranked_0 is the\n"
            "# lone Amber-relaxed structure -- see run_alphafold2.sh header).\n"
        )
        for row in rows:
            f.write(Path(row["path_to_unrelaxed_pdb"]).name + "\n")

    summary_path = target_dir / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(rows, f, indent=2)

    deleted = []
    for pkl in sorted(target_dir.glob("result_*.pkl")):
        pkl.unlink()
        deleted.append(pkl.name)

    print(f"OK: wrote {manifest_path.name} and {summary_path.name} "
          f"(25 predictions). Deleted {len(deleted)} result_*.pkl files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
