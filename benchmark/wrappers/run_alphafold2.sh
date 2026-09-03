#!/bin/bash --norc
# run_alphafold2.sh
#
# Reusable wrapper to run a single AlphaFold2-Multimer structure prediction
# for a peptide-MHC class I complex, using the cluster's validated AF2 2.3.2
# install. Modeled on the cluster's own working driver script
# /sb/apps/alphafold232/scripts/run_alphafold_ACCRE_Ampere.sb (env setup,
# especially the LD_LIBRARY_PATH export, is copied from there), but this
# wrapper targets the MULTIMER model preset, not that script's monomer
# default, because a peptide-MHC complex is a two-chain hetero-complex.
#
# Install:
#   AF2_MINICONDA=/sb/apps/alphafold232/miniconda3   (conda env `af232`)
#   AF2_DATADIR=/sb/apps/alphafold-data.230
#   AF2_REPO=/sb/apps/alphafold232/alphafold
#
# PRESET: --model_preset=multimer (decided; not re-litigated here).
#   A pMHC complex is a two-chain hetero-complex, which is exactly what the
#   multimer preset is built for (monomer folds a single chain in
#   isolation and cannot model the peptide-MHC interface at all).
#
# WEIGHTS / TRAINING CUTOFF: multimer uses model_{1..5}_multimer_v3 weights
#   (verified present: params_model_{1,2,3,4,5}_multimer_v3.npz under
#   $AF2_DATADIR/params/), trained to 2021-09-30 -- NOT the 2018-04-30 cutoff
#   of the monomer weights. 2021-09-30 is the number for Methods. It is
#   still comfortably before the 2024-01-01 benchmark holdout.
#
# TEMPLATES / LEAKAGE: the multimer path does NOT use PDB70 (that database
#   is monomer-only and pinned at 2020-04-01, which is why the monomer path
#   was leakage-safe "by accident"). Multimer instead searches
#   pdb_seqres_database_path (hmmsearch) + pdb_mmcif for templates, and
#   pdb_seqres/uniprot are NOT date-pinned the way PDB70 is -- so
#   --max_template_date is now doing REAL work here, not a redundant
#   belt-and-suspenders flag. This wrapper hard-codes
#   --max_template_date=2023-12-31 (NOT 9999-12-31, unlike the reference
#   driver script) so that no structure released on/after the 2024-01-01
#   benchmark holdout date can ever be used as a template, regardless of
#   what is currently in RCSB. Do not change this date without re-review.
#   NOTE: passing --pdb70_database_path under model_preset=multimer is a
#   hard error in this install (run_alphafold.py:345-348 asserts
#   pdb70_database_path is UNSET for multimer and pdb_seqres_database_path /
#   uniprot_database_path ARE set) -- so pdb70 is intentionally omitted
#   below and pdb_seqres_database_path / uniprot_database_path are passed
#   instead. Verified present on disk:
#     $AF2_DATADIR/pdb_seqres/pdb_seqres.txt
#     $AF2_DATADIR/uniprot/uniprot.fasta
#
# MSA: full_dbs preset (uniref90 + mgnify + uniref30 + bfd), i.e. MSAs ON,
#   per the benchmark's primary condition. AF2 has no "remote MSA server"
#   option like Boltz/Chai -- its MSAs always come from local HMMER/HHblits
#   searches against the databases below, which is its native/only MSA
#   mechanism, so nothing further is required to turn MSAs "on" here.
#
# SAMPLES: --num_multimer_predictions_per_model=5. This install ships 5
#   multimer_v3 model weights (model_1..model_5), and this flag multiplies:
#   5 models x 5 predictions-per-model (each a different random seed) = 25
#   total predictions per target, matching the benchmark's 25-samples
#   requirement using AF2's own native multi-sample mechanism. (Confirmed
#   in run_alphafold.py:119-123 this flag ONLY takes effect when
#   model_preset=multimer -- under monomer it is silently inert.)
#
# RELAXATION / SCORED SET: measured on the first test run that AF2 Amber-
#   relaxes only the top-ranked prediction (ranked_0.pdb was byte-identical
#   to relaxed_model_3_multimer_v3_pred_2.pdb; ranked_1..24 matched the
#   unrelaxed files). Scoring ranked_* would mean 1 of every 25 samples was
#   geometry-minimised and 24 were not -- a systematic inconsistency that
#   specifically advantages the best-ranked sample and does not match any
#   other model in the panel (none of which relax). So this wrapper passes
#   --models_to_relax=none (flag defined in run_alphafold.py, default
#   `best`; confirmed via `flags.DEFINE_enum_class('models_to_relax', ...)`
#   in run_alphafold.py) to skip Amber relaxation entirely -- this also
#   saves the relax wall-clock time and means no relaxed_*.pdb or ranked_*
#   file is ever produced, so there is nothing to accidentally pick up.
#   The canonical scored set is the 25
#   unrelaxed_model_{1..5}_multimer_v3_pred_{0..4}.pdb files; the wrapper
#   calls af2_postprocess.py after a successful run to write a manifest
#   (MANIFEST_unrelaxed.txt) and a per-prediction summary (summary.json:
#   model_name, ranking_confidence, path_to_unrelaxed_pdb) naming exactly
#   those 25 files, then deletes result_*.pkl (see that script's docstring
#   for the safety check it runs first).
#
# Usage:
#   run_alphafold2.sh <peptide_seq> <mhc_seq> <output_dir>
#
# IMPORTANT: this must be run on a GPU node (via sbatch or salloc), never
# on the login node. This script activates its own env, so it can be
# sourced/run directly inside an salloc GPU shell, or wrapped in a small
# sbatch script that just calls it.

set -eo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <peptide_seq> <mhc_seq> <output_dir>" >&2
    exit 1
fi

PEPTIDE_SEQ="$1"
MHC_SEQ="$2"
OUT_DIR="$3"

AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/sb/apps/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold

source "$AF2_MINICONDA/bin/activate" af232
export LD_LIBRARY_PATH="$AF2_MINICONDA/envs/af232/lib:$LD_LIBRARY_PATH"

mkdir -p "$OUT_DIR"
FASTA="$OUT_DIR/pmhc_input.fasta"

# Multimer input: BOTH chains as separate FASTA records in one file (AF2
# multimer detects >1 sequence in the file and folds them as a complex).
# Chain order does not matter to AF2, but is kept MHC-then-peptide to match
# the boltz1/boltz2/chai wrappers' A=MHC, B=peptide convention.
cat > "$FASTA" << EOF
>mhc
${MHC_SEQ}
>peptide
${PEPTIDE_SEQ}
EOF

echo "Wrote input FASTA to $FASTA"
echo "Running AlphaFold2-Multimer (model_preset=multimer, 5 models x 5" \
     "predictions/model = 25 samples, max_template_date=2023-12-31) -> $OUT_DIR"

python "$AF2_REPO/run_alphafold.py" \
    --fasta_paths="$FASTA" \
    --data_dir="$AF2_DATADIR" \
    --output_dir="$OUT_DIR" \
    --model_preset=multimer \
    --num_multimer_predictions_per_model=5 \
    --max_template_date=2023-12-31 \
    --models_to_relax=none \
    --uniref90_database_path="$AF2_DATADIR/uniref90/uniref90.fasta" \
    --mgnify_database_path="$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa" \
    --uniref30_database_path="$AF2_DATADIR/uniref30/UniRef30_2021_03" \
    --bfd_database_path="$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt" \
    --pdb_seqres_database_path="$AF2_DATADIR/pdb_seqres/pdb_seqres.txt" \
    --uniprot_database_path="$AF2_DATADIR/uniprot/uniprot.fasta" \
    --template_mmcif_dir="$AF2_DATADIR/pdb_mmcif/mmcif_files" \
    --obsolete_pdbs_path="$AF2_DATADIR/pdb_mmcif/obsolete.dat"

TARGET_DIR="$OUT_DIR/$(basename "$FASTA" .fasta)"
echo "Done. Results under $TARGET_DIR/"

# Post-run manifest + storage reduction (see af2_postprocess.py docstring).
# Confirms ranking_debug.json has all 25 prediction keys and their unrelaxed
# pdb files exist, writes MANIFEST_unrelaxed.txt + summary.json naming the
# canonical scored set, and only then deletes the large result_*.pkl files.
# If the check fails, nothing is deleted and this script exits non-zero.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$SCRIPT_DIR/af2_postprocess.py" "$TARGET_DIR"
