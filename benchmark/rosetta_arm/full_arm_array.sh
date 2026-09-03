#!/bin/bash
#SBATCH --job-name=rosetta_arm
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=02:30:00
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-170%100
#SBATCH --output=/data/p_csb_meiler/huntek1/benchmark/rosetta_arm/logs/arm-%A_%a.log
#
# Full 174-target Rosetta benchmark arm (170 queued; 4 skipped, see
# arm_skipped.csv). One array task = one target: thread the reference
# peptide onto its resolved-allele receptor (same pipeline/protocol as the
# validated 3-target pilot), then FlexPepDock nstruct 25 on the threaded
# input, using the same Rosetta build as the pilot
# (/data/p_csb_meiler/huntek1/dock_validation/shim ->
#  rosetta.binary.linux.release-371).
#
# IEDBTestPipeline_ACCRE.py's thread_all() uses os.mkdir (not makedirs),
# and safe_thread_all() swallows any exception from it and returns None,
# so a rerun over an existing target dir can silently produce zero new
# decoys while SLURM still reports COMPLETED. This script does NOT trust
# that exit code: it counts actual decoy PDBs at the end and fails loudly
# (nonzero exit) if the count is not exactly 25.
#
# Usage:
#   sbatch /data/p_csb_meiler/huntek1/benchmark/rosetta_arm/scripts/full_arm_array.sh <job_list_file>
#
# <job_list_file>: one '|'-delimited record per line (1-indexed), of the form
#   pdb_id|resolved_allele|peptide_seq|target_dir
# target_dir must already contain fasta/<pdb_id>.fasta and an (empty)
# docking/ subdir, as produced by build_arm.py.

set -uo pipefail

SHIM="/data/p_csb_meiler/huntek1/dock_validation/shim"
FLEXPEPDOCK="${SHIM}/source/bin/FlexPepDocking.linuxgccrelease"
ROSETTA_DB="${SHIM}/database"
PY="/data/p_csb_meiler/huntek1/envs/af2/bin/python"
PIPELINE="/home/huntek1/main_project/scripts/IEDBTestPipeline_ACCRE.py"
N_EXPECTED_DECOYS=25

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-<none>}"
echo "Job list file: $1"

JOB_LIST="$1"
if [[ -z "$JOB_LIST" || ! -f "$JOB_LIST" ]]; then
    echo "ERROR: job list file '$JOB_LIST' not found." >&2
    exit 1
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOB_LIST")
if [[ -z "$LINE" ]]; then
    echo "ERROR: no record found at line ${SLURM_ARRAY_TASK_ID} of ${JOB_LIST}" >&2
    exit 1
fi

IFS='|' read -r pdb_id resolved_allele peptide_seq target_dir <<< "$LINE"

echo "Target: ${pdb_id}  allele: ${resolved_allele}  peptide: ${peptide_seq}"
echo "Target dir: ${target_dir}"

FASTA="${target_dir}/fasta/${pdb_id}.fasta"
DOCK_DIR="${target_dir}/docking"

if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: fasta ${FASTA} not found." >&2
    exit 1
fi

# --- Threading (same invocation as the pilot) ---
THREAD_CMD=("$PY" "$PIPELINE" --IEDBquery skip --setAllele "$resolved_allele" \
    --buildFasta "$FASTA" 0 --rosetta "$SHIM")

echo "Running threading command:"
echo "${THREAD_CMD[@]}"
"${THREAD_CMD[@]}"
THREAD_STATUS=$?
echo "Threading command exited with status ${THREAD_STATUS}"

if [[ $THREAD_STATUS -ne 0 ]]; then
    echo "ERROR: threading failed for ${pdb_id} (exit ${THREAD_STATUS})." >&2
    exit 1
fi

INPUT_PDB="${target_dir}/fasta/output/${pdb_id}_batch1/${peptide_seq}/${peptide_seq}_input.pdb"
if [[ ! -f "$INPUT_PDB" ]]; then
    echo "ERROR: expected threaded input ${INPUT_PDB} not found after threading." >&2
    exit 1
fi

mkdir -p "$DOCK_DIR"

# --- FlexPepDock refinement (same invocation as the pilot) ---
DOCK_CMD=("$FLEXPEPDOCK" -in:file:s "$INPUT_PDB" -database "$ROSETTA_DB" \
    -pep_refine -nstruct "$N_EXPECTED_DECOYS" -ex1 -ex2aro \
    -out:path:all "$DOCK_DIR")

echo "Running docking command:"
echo "${DOCK_CMD[@]}"
"${DOCK_CMD[@]}"
DOCK_STATUS=$?
echo "FlexPepDocking exited with status ${DOCK_STATUS}"

# --- Hard-fail check: do not trust the exit code alone ---
N_DECOYS=$(find "$DOCK_DIR" -maxdepth 1 -name "${peptide_seq}_input_[0-9][0-9][0-9][0-9].pdb" | wc -l)
echo "Decoys found for ${pdb_id}: ${N_DECOYS} (expected ${N_EXPECTED_DECOYS})"

if [[ "$N_DECOYS" -ne "$N_EXPECTED_DECOYS" ]]; then
    echo "ERROR: ${pdb_id} ended with ${N_DECOYS} decoys, expected ${N_EXPECTED_DECOYS}. Failing loudly." >&2
    exit 1
fi

echo "${pdb_id} completed successfully with ${N_DECOYS} decoys."
exit 0
