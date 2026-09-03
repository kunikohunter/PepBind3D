#!/bin/bash
#SBATCH --job-name=hlac-dock
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --output=/data/p_csb_meiler/huntek1/hlac_production/logs/dock_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/huntek1/hlac_production/logs/dock_%A_%a.err
#SBATCH --time=02:30:00
#SBATCH --mem=1G   # matches production flexpep_array_job.sh
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=1-1861%100

set -e
SHIM="/data/p_csb_meiler/huntek1/dock_validation/shim"
TASKLIST="/data/p_csb_meiler/huntek1/hlac_production/jobs/dock_tasklist.txt"
STAGE_ROOT="/data/p_csb_meiler/huntek1/hlac_production/pdb"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKLIST")
COMPACT=$(echo "$LINE" | cut -f1)
PEPTIDE=$(echo "$LINE" | cut -f2)

if [ -z "$COMPACT" ] || [ -z "$PEPTIDE" ]; then
  echo "FATAL: could not resolve task ${SLURM_ARRAY_TASK_ID} from $TASKLIST" >&2
  exit 1
fi

OUTDIR="$STAGE_ROOT/$COMPACT/$PEPTIDE"
INPUT_PDB="$OUTDIR/${PEPTIDE}_input.pdb"

echo "=== Docking task ${SLURM_ARRAY_TASK_ID}: allele=$COMPACT peptide=$PEPTIDE ==="

if [ ! -f "$INPUT_PDB" ]; then
  echo "FATAL: missing threaded input $INPUT_PDB (threading stage for $COMPACT did not stage this peptide)" >&2
  exit 1
fi

cd "$OUTDIR"

"$SHIM/source/bin/FlexPepDocking.linuxgccrelease" \
  -in:file:s "${PEPTIDE}_input.pdb" \
  -database "$SHIM/database" \
  -pep_refine \
  -nstruct 25 \
  -ex1 -ex2aro \
  -out:path:all "$OUTDIR"

# Verify the expected 25 decoys + score.sc landed, matching the A/B convention
NDECOY=$(ls "$OUTDIR"/${PEPTIDE}_input_00*.pdb 2>/dev/null | wc -l)
if [ "$NDECOY" -ne 25 ] || [ ! -f "$OUTDIR/score.sc" ]; then
  echo "FATAL: expected 25 decoys + score.sc in $OUTDIR, found $NDECOY decoys, score.sc present=$([ -f "$OUTDIR/score.sc" ] && echo yes || echo no)" >&2
  exit 1
fi
echo "=== OK: $PEPTIDE ($COMPACT) docked, 25 decoys + score.sc staged ==="
