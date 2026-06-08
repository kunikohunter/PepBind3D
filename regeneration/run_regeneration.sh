#!/bin/bash
# Regenerate FlexPepDock ensembles for the 52 validation pairs using
# --ignore_epitope_match so that self-matching templates are excluded.
#
# Run this from /home/huntek1/main_project/scripts/IEDB_validation/regeneration/
# (or any working directory; paths below are absolute).
#
# Prerequisites:
#   1. The patched HLA_db.py is in /home/huntek1/main_project/scripts/
#   2. build_validation_fastas.py has been run to produce per-allele FASTAs.
#   3. The Python environment used has biopython, pandas, etc. installed
#      (the same one used for the validation notebooks).

set -euo pipefail

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
PIPELINE_SCRIPT=/home/huntek1/main_project/scripts/IEDBTestPipeline_ACCRE.py
ROSETTA_PATH=/sb/meilerapps/rosetta/rosetta-3.15/main
FASTA_DIR=/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration/fastas
MANIFEST=/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration/allele_manifest.csv
LOG_DIR=/home/huntek1/main_project/data/IEDB_data_clean/IEDB_validation/regeneration/logs

THREADS=28

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
mkdir -p "$LOG_DIR"

if [ ! -f "$MANIFEST" ]; then
    echo "Manifest not found: $MANIFEST"
    echo "Run build_validation_fastas.py first."
    exit 1
fi

echo "Regeneration starting at $(date)"
echo "Threads:        $THREADS"
echo "Rosetta:        $ROSETTA_PATH"
echo "Pipeline:       $PIPELINE_SCRIPT"
echo

# -----------------------------------------------------------------------------
# Iterate over alleles in the manifest
# -----------------------------------------------------------------------------
# Skip the header row, then for each row pull the allele and FASTA path.
tail -n +2 "$MANIFEST" | while IFS=, read -r ALLELE ALLELE_DIR N_PEPTIDES FASTA_PATH; do
    ALLELE=$(echo "$ALLELE" | tr -d '"')
    FASTA_PATH=$(echo "$FASTA_PATH" | tr -d '"')

    LOG_FILE="$LOG_DIR/${ALLELE_DIR}.log"

    echo "============================================================"
    echo "Allele:    $ALLELE  ($N_PEPTIDES peptides)"
    echo "FASTA:     $FASTA_PATH"
    echo "Log:       $LOG_FILE"
    echo "Starting:  $(date)"

    python "$PIPELINE_SCRIPT" \
        --IEDBquery skip \
        --buildFasta "$FASTA_PATH" 0 \
        --setAllele "$ALLELE" \
        --rosetta "$ROSETTA_PATH" \
        --threads "$THREADS" \
        --ignore_epitope_match \
        2>&1 | tee "$LOG_FILE"

    echo "Completed: $(date)"
done

echo
echo "All alleles done at $(date)"
echo "Threaded inputs in: $FASTA_DIR/output"