#!/bin/bash --norc
# run_boltz2.sh
#
# Reusable wrapper to run a single Boltz-2 structure prediction for a
# peptide-MHC class I complex. Mirrors
# /data/p_csb_meiler/huntek1/benchmark/boltz1/run_boltz1.sh, but selects the
# Boltz-2 model explicitly instead of Boltz-1.
#
# Boltz-2 is provided by the SAME `boltz` PyPI package used for Boltz-1
# (package version 2.2.1 at /data/p_csb_meiler/apps/miniforge3/envs/boltz).
# `boltz predict --help` documents `--model [boltz1|boltz2]` with
# "Default is boltz2" -- i.e. Boltz-2 is already the tool's default model,
# but this wrapper passes `--model boltz2` EXPLICITLY so the choice is
# visible in the invocation and cannot silently change if the package's
# own default ever changes in a future version.
#
# Usage:
#   run_boltz2.sh <peptide_seq> <mhc_seq> <output_dir> [extra boltz predict args]
#
# Example:
#   run_boltz2.sh GILGFVFTL \
#     SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRT \
#     /data/p_csb_meiler/huntek1/benchmark/wrappers/test_boltz2/out
#
# MSA: `--use_msa_server` (MMSeqs2 remote server) is ON, per the benchmark's
# primary condition (MSAs on for every model).
#
# Templates: Boltz predict has no template-input flag in this invocation --
# no `--template*` argument is passed anywhere in this wrapper, so no
# structural template is ever supplied. Templates stay OFF by omission.
#
# Samples: `--diffusion_samples 25` is Boltz's own multi-sample mechanism
# (verified in `boltz predict --help`: "The number of diffusion samples to
# use for prediction. Default is 1."), giving exactly 25 samples per
# target in one invocation.
#
# SBATCH GOTCHA (found during testing): boltz's Trainer is PyTorch
# Lightning, which auto-detects a SLURM cluster environment and expects
# distributed-training-shaped resource flags. An sbatch script that sets
# `#SBATCH --ntasks=N` (rather than `--ntasks-per-node=N`) makes Lightning
# raise "You set --ntasks=N ... this variable is not supported" and the
# job dies before predicting anything. Use `--ntasks-per-node=N` (or
# `--ntasks=1`) in any sbatch wrapper around this script.
#
# IMPORTANT: this must be run on a GPU node (via sbatch or salloc), never
# on the login node. Activation command (must be sourced before boltz is
# on PATH):
#   source /data/p_csb_meiler/apps/miniforge3/bin/activate boltz
#
# This script activates the env itself, so it can also be sourced/run
# directly inside an salloc GPU shell.

set -eo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <peptide_seq> <mhc_seq> <output_dir> [extra boltz predict args...]" >&2
    exit 1
fi

PEPTIDE_SEQ="$1"
MHC_SEQ="$2"
OUT_DIR="$3"
shift 3
EXTRA_ARGS=("$@")

BOLTZCONDA=/data/p_csb_meiler/apps/miniforge3
# The conda activation scripts reference some variables (e.g. ADDR2LINE)
# without defaults, which trips `set -u`. Activate with -u relaxed, then
# re-enable it for the rest of the script.
source "$BOLTZCONDA/bin/activate" boltz
set -u

mkdir -p "$OUT_DIR"
FASTA="$OUT_DIR/pmhc_input.fasta"

cat > "$FASTA" << EOF
>A|protein
${MHC_SEQ}
>B|protein
${PEPTIDE_SEQ}
EOF

echo "Wrote input FASTA to $FASTA"
echo "Running Boltz-2 prediction (model=boltz2, 25 diffusion samples) -> $OUT_DIR/predictions"

boltz predict "$FASTA" \
    --out_dir "$OUT_DIR" \
    --model boltz2 \
    --use_msa_server \
    --diffusion_samples 25 \
    "${EXTRA_ARGS[@]}"

echo "Done. Results under $OUT_DIR/predictions"
