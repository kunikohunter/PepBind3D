#!/bin/bash --norc
# run_protenix.sh
#
# Reusable wrapper to run a single Protenix structure prediction for a
# peptide-MHC class I complex. Mirrors run_boltz2.sh and run_chai1_msa.py so
# that the fifth co-folding arm is built under the same conditions as the
# other four.
#
# Usage:
#   run_protenix.sh <peptide_seq> <mhc_seq> <output_dir> [extra protenix args]
#
# MODEL / LEAKAGE CONTROL  -- read before changing anything here.
#
# The model name is pinned to `protenix_base_default_v1.0.0`, which ByteDance
# trained on a 2021-09-30 wwPDB cutoff. Verified 2026-09-08 from the installed
# package itself:
#   configs/configs_model_type.py:24
#       "models are trained based on the 2021-09-30 wwPDB cutoff"
#   configs/configs_data.py:147
#       indices/weightedPDB_indices_before_2021-09-30_wo_posebusters_...csv.gz
#   protenix/data/core/parser.py:2098
#       "our training set which had a maximum release date of 2021-09-30"
# All 174 benchmark targets were released 2024-01-10..2026-08-26, so the
# margin is >2.2 years. This is the same cutoff as AF2 multimer_v3.
#
# DO NOT substitute `protenix_base_20250630_v1.0.0`. It has the same parameter
# count and the same v1.0.0 suffix, but it is trained to a 2025-06-30 cutoff,
# which would put nearly every benchmark target inside its training set and
# silently invalidate the arm. ByteDance's own docs say to use `default` for
# fair cross-version benchmarking. The name is asserted below rather than left
# to the CLI default, so a future change to that default cannot leak in
# unnoticed.
#
# MSA: ON. Mode is `protenix` (ByteDance's own MSA service) by default.
#
# DO NOT set this to `colabfold` for a two-chain pMHC input. That was tried on
# 2026-09-08 (job 13883076) for parity with Boltz-2 and Chai-1, which both query
# api.colabfold.com, and 160 of 172 targets FAILED. The colabfold branch in
# protenix/web_service/colab_request_parser.py:279 loops over chains with
# use_pairing=False, then issues a third, paired/complex request that
# api.colabfold.com answers with "Server didn't reply with json: Not Found",
# and the run dies with 'msa search failed'. That branch also never writes the
# pairing.a3m / non_pairing.a3m files that runner/msa_search.py:155 update_seq_msa
# then looks for, so colabfold mode is not a working production path here even
# when the network call succeeds.
#
# Consequence for the paper: Protenix's MSAs come from a different MMseqs2
# service than Boltz-2's and Chai-1's. Both search UniRef plus environmental
# databases, but this is a genuine methodological difference between arms and
# must be stated in the methods, not glossed. Protenix warns that performance
# "might degrade significantly" without MSAs, so MSA-off is not a safe
# alternative -- it would be a rigged comparison, not a conservative one.
#
# TEMPLATES: OFF (`--use_template False`). v1.0.0 is the first Protenix release
# that supports templates; every arm in this benchmark runs templates off, and
# templates-off is also a second, independent leakage barrier.
#
# SAMPLES: 25 diffusion samples from a single seed (`-e 25 -s 101`), matching
# Boltz-2's 25 diffusion samples. Chai-1 reaches 25 as 5 trunks x 5 samples.
#
# CYCLES/STEPS: `-c 10 -p 200` are the package defaults for this checkpoint,
# passed explicitly for the same reason the model name is: so the invocation
# records them. (The Sep-7 smoke test used -c 4 -p 20; those were fast-test
# values and must NOT be used for production.)
#
# PROTENIX_ROOT_DIR: Protenix defaults this to $HOME and will re-download the
# 1.5 GB checkpoint and ~650 MB of data caches there on first use. Pointed at
# /data below, where both are already staged, so nothing is downloaded and
# nothing lands against the home quota.
#
# Must run on a GPU node. ~30 s/target at these settings on an A6000.

set -eo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <peptide_seq> <mhc_seq> <output_dir> [extra protenix args...]" >&2
    exit 1
fi

PEPTIDE_SEQ="$1"
MHC_SEQ="$2"
OUT_DIR="$3"
shift 3
EXTRA_ARGS=("$@")

MODEL_NAME="protenix_base_default_v1.0.0"
# Override only for deliberate experiments; see the MSA note above.
MSA_MODE="${PROTENIX_MSA_MODE:-protenix}"
PROTENIX_BIN=/data/p_csb_meiler/huntek1/envs/protenix/bin/protenix

export PYTHONNOUSERSITE=1
export PROTENIX_ROOT_DIR=/data/p_csb_meiler/huntek1/envs/protenix_root

# Fail loudly if the pinned checkpoint is not already staged, rather than
# letting Protenix silently download something over the network.
CKPT="$PROTENIX_ROOT_DIR/checkpoint/${MODEL_NAME}.pt"
if [[ ! -f "$CKPT" ]]; then
    echo "ERROR: pinned checkpoint not found at $CKPT" >&2
    echo "Refusing to run: an auto-download could fetch a different cutoff." >&2
    exit 1
fi

# Guard against a caller sneaking a different model in through EXTRA_ARGS.
for a in "${EXTRA_ARGS[@]}"; do
    if [[ "$a" == *"20250630"* ]]; then
        echo "ERROR: refusing to run the 2025-06-30-cutoff checkpoint." >&2
        echo "All 174 benchmark targets postdate that cutoff -- the arm would be leaked." >&2
        exit 1
    fi
done

mkdir -p "$OUT_DIR"
INPUT_JSON="$OUT_DIR/pmhc_input.json"

# Chain layout matches every other arm: MHC heavy chain as A, peptide as B,
# no beta-2-microglobulin. Verified against the Boltz-2 and Chai-1 inputs.
cat > "$INPUT_JSON" << EOF
[
  {
    "name": "pmhc_input",
    "sequences": [
      {"proteinChain": {"sequence": "${MHC_SEQ}", "count": 1, "id": ["A"]}},
      {"proteinChain": {"sequence": "${PEPTIDE_SEQ}", "count": 1, "id": ["B"]}}
    ]
  }
]
EOF

echo "Wrote input JSON to $INPUT_JSON"
echo "Running Protenix (model=$MODEL_NAME, 25 samples, MSA=$MSA_MODE, templates=off) -> $OUT_DIR"

"$PROTENIX_BIN" pred \
    -i "$INPUT_JSON" \
    -o "$OUT_DIR" \
    --model_name "$MODEL_NAME" \
    -s 101 \
    -c 10 \
    -p 200 \
    -e 25 \
    --use_msa True \
    --msa_server_mode "$MSA_MODE" \
    --use_template False \
    "${EXTRA_ARGS[@]}"

echo "Done. Results under $OUT_DIR/pmhc_input/seed_101/predictions"
