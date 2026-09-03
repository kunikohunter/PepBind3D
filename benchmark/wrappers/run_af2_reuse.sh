#!/bin/bash --norc
# run_af2_reuse.sh
#
# AF2-multimer pMHC wrapper that reuses a cached MHC-chain MSA across
# benchmark targets sharing the same HLA allele, instead of recomputing it
# per target. Built on top of, and does NOT modify,
# /data/p_csb_meiler/huntek1/benchmark/wrappers/run_alphafold2.sh.
#
# WHY THIS EXISTS:
# alphafold/data/pipeline_multimer.py assigns each FASTA record its own
# chain letter (A, B, ...) in input order and runs the ENTIRE monomer
# search pipeline (jackhmmer vs uniref90/mgnify, hhblits vs bfd+uniref30,
# hmmsearch vs pdb_seqres, jackhmmer vs uniprot for pairing) independently
# per chain, writing results under <output_dir>/msas/<chain_id>/. Nothing
# in that pipeline is pair- or complex-specific -- see pipeline_multimer.py
# DataPipeline._process_single_chain() and _all_seq_msa_features(): both
# operate on a single chain's own sequence and write into that chain's own
# msa_output_dir. So a chain's MSA depends only on that chain's sequence,
# not on what it is complexed with.
#
# run_alphafold.py exposes:
#   flags.DEFINE_boolean('use_precomputed_msas', False, 'Whether to read MSAs
#   that have been written to disk instead of running the MSA tools. The MSA
#   files are looked up in the output directory, so it must stay the same
#   between multiple runs that are to reuse the MSAs. WARNING: This will not
#   check if the sequence, database or configuration have changed.')
# and alphafold/data/pipeline.py run_msa_tool():
#   if not use_precomputed_msas or not os.path.exists(msa_out_path):
#     <run the search tool>
#   else:
#     <read msa_out_path from disk>
# i.e. if the expected file already exists at <output_dir>/msas/<chain>/...,
# the search is skipped entirely and the file is read as-is. Because it does
# NOT check the sequence, WE are responsible for only staging a cached chain
# A directory when the MHC sequence actually matches -- done below via
# chain_id_map.json comparison.
#
# In our wrapper convention (matching run_alphafold2.sh) the FASTA always
# lists MHC first, so MHC is always chain A and peptide is always chain B.
# This script:
#   1. Looks up (or creates) a cache slot for the allele under
#      $MSA_CACHE_DIR/<allele>/msas/A/.
#   2. If a cache hit exists and its cached MHC sequence matches mhc_seq,
#      stages those files into <output_dir>/msas/A/ before running AF2 with
#      --use_precomputed_msas=true. Chain B (peptide) has no staged files,
#      so AF2 computes it fresh, as it must (peptide differs per target).
#   3. If no cache hit, runs normally (both chains computed), then copies
#      the resulting msas/A/ into the cache for future targets sharing this
#      allele.
#
# CAVEATS (see the correctness proof in af2_msa/test/):
#  - use_precomputed_msas does NOT itself verify the sequence matches the
#    cached files -- this wrapper does that check (via chain_id_map.json)
#    before staging, and refuses to stage on mismatch. Do not bypass this.
#  - Only the MHC chain (heavy-chain-only, single sequence, single allele)
#    is reused. The 8-15aa peptide chain, and the paired (uniprot) search,
#    are recomputed every run. Peptide search is NOT free even though the
#    peptide is short: search wall time is dominated by database size, not
#    query length (measured ~21.5 min for a 9-mer peptide chain B vs ~35 min
#    for the ~180aa MHC chain A in the reference test run -- see REPORT
#    notes). So expect a partial, not proportional, speedup.
#
# Usage:
#   run_af2_reuse.sh <peptide_seq> <mhc_seq> <allele> <output_dir> <msa_cache_dir>

set -eo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <peptide_seq> <mhc_seq> <allele> <output_dir> <msa_cache_dir>" >&2
    exit 1
fi

PEPTIDE_SEQ="$1"
MHC_SEQ="$2"
ALLELE="$3"
OUT_DIR="$4"
MSA_CACHE_DIR="$5"

AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/sb/apps/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold

source "$AF2_MINICONDA/bin/activate" af232
export LD_LIBRARY_PATH="$AF2_MINICONDA/envs/af232/lib:$LD_LIBRARY_PATH"

mkdir -p "$OUT_DIR"
FASTA="$OUT_DIR/pmhc_input.fasta"

# MHC first (-> chain A), peptide second (-> chain B): matches
# run_alphafold2.sh and pipeline_multimer.py's in-order A/B/... assignment.
cat > "$FASTA" << EOF
>mhc
${MHC_SEQ}
>peptide
${PEPTIDE_SEQ}
EOF

RUN_NAME="$(basename "$FASTA" .fasta)"
RUN_OUT_DIR="$OUT_DIR/$RUN_NAME"
MSAS_DIR="$RUN_OUT_DIR/msas"
CACHE_ALLELE_DIR="$MSA_CACHE_DIR/$ALLELE"
CACHE_A_DIR="$CACHE_ALLELE_DIR/msas/A"

mkdir -p "$MSAS_DIR"

CACHE_HIT=0
if [[ -f "$CACHE_A_DIR/uniref90_hits.sto" && -f "$CACHE_ALLELE_DIR/chain_id_map.json" ]]; then
    CACHED_MHC_SEQ="$(python -c "
import json, sys
with open('$CACHE_ALLELE_DIR/chain_id_map.json') as f:
    d = json.load(f)
print(d['A']['sequence'])
")"
    if [[ "$CACHED_MHC_SEQ" == "$MHC_SEQ" ]]; then
        CACHE_HIT=1
    else
        echo "WARNING: cached MHC sequence for allele $ALLELE does not match" \
             "mhc_seq argument -- ignoring cache, will recompute." >&2
    fi
fi

if [[ "$CACHE_HIT" -eq 1 ]]; then
    echo "Cache HIT for allele $ALLELE: staging precomputed chain-A (MHC)" \
         "MSAs from $CACHE_A_DIR into $MSAS_DIR/A"
    mkdir -p "$MSAS_DIR/A"
    cp -p "$CACHE_A_DIR"/*.sto "$CACHE_A_DIR"/*.a3m "$MSAS_DIR/A/" 2>/dev/null || true
    USE_PRECOMPUTED="true"
else
    echo "Cache MISS for allele $ALLELE: chain A (MHC) MSA will be computed" \
         "fresh and cached afterward at $CACHE_A_DIR"
    USE_PRECOMPUTED="false"
fi

echo "Running AlphaFold2-Multimer -> $OUT_DIR (use_precomputed_msas=$USE_PRECOMPUTED)"

python "$AF2_REPO/run_alphafold.py" \
    --fasta_paths="$FASTA" \
    --data_dir="$AF2_DATADIR" \
    --output_dir="$OUT_DIR" \
    --model_preset=multimer \
    --num_multimer_predictions_per_model=5 \
    --max_template_date=2023-12-31 \
    --use_gpu_relax \
    --use_precomputed_msas="$USE_PRECOMPUTED" \
    --uniref90_database_path="$AF2_DATADIR/uniref90/uniref90.fasta" \
    --mgnify_database_path="$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa" \
    --uniref30_database_path="$AF2_DATADIR/uniref30/UniRef30_2021_03" \
    --bfd_database_path="$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt" \
    --pdb_seqres_database_path="$AF2_DATADIR/pdb_seqres/pdb_seqres.txt" \
    --uniprot_database_path="$AF2_DATADIR/uniprot/uniprot.fasta" \
    --template_mmcif_dir="$AF2_DATADIR/pdb_mmcif/mmcif_files" \
    --obsolete_pdbs_path="$AF2_DATADIR/pdb_mmcif/obsolete.dat"

echo "Done. Results under $RUN_OUT_DIR/"

# Populate/refresh the per-allele cache from this run's chain-A output so
# later targets sharing this allele get a cache hit. Safe to re-run: only
# copies chain A (MHC), never chain B (peptide).
if [[ -d "$MSAS_DIR/A" ]]; then
    mkdir -p "$CACHE_A_DIR"
    cp -p "$MSAS_DIR/A"/*.sto "$MSAS_DIR/A"/*.a3m "$CACHE_A_DIR/" 2>/dev/null || true
    cp -p "$MSAS_DIR/chain_id_map.json" "$CACHE_ALLELE_DIR/chain_id_map.json"
    echo "Cache refreshed for allele $ALLELE at $CACHE_ALLELE_DIR"
fi
