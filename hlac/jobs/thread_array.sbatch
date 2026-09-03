#!/bin/bash
#SBATCH --job-name=hlac-thread
#SBATCH --output=/data/p_csb_meiler/huntek1/hlac_production/logs/thread_%A_%a.log
#SBATCH --array=1-81%20
#SBATCH --time=04:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
# Concurrency is SAFE here: set_mhc_templates() returns early without writing
# allele_list.info when Templates are already populated (HLA_db.py, ~line 37-39).
# All 11 HLA-C alleles are now populated, so threading is read-only on that file.
PY=/data/p_csb_meiler/huntek1/envs/af2/bin/python
cd /home/huntek1/main_project/scripts || exit 1
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /data/p_csb_meiler/huntek1/hlac_production/thread_tasklist.txt)
echo "[task ${SLURM_ARRAY_TASK_ID}] $CMD"
eval "${CMD/#python/$PY}"
