#!/bin/bash
#SBATCH --job-name=hlac-stage
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --output=/data/p_csb_meiler/huntek1/hlac_production/logs/stage.%j.out
#SBATCH --error=/data/p_csb_meiler/huntek1/hlac_production/logs/stage.%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

set -e
/data/p_csb_meiler/huntek1/envs/af2/bin/python /data/p_csb_meiler/huntek1/hlac_production/jobs/stage_and_verify.py
