#!/usr/bin/env bash
#SBATCH --job-name=temporal
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=120
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/temporal_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/temporal_holdout_novelty.py --cutoff 2016-01-01 --out-dir results/temporal_holdout
echo "### DONE $(date)"
