#!/usr/bin/env bash
#SBATCH --job-name=gi5
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=60
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/gi5_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/analyse_gi5_ptype.py --out-dir results/gi5_ptype
echo "### DONE $(date)"
