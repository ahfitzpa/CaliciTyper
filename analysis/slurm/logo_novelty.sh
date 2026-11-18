#!/usr/bin/env bash
#SBATCH --job-name=logo_novel
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=120
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/logo_novel_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/logo_novelty.py --out results/logo_novelty.csv
echo "### DONE $(date)"
