#!/usr/bin/env bash
#SBATCH --job-name=calib
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=120
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/calib_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/calibrate_thresholds.py --out-dir results/calibration
echo "### DONE $(date)"
