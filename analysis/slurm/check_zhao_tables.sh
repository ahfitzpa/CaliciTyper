#!/usr/bin/env bash
#SBATCH --job-name=zhao_check
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --time=30
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/zhao_check_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/check_against_zhao_tables.py
echo "### exit $? $(date)"
