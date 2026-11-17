#!/usr/bin/env bash
#SBATCH --job-name=clust_assign
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=60
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/clust_assign_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/test_cluster_assignment.py
echo "### DONE $(date)"
