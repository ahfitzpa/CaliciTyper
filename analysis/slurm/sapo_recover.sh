#!/usr/bin/env bash
#SBATCH --job-name=sapo_recover
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/sapo_recover_%j.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export PYTHONPATH=/home/people/fitzpatria/scratch/calicityper
export PYTHONUNBUFFERED=1
cd /home/people/fitzpatria/scratch/calicityper
B=/home/people/fitzpatria/scratch/calicityper_runs/sapovirus_build
echo "### host $(hostname) $(date)"
python scripts/recover_sapovirus_rdrp.py \
    --genbank "$B/sapovirus.gb" --out-dir "$B/recovery"
echo "### DONE $(date)"
