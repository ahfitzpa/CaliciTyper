#!/usr/bin/env bash
#SBATCH --job-name=recomb_t3
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --array=0-12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/recomb_t3_%A_%a.out
#
# Recombination validation against Table 3 of Zhao et al. 2024: thirteen real
# sapovirus strains whose VP1 and RdRp types disagree. Twelve were confirmed
# recombinant at the RdRp/VP1 junction by RDP4 and Simplot; the sea lion strain
# JN420370 was NOT, and serves as the negative control. Until now the
# recombination path had only been tested on synthetic alignments.
#
# One array task per strain, because OpenRDP takes tens of minutes per run.
#     sbatch recomb_table3.sh
# then, once all tasks finish:
#     python scripts/validate_recombination_table3.py --summarise --out-dir results/recomb_table3
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export CALICITYPER_OPENRDP=/home/people/fitzpatria/scratch/openrdp_env/bin/openrdp
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export PYTHONUNBUFFERED=1
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | task ${SLURM_ARRAY_TASK_ID} | commit $(git rev-parse --short HEAD)"
# stagger the GenBank fetches so thirteen tasks do not hit NCBI at once
sleep $(( SLURM_ARRAY_TASK_ID * 4 ))
python scripts/validate_recombination_table3.py \
    --index "${SLURM_ARRAY_TASK_ID}" \
    --out-dir results/recomb_table3 \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
echo "### exit $? $(date)"
