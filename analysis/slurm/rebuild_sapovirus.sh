#!/usr/bin/env bash
#SBATCH --job-name=sapo_rebuild
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=180
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/sapo_rebuild_%j.out
#
# Regenerate the sapovirus taxonomy tables, rebuild the reference sets and the
# genogroup BLAST database, then re-run every sapovirus benchmark, so the
# effect of a change to the tables is measured rather than assumed.
#
# Run this after ANY change to scripts/build_sapovirus_taxonomy.py -- a RELABEL,
# a correction, a retirement. Pass a short reason so the log says why it ran:
#
#     sbatch rebuild_sapovirus.sh "MH003837 GNA4 -> GXV and the GI.P5/GI.P6 split"
#
# Replaces three near-identical scripts (rebuild_and_rebenchmark_sapovirus.sh,
# rebuild_after_gna_swap.sh, rebuild_after_gxv.sh) that differed only in their
# job name and header comment.
set -uo pipefail
REASON="${1:-no reason given}"
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
B=/home/people/fitzpatria/scratch/calicityper_runs/sapovirus_build
T="${SLURM_CPUS_PER_TASK:-16}"
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
echo "### reason: $REASON"

echo "### 0. regenerate the taxonomy tables from Supplemental Table 1"
python scripts/build_sapovirus_taxonomy.py

echo "### 1. rebuild reference sets"
python scripts/build_sapovirus_refsets.py --regions "$B/regions" \
    --recovered "$B/recovery/recovered_rdrp.fasta" --threads "$T" | tail -6

echo "### 2. rebuild the genogroup BLAST database"
python scripts/build_sapovirus_genogroup_db.py \
    --genbank "$B/sapovirus_with_prototypes.gb" --apply 2>&1 | tail -3

echo "### 3. diagnostics (within-type spread; is one member driving it?)"
python scripts/diagnose_sapovirus_rdrp.py 2>&1 | sed -n '/B. WITHIN-TYPE/,/C. WINDOW/p' | head -20

echo "### 4. nearest-neighbour vs cluster-mean"
python scripts/test_cluster_assignment.py 2>&1 | grep -vE "^$"

echo "### 5. truncation leave-one-out"
python scripts/truncation_loo_sapovirus.py --out results/truncation_loo_sapovirus.csv \
    --bootstraps 20 --threads "$T" 2>&1 | grep -E "references,|full-length|singletons|floor|^ 0\.|wrote"

echo "### 6. benchmark (published ranges + RIVM)"
python scripts/benchmark_sapovirus.py \
    --rivm-fasta /home/people/fitzpatria/scratch/calicityper_runs/rivm_genogroup/rivm_all.fasta \
    --out-dir results/benchmark_sapovirus --loo results/truncation_loo_sapovirus.csv \
    --threads "$T" 2>&1 | grep -vE "^\[type\]|^Loading"

echo "### 7. cross-check against both published tables"
python scripts/check_against_zhao_tables.py 2>&1 | sed -n "/SUMMARY/,999p"

echo "### DONE $(date)"
