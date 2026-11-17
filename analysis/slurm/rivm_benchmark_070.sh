#!/usr/bin/env bash
#SBATCH --job-name=noro_bench070
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=600
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/noro_bench070_%j.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
set -euo pipefail   # -e: a failed step must fail the job (1110626 crashed in typing but showed COMPLETED)
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
# BLAS pinned to 1 thread per worker: 62.9x faster, byte-identical calls (job 1100271, 2026-11-10).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export NCBI_EMAIL="calicityper@example.org"
cd /home/people/fitzpatria/scratch/calicityper
W=/home/people/fitzpatria/scratch/calicityper_runs/rivm_genogroup_0.7.0
mkdir -p "$W" results
FA="$W/rivm_all.fasta"
# seed from the already-fetched 4,392 so only the ~1,594 extra are downloaded
[ -s "$FA" ] || cp /home/people/fitzpatria/scratch/calicityper_runs/rivm_download/rivm_test_set.fasta "$FA" 2>/dev/null || true
echo "### host $(hostname)  $(date)"
echo "### step 1: fetch ALL 5,986 RIVM records (resumable)"
python3 scripts/fetch_rivm_test_set.py "$FA" --all --batch 50
N=$(grep -c "^>" "$FA" 2>/dev/null || echo 0); echo "sequences: $N"
[ "$N" -lt 5000 ] && { echo "FATAL: only $N fetched, expected ~5986"; exit 1; }
echo "### step 2: typing  $(date)"
# No --report: one alignment panel per call would make a ~240 MB HTML for
# 5,986 records. The TSV carries every call and novelty flag.
calicityper type "$FA" -o "$W/all_typed.tsv" --threads "${SLURM_CPUS_PER_TASK:-16}"
echo "typed rows: $(( $(wc -l < "$W/all_typed.tsv") - 1 ))"
echo "### step 3: compare, genogroup + genotype  $(date)"
python3 scripts/validate_vs_rivm.py "$W/all_typed.tsv" --out results/rivm_comparison_all_0.7.0.csv
echo "### summary"; cat results/rivm_comparison_all_0.7.0.csv
echo "### DONE $(date)"
