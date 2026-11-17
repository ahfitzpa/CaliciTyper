#!/usr/bin/env bash
#SBATCH --job-name=noro_ab
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=120
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/noro_ab_%j.out
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
cd /home/people/fitzpatria/scratch/calicityper
python - <<'PY'
from Bio import SeqIO
rs=list(SeqIO.parse("/home/people/fitzpatria/scratch/calicityper_runs/rivm_download/rivm_test_set.fasta","fasta"))
SeqIO.write(rs[:100],"/tmp/noro_ab100.fa","fasta")
PY
echo "### A: default BLAS threading (what the slow runs do)"
unset OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS NUMEXPR_NUM_THREADS
S=$SECONDS; calicityper type /tmp/noro_ab100.fa -o /tmp/noro_ab_A.tsv --threads 16; A=$((SECONDS-S))
echo "  A wall: ${A}s"

echo "### B: BLAS pinned to 1 thread per worker"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
S=$SECONDS; calicityper type /tmp/noro_ab100.fa -o /tmp/noro_ab_B.tsv --threads 16; B=$((SECONDS-S))
echo "  B wall: ${B}s"

echo "### result"
python - <<PY
a,b=$A,$B
print(f"  default BLAS : {a}s")
print(f"  pinned BLAS  : {b}s")
print(f"  speedup      : {a/b:.1f}x" if b else "  n/a")
print(f"  projected 4391 seqs, pinned: {b/100*4391/60:.1f} min")
PY
echo "### identical calls?"
if diff -q /tmp/noro_ab_A.tsv /tmp/noro_ab_B.tsv >/dev/null; then echo "  YES - byte-identical"; else echo "  NO - DIFFER:"; diff /tmp/noro_ab_A.tsv /tmp/noro_ab_B.tsv | head; fi
echo "=== DONE ==="
