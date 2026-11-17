#!/usr/bin/env bash
#SBATCH --job-name=noro_bench
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=120
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/noro_bench_%j.out
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
cd /home/people/fitzpatria/scratch/calicityper
python - <<'PY'
from Bio import SeqIO
rs=list(SeqIO.parse("/home/people/fitzpatria/scratch/calicityper_runs/rivm_download/rivm_test_set.fasta","fasta"))
print("total fetched:",len(rs))
for n in (100,400):
    SeqIO.write(rs[:n], f"/tmp/noro_s{n}.fa","fasta")
PY
python -c "from calicityper import tree; print('HAVE_NUMBA:', tree.HAVE_NUMBA)"
for n in 100 400; do
  echo "### $n sequences, 16 threads"
  S=$SECONDS
  norotyper /tmp/noro_s$n.fa -o /tmp/noro_s$n.tsv --threads 16
  E=$((SECONDS-S))
  echo "  wall: ${E}s  -> $(python -c "print(f'{$E/$n:.2f}')") s/seq  -> projected 4391 seqs: $(python -c "print(f'{$E/$n*4391/60:.1f}')") min"
  echo "  rows: $(( $(wc -l < /tmp/noro_s$n.tsv) - 1 ))"
  python - <<PY
import csv,collections
r=list(csv.DictReader(open("/tmp/noro_s$n.tsv"),delimiter="\t"))
print("  genogroups:",dict(collections.Counter(x["genogroup"] for x in r)))
PY
done
echo "=== DONE ==="
