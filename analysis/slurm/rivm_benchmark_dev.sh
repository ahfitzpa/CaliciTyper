#!/usr/bin/env bash
#SBATCH --job-name=noro_benchdev
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=600
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/noro_benchdev_%j.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
#
# Norovirus regression benchmark for the current working tree, against 0.7.0.
# Same 5,986 RIVM records and same comparison as rivm_benchmark_070.sh, but
# with no download (the 0.7.0 FASTA is reused) and a row-by-row diff against
# the 0.7.0 typing output. Checks two things the TODO asks for:
#   - the rebuilt calicityper_env (Python 3.15, numpy 2.5, IQ-TREE 3) gives
#     the benchmarked numbers, and
#   - the 2026-11-16 changes (sapovirus genogroups in the BLAST db, pinned
#     -dbsize, 100 nt minimum genogroup hit) move no norovirus GI/GII call.
# Default settings: sapovirus genotyping stays off, as for a user.
#
#     sbatch rivm_benchmark_dev.sh
set -euo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
# BLAS pinned to 1 thread per worker: 62.9x faster, byte-identical calls (job 1100271).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /home/people/fitzpatria/scratch/calicityper
OLD=/home/people/fitzpatria/scratch/calicityper_runs/rivm_genogroup_0.7.0
W=/home/people/fitzpatria/scratch/calicityper_runs/rivm_genogroup_dev
mkdir -p "$W" results
FA="$OLD/rivm_all.fasta"
echo "### host $(hostname)  $(date)"
echo "### commit $(git rev-parse --short HEAD)$(git diff --quiet || echo ' (+uncommitted changes)')"
echo "### sequences: $(grep -c '^>' "$FA")"

echo "### step 1: typing  $(date)"
calicityper type "$FA" -o "$W/all_typed.tsv" --threads "${SLURM_CPUS_PER_TASK:-16}"
echo "typed rows: $(( $(wc -l < "$W/all_typed.tsv") - 1 ))"

echo "### step 2: compare with RIVM  $(date)"
python3 scripts/validate_vs_rivm.py "$W/all_typed.tsv" --out results/rivm_comparison_all_dev.csv
echo "--- dev";   cat results/rivm_comparison_all_dev.csv
echo "--- 0.7.0"; cat /home/people/fitzpatria/scratch/calicityper_runs/benchmarking_wip_0.7.0/results/rivm_comparison_all_0.7.0.csv

echo "### step 3: row-by-row diff against 0.7.0  $(date)"
python3 - "$OLD/all_typed.tsv" "$W/all_typed.tsv" "$W/changed_vs_0.7.0.tsv" <<'EOF'
import sys
import pandas as pd
old = pd.read_csv(sys.argv[1], sep="\t", dtype=str).fillna("").set_index("query")
new = pd.read_csv(sys.argv[2], sep="\t", dtype=str).fillna("").set_index("query")
cols = ["genogroup", "orf1_type", "orf2_type", "orf1_accepted", "orf2_accepted"]
both = old.index.intersection(new.index)
print(f"queries: 0.7.0 {len(old)}, dev {len(new)}, shared {len(both)}")
rows = []
for q in both:
    ch = {c: (old.at[q, c], new.at[q, c]) for c in cols
          if c in old.columns and c in new.columns and old.at[q, c] != new.at[q, c]}
    if ch:
        rows.append(dict(query=q, **{f"{c}_0.7.0": a for c, (a, b) in ch.items()},
                         **{f"{c}_dev": b for c, (a, b) in ch.items()}))
d = pd.DataFrame(rows)
d.to_csv(sys.argv[3], sep="\t", index=False)
noro = lambda g: g in ("GI", "GII")
print(f"queries with any change: {len(d)}")
for c in cols:
    a, b = f"{c}_0.7.0", f"{c}_dev"
    if a in d.columns:
        sub = d[d[a].notna() & (d[a] != d[b])]
        gi_gii = sub[[noro(old.at[q, 'genogroup']) or noro(new.at[q, 'genogroup']) for q in sub['query']]]
        print(f"  {c}: {len(sub)} changed, {len(gi_gii)} of them GI/GII")
        for _, r in gi_gii.head(15).iterrows():
            print(f"      {r['query']}: {r[a]!r} -> {r[b]!r}")
print(f"wrote {sys.argv[3]}")
EOF
echo "### DONE $(date)"
