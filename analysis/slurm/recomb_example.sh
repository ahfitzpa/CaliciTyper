#!/usr/bin/env bash
#SBATCH --job-name=recomb_example
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=180
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/recomb_example_%j.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
#
# Worked example for `calicityper recomb` (TODO section 4: `type` has one in
# examples/, `recomb` has none). Runs one norovirus and one sapovirus
# synthetic recombinant from tests/synthetic.py -- the planted breakpoint is
# known, so the example shows the tool recovering it -- with the CLI exactly
# as a user would. Outputs go to calicityper_runs/recomb_example/ for review;
# nothing is written into the repo. Copy the small files into examples/ once
# they have been looked at.
#
#     sbatch recomb_example.sh
set -euo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export CALICITYPER_OPENRDP=/home/people/fitzpatria/scratch/openrdp_env/bin/openrdp
cd /home/people/fitzpatria/scratch/calicityper
OUT=/home/people/fitzpatria/scratch/calicityper_runs/recomb_example
mkdir -p "$OUT"
echo "### host $(hostname)  $(date)"
echo "### commit $(git rev-parse --short HEAD)$(git diff --quiet || echo ' (+uncommitted changes)')"

for spec in "norovirus clear_recombinant" "sapovirus sapovirus_clear_recombinant"; do
  set -- $spec; virus=$1; name=$2
  aln="$OUT/$name.fasta"
  # the generator writes the alignment; the query is its last record, as in tests/test_recomb.py
  query=$(python3 - "$name" "$aln" <<'EOF'
import sys
sys.path.insert(0, "tests")
import synthetic
from Bio import SeqIO
getattr(synthetic, sys.argv[1])(sys.argv[2])
print(list(SeqIO.parse(sys.argv[2], "fasta"))[-1].id)
EOF
)
  echo "### $virus: $name, query $query  $(date)"
  calicityper recomb --alignment "$aln" --query "$query" --virus "$virus" \
      -o "$OUT/$name" --threads "${SLURM_CPUS_PER_TASK:-8}" --verbose
  ls -la "$OUT/$name"
  python3 - "$OUT/$name/recombination.json" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
res = d["results"] if isinstance(d, dict) and "results" in d else d
for r in (res if isinstance(res, list) else [res]):
    rec = r.get("recombination", r)
    print("status:", rec.get("status"), "| score:", rec.get("score"), "| confidence:", rec.get("confidence"))
    for e in rec.get("events", []):
        bp = (e.get("breakpoints") or [{}])[0]
        print(f"  {e.get('event_id')}: {bp.get('alignment_start')}-{bp.get('alignment_end')} "
              f"| {e.get('biological_region')} | {','.join(e.get('supporting_methods', []))}")
EOF
done
echo "### DONE $(date)"
