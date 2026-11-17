#!/usr/bin/env bash
#SBATCH --job-name=zen_fetch
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --time=60
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/zen_fetch_%j.out
#
# Fetch Zenodo 7271067 (Fitzpatrick et al. 2022, "Benchmarking bioinformatic
# tools for amplicon-based sequencing of norovirus", CC BY 4.0), the simulated
# amplicon benchmark with KNOWN expected composition, to replace the inherited
# RIVM panel as the primary norovirus capsid benchmark.
set -uo pipefail
D=/home/people/fitzpatria/scratch/calicityper_runs/zenodo_7271067
mkdir -p "$D"; cd "$D"
echo "### $(hostname) $(date)"
if [ ! -s record.json ]; then
    curl -sSL -o record.json "https://zenodo.org/api/records/7271067"
fi
python3 - <<'PY'
import json, urllib.request, os, hashlib
rec = json.load(open("record.json"))
print("title:", rec["metadata"]["title"])
print("doi  :", rec.get("doi"), "| licence:", rec["metadata"].get("license"))
for f in rec["files"]:
    name = f.get("key") or f.get("filename")
    url = f["links"].get("self") or f["links"].get("download")
    if os.path.exists(name) and os.path.getsize(name) == f.get("size", -1):
        continue
    urllib.request.urlretrieve(url, name)
    print(f"  fetched {name} ({os.path.getsize(name)} bytes)")
PY
echo "### files:"; ls -la | head -32
echo "### expected_sequences.001.fasta head:"; head -4 expected_sequences.001.fasta
echo "### records per simulation FASTA:"; for f in expected_sequences.*.fasta; do printf "%s %s\n" "$f" "$(grep -c '^>' $f)"; done
echo "### expected_composition.001.tsv head:"; head -3 expected_composition.001.tsv | cut -c1-300
echo "### amplicon length distribution (sim 001):"
python3 -c "
from itertools import groupby
import statistics as st
L=[]
for k,g in groupby(open('expected_sequences.001.fasta'), lambda l: l.startswith('>')):
    if not k: L.append(len(''.join(x.strip() for x in g)))
print('n',len(L),'min',min(L),'median',int(st.median(L)),'max',max(L))"
echo "### DONE $(date)"
