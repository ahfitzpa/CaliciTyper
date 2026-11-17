#!/usr/bin/env bash
#SBATCH --job-name=sapo_fetch
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --time=90
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/sapo_fetch_%j.out
set -euo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
cd /home/people/fitzpatria/scratch/calicityper
W=/home/people/fitzpatria/scratch/calicityper_runs/sapovirus_build
mkdir -p "$W"
python3 - <<'PY'
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path
W = Path("/home/people/fitzpatria/scratch/calicityper_runs/sapovirus_build")
accs = {}
for f, col in (("validation/sapovirus_vp1_taxonomy.csv", "vp1"),
               ("validation/sapovirus_rdrp_taxonomy.csv", "rdrp")):
    for line in open(f):
        if line.startswith("#"): continue
        p = [x.split("#")[0].strip() for x in line.strip().split(",")]
        if len(p) < 2 or p[1] in ("accession", ""): continue
        accs.setdefault(p[1], {})[col] = p[0]
accs = {k: v for k, v in accs.items() if k}
print("unique accessions:", len(accs))
json.dump(accs, open(W / "accession_labels.json", "w"), indent=2)

# GenBank flatfile, not FASTA: the ORF coordinates are needed and only the
# flatfile carries the CDS features.
out = W / "sapovirus.gb"
ids = sorted(accs)
with open(out, "w") as fh:
    for i in range(0, len(ids), 40):
        chunk = ids[i:i+40]
        q = {"db": "nuccore", "id": ",".join(chunk), "rettype": "gbwithparts",
             "retmode": "text", "tool": "calicityper", "email": "calicityper@example.org"}
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(q)
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=180) as r:
                    fh.write(r.read().decode("utf-8", "replace"))
                break
            except Exception as e:
                print(f"  retry {attempt+1}: {type(e).__name__}", flush=True)
                time.sleep(2 ** attempt)
        print(f"  [{i//40+1}] {len(chunk)} records", flush=True)
        time.sleep(0.4)
print("wrote", out)
PY
grep -c "^LOCUS" "$W/sapovirus.gb" && echo "records fetched"
echo "=== DONE ==="
