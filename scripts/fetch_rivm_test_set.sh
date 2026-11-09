#!/usr/bin/env bash
# Download the GenBank sequences that were typed by the RIVM tool (validation/results_rivm.csv) so norotyper can be
# compared with RIVM calls. Needs entrez-direct and internet access to NCBI. Output: rivm_test_set.fasta
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-rivm_test_set.fasta}"
python3 - <<'PY' > .rivm_accessions.txt
import pandas as pd
r = pd.read_csv("validation/results_rivm.csv")
t = r[r["capsid type"].notna() | r["polymerase type"].notna()]
print("\n".join(sorted(set(t.name.str.split(".").str[0]))))
PY
echo "accessions: $(wc -l < .rivm_accessions.txt)"
: > "$OUT"
split -l 200 -d .rivm_accessions.txt .rivm_batch_
for b in .rivm_batch_*; do
  efetch -db nuccore -id "$(paste -sd, "$b")" -format fasta >> "$OUT"; sleep 0.5
done
rm -f .rivm_batch_* .rivm_accessions.txt
echo "wrote $OUT ($(grep -c '>' "$OUT") sequences)"
