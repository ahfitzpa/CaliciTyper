#!/usr/bin/env python3
"""Compare norotyper output with the RIVM Norovirus Typing Tool results in validation/results_rivm.csv.
Usage: norotyper rivm_test_set.fasta -o norotyper_out.tsv --threads 8
       python scripts/validate_vs_rivm.py norotyper_out.tsv"""
import argparse, re
import pandas as pd

def norm(x):
    if pd.isna(x) or str(x).strip() in ("", "Could not assign", "unassigned"):
        return "unassigned"
    return re.sub(r"\s*\(.*\)$", "", str(x).strip())      # 'GI.P10 (GI.Pa)' -> 'GI.P10'

ap = argparse.ArgumentParser()
ap.add_argument("norotyper_tsv"); ap.add_argument("--rivm", default="validation/results_rivm.csv")
ap.add_argument("--out", default="results/rivm_comparison.csv")
a = ap.parse_args()
mine = pd.read_csv(a.norotyper_tsv, sep="\t", dtype=str)
mine["acc"] = mine["query"].str.split(".").str[0]
rv = pd.read_csv(a.rivm); rv["acc"] = rv["name"].str.split(".").str[0]
m = mine.merge(rv, on="acc", suffixes=("_nt", "_rivm"))
rows = []
for label, mc, rc, ms, rs in (("ORF2 capsid", "orf2_type", "capsid type", "orf2_support", "capsid type support"),
                              ("ORF1 P-type", "orf1_type", "polymerase type", "orf1_support", "polymerase type support")):
    d = m.copy()
    d["nt"] = d[mc].map(norm); d["rivm"] = d[rc].map(norm)
    d = d[(d.nt != "unassigned") | (d.rivm != "unassigned")]
    both = d[(d.nt != "unassigned") & (d.rivm != "unassigned")]
    row = dict(region=label, n=len(d), both_assigned=len(both),
               agree_of_both=round(100 * (both.nt == both.rivm).mean(), 2) if len(both) else None,
               rivm_only=int(((d.nt == "unassigned") & (d.rivm != "unassigned")).sum()),
               norotyper_only=int(((d.nt != "unassigned") & (d.rivm == "unassigned")).sum()),
               agree_overall=round(100 * (d.nt == d.rivm).mean(), 2))
    rows.append(row)
    d[d.nt != d.rivm][["query", "nt", "rivm", ms, rs]].to_csv(a.out.replace(".csv", f"_{label.split()[0]}_discordant.csv"), index=False)
s = pd.DataFrame(rows); s.to_csv(a.out, index=False); print(s.to_string(index=False))
