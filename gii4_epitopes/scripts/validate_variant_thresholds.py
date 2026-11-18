#!/usr/bin/env python3
"""F1/F2: LOO sweep of candidate distance thresholds, using the AA alignment.

For each reference, remove it, compute uncorrected p-distance to the nearest
remaining reference of every variant, and record whether the nearest variant
matches the held-out label. Sweep the candidate within/between percentages
from inputs/variant_thresholds.tsv and report accuracy and the calls that
would be declined as ambiguous.

This script produces EVIDENCE for choosing a threshold. It does not choose one
and does not alter the tool.
"""
import csv
import sys
from pathlib import Path

from common import INPUTS, RESULTS, p_distance, read_fasta, read_tsv, write_tsv

ALN = RESULTS / "gii4_vp1_aa_279.aln"
FASTA = RESULTS / "gii4_vp1_aa_279.fasta"


def load_index():
    return read_tsv(RESULTS / "gii4_variant_index.tsv", required=["accession", "variant"])


def main():
    path = ALN if ALN.exists() else FASTA
    if not path.exists():
        raise SystemExit(f"missing {path}; run build_gii4_aa_reference.py first")
    seqs = read_fasta(path)
    index = load_index()
    labels = {f"{r['accession']}|{r['variant']}": r["variant"] for r in index}
    seqs = {k: v for k, v in seqs.items() if k in labels}
    if len(seqs) < 2:
        raise SystemExit("need at least two aligned references")

    thresholds = {}
    for r in read_tsv(INPUTS / "variant_thresholds.tsv"):
        if r["param"] in ("variant_within_pct", "variant_between_pct"):
            thresholds[r["param"]] = float(r["value"])
    within = thresholds.get("variant_within_pct", 0.05)
    between = thresholds.get("variant_between_pct", 0.05)

    rows = []
    for held, held_seq in seqs.items():
        truth = labels[held]
        best_var, best_d, second_d = None, float("inf"), float("inf")
        for other, other_seq in seqs.items():
            if other == held:
                continue
            d, n = p_distance(held_seq, other_seq)
            v = labels[other]
            if d < best_d:
                second_d, best_d, best_var = best_d, d, v
            elif d < second_d and v != best_var:
                second_d = d
        if best_var is None:
            continue
        ambiguous = (best_d > within) or (second_d <= between and best_var != truth)
        rows.append({
            "held_out": held, "true_variant": truth, "nearest_variant": best_var,
            "nearest_distance": f"{best_d:.4f}", "second_distance": f"{second_d:.4f}",
            "correct": int(best_var == truth),
            "would_flag_divergent": int(best_d > within),
            "ambiguous_at_candidate_threshold": int(ambiguous),
        })
    write_tsv(RESULTS / "f1_variant_loo.csv", rows,
              ["held_out", "true_variant", "nearest_variant", "nearest_distance",
               "second_distance", "correct", "would_flag_divergent",
               "ambiguous_at_candidate_threshold"])
    n = len(rows)
    if n:
        print(f"LOO n={n}: correct={sum(r['correct'] for r in rows)}/{n}, "
              f"flagged_divergent={sum(r['would_flag_divergent'] for r in rows)}/{n}")
    print("EVIDENCE ONLY — review results before adopting any threshold or changing the tool.")


if __name__ == "__main__":
    main()
