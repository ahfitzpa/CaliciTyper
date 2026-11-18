#!/usr/bin/env python3
"""F3: choose up to two references per genotype from DISTINCT intra-genotypic
clusters, and compare with the current 'first two in the alignment' rule.

Reads inputs/intragenotypic_clusters.tsv (genotype, cluster, accession). For
each genotype, picks one representative from each distinct cluster in stable
accession order, capped at two. Emits the proposed table plus a diff against
the current rule, for review. Does not change align.py or tree.py.
"""
from common import INPUTS, RESULTS, read_tsv, write_tsv


def main():
    rows = read_tsv(INPUTS / "intragenotypic_clusters.tsv",
                    required=["genotype", "cluster", "accession"])
    by_geno = {}
    for r in rows:
        if r["accession"].startswith("FILL"):
            continue
        by_geno.setdefault(r["genotype"], []).append(r)

    proposed, current = [], []
    for geno in sorted(by_geno):
        items = sorted(by_geno[geno], key=lambda r: (r["cluster"], r["accession"]))
        seen_clusters, picked = set(), []
        for r in items:
            if r["cluster"] in seen_clusters:
                continue
            seen_clusters.add(r["cluster"])
            picked.append(r)
            if len(picked) == 2:
                break
        for r in picked:
            proposed.append({"genotype": geno, "cluster": r["cluster"],
                             "accession": r["accession"], "rule": "distinct_cluster"})
        for r in items[:2]:
            current.append({"genotype": geno, "cluster": r["cluster"],
                            "accession": r["accession"], "rule": "first_two"})

    pkeys = {(r["genotype"], r["accession"]) for r in proposed}
    ckeys = {(r["genotype"], r["accession"]) for r in current}
    for r in proposed:
        r["changed_vs_current"] = int((r["genotype"], r["accession"]) not in ckeys)
    for r in current:
        r["changed_vs_current"] = int((r["genotype"], r["accession"]) not in pkeys)

    write_tsv(RESULTS / "f3_proposed_representatives.tsv", proposed,
              ["genotype", "cluster", "accession", "rule", "changed_vs_current"])
    write_tsv(RESULTS / "f3_current_representatives.tsv", current,
              ["genotype", "cluster", "accession", "rule", "changed_vs_current"])
    changed = sum(r["changed_vs_current"] for r in proposed)
    print(f"{len(proposed)} proposed representatives across {len(by_geno)} genotypes; "
          f"{changed} differ from the current first-two rule.")
    if not proposed:
        print("No usable rows (all accessions are FILL placeholders) — fill intragenotypic_clusters.tsv first.")


if __name__ == "__main__":
    main()
