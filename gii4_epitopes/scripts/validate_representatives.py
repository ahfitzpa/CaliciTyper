#!/usr/bin/env python3
"""F3 validation: LOO with current vs proposed representatives.

For each reference in inputs/intragenotypic_clusters.tsv that has a same-type
partner, remove it and compute its p-distance to the retained representatives
under (a) the current first-two rule and (b) the proposed distinct-cluster
rule. The proposal is acceptable only if it does not reduce recovery or
increase the distance to the nearest retained representative. This is a
screening check; a full re-typing run on the HPC is still required.
"""
from pathlib import Path

from common import RESULTS, p_distance, read_fasta, read_tsv, write_tsv

CURRENT = RESULTS / "f3_current_representatives.tsv"
PROPOSED = RESULTS / "f3_proposed_representatives.tsv"


def main():
    if not CURRENT.exists() or not PROPOSED.exists():
        raise SystemExit("run select_representatives.py first")
    data_dir = None
    try:
        from common import check_data_dir
        data_dir = check_data_dir()
    except SystemExit as e:
        print(str(e), file="__import__('sys').stderr")

    vp1 = None
    if data_dir:
        cand = Path(data_dir) / "reference_sequences_VP1" / "GII.fasta"
        if cand.exists():
            vp1 = read_fasta(cand)
    if vp1 is None:
        raise SystemExit("Need verified reference set with GII VP1 FASTA for LOO distances")

    cur = read_tsv(CURRENT, required=["genotype", "accession"])
    prop = read_tsv(PROPOSED, required=["genotype", "accession"])
    clusters = read_tsv(__import__("common").INPUTS / "intragenotypic_clusters.tsv")
    acc_to_geno = {r["accession"]: r["genotype"] for r in clusters}

    def nearest(acc, reps):
        a = vp1.get(acc) or vp1.get(acc.split(".")[0])
        if a is None:
            return None
        ds = []
        for r in reps:
            if r["accession"] == acc:
                continue
            b = vp1.get(r["accession"]) or vp1.get(r["accession"].split(".")[0])
            if b is None:
                continue
            d, n = p_distance(a, b)
            if n >= 100:
                ds.append(d)
        return min(ds) if ds else None

    rows = []
    for acc, geno in acc_to_geno.items():
        a = nearest(acc, [r for r in cur if r["genotype"] == geno])
        b = nearest(acc, [r for r in prop if r["genotype"] == geno])
        if a is None or b is None:
            continue
        rows.append({"genotype": geno, "accession": acc,
                     "nearest_current": f"{a:.4f}", "nearest_proposed": f"{b:.4f}",
                     "delta": f"{b - a:+.4f}",
                     "proposed_worse": int(b > a + 1e-9)})
    write_tsv(RESULTS / "f3_representative_loo.csv", rows,
              ["genotype", "accession", "nearest_current", "nearest_proposed",
               "delta", "proposed_worse"])
    worse = sum(r["proposed_worse"] for r in rows)
    print(f"LOO per-reference n={len(rows)}; proposed increases nearest distance in {worse}/{len(rows)}. "
          "Any increase requires full re-typing before adoption.")


if __name__ == "__main__":
    main()
