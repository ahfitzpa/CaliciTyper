#!/usr/bin/env python3
"""Leave-one-out check of the novelty screen (per-cluster 2x SD rule) on the
bundled references.

Each reference is screened as a query, with the cluster bounds recomputed
without it, in three scenarios:

  own genotype present   its genotype's other members remain   -> expect `known`
  own genotype removed   every member of its genotype removed  -> expect a flag
  own genogroup removed  (VP1 only) its whole genogroup removed -> expect
                         `putative_new_genogroup`

RdRp distances exist only within one genogroup's alignment, so the genogroup
scenario is VP1-only. Run time: about a minute.

    python tests/novelty_validation.py
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from calicityper import novelty as nv  # noqa: E402


def table_without(region, blocks, drop):
    """Bound table in bounds() format, built from references not in `drop`
    (a set of (block index, row index))."""
    t = {}
    for level in ("genotype", "genogroup"):
        pairs, size = {}, {}
        for b, (meta, D) in enumerate(blocks):
            keep = [i for i in range(len(meta)) if (b, i) not in drop]
            key = (lambda m: m[2]) if level == "genotype" else (lambda m: m[1])
            for x, i in enumerate(keep):
                pairs.setdefault(key(meta[i]), [])
                size[key(meta[i])] = size.get(key(meta[i]), 0) + 1
                for j in keep[x + 1:]:
                    if key(meta[j]) == key(meta[i]):
                        pairs[key(meta[i])].append(D[i, j])
        stats, (pm, psd) = nv.cluster_stats(pairs, size, nv.OWN_LEVEL[level])
        for c, s in stats.items():
            t[(region, level, c)] = s["bound"]
        t[(region, level, "*pooled*")] = pm + 2 * psd
    return t


def run(region, blocks, tally, genogroup_scenario):
    for b, (meta, D) in enumerate(blocks):
        for i, (acc, g, t) in enumerate(meta):
            scen = [("own genotype present", {(b, i)})]
            scen.append(("own genotype removed", {(b, j) for j, m in enumerate(meta) if m[2] == t}))
            if genogroup_scenario:
                scen.append(("own genogroup removed", {(b, j) for j, m in enumerate(meta) if m[1] == g}))
            for name, drop in scen:
                if name == "own genotype present" and sum(m[2] == t for m in meta) < 2:
                    continue
                keep = [j for j in range(len(meta)) if (b, j) not in drop]
                if not keep:
                    continue
                d = D[i, keep]
                type_avg = nv._averages(d, [meta[j][2] for j in keep])
                group_avg = nv._averages(d, [meta[j][1] for j in keep])
                flag = nv.classify(region, type_avg, group_avg, table_without(region, blocks, drop))[0]
                tally[(region, name)][flag] += 1


tally = defaultdict(Counter)
run("orf2", [nv.vp1_distance_matrix()], tally, genogroup_scenario=True)
run("orf1", list(nv.rdrp_distance_blocks()), tally, genogroup_scenario=False)

order = [nv.KNOWN, nv.NEW_TYPE, nv.NEW_GROUP]
print("region\tscenario\tn\t" + "\t".join(order))
for (region, scen), c in tally.items():
    n = sum(c.values())
    print(f"{region}\t{scen}\t{n}\t" + "\t".join(f"{c[o]} ({c[o] / n:.0%})" for o in order))
