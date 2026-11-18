#!/usr/bin/env python3
"""Leave-one-genotype-out evaluation of the putative-novel screen.

Three different questions have been confused in discussing this screen, and
this script answers the third.

  leave-one-out       removes one strain. Its genotype remains represented, so
                      the right answer is to recover the type, not to flag it.
  singleton test      removes the only strain of a type. Correct answer is a
                      decline, but it covers only types that happen to have
                      one member.
  leave-one-genotype-out (here)
                      removes EVERY member of a genotype, then screens those
                      members against a reference set from which their type is
                      genuinely absent. That is the situation a novel genotype
                      presents, and it can be set up for any multi-member type
                      rather than only for singletons.

Scoring, per round (one genotype removed):
  removed members  -> expect a novelty flag   (sensitivity)
  retained members -> expect no flag          (specificity)

Specificity is pooled over (round, retained sequence) pairs, so each retained
sequence contributes once per round; the denominator is reported so the figure
is not mistaken for a per-sequence rate.

The bounds are recomputed from the retained references in every round, so no
information about the removed genotype leaks into its own test.

    python scripts/logo_novelty.py --out results/logo_novelty.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import putative_novel as nv   # noqa: E402


def bounds_from(blocks, keep_idx, region):
    """The 2 x SD table rebuilt from only the sequences in keep_idx."""
    table = {}
    for level in ("genotype", "genogroup"):
        pairs, size = {}, {}
        key = (lambda m: m[2]) if level == "genotype" else (lambda m: m[1])
        for b, (meta, D) in enumerate(blocks):
            keep = [i for i in range(len(meta)) if (b, i) in keep_idx]
            for x, i in enumerate(keep):
                c = key(meta[i])
                pairs.setdefault(c, [])
                size[c] = size.get(c, 0) + 1
                for j in keep[x + 1:]:
                    if key(meta[j]) == c:
                        pairs[c].append(D[i, j])
        stats, (pm, psd) = nv.cluster_stats(pairs, size, nv.OWN_LEVEL[level])
        for c, s in stats.items():
            table[(region, level, c)] = s["bound"]
        table[(region, level, "*pooled*")] = pm + 2 * psd
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="results/logo_novelty.csv")
    ap.add_argument("--measure", default="ml", choices=["ml", "nj", "p"])
    ap.add_argument("--min-members", type=int, default=2,
                    help="smallest genotype to remove; 1 reproduces the singleton test")
    a = ap.parse_args()

    aa = dict(nv.vp1_aa_alignment()[0])
    _vp1 = nv.vp1_distance_matrix()
    blocks = {
        "orf2": [(_vp1[0], nv.as_patristic(_vp1[1], a.measure,
                                           records=[(m[0], aa[m[0]]) for m in _vp1[0]],
                                           alphabet="aa"))],
        "orf1": [(meta, nv.as_patristic(D, a.measure,
                                        records=nv.rdrp_alignment(meta[0][1]), alphabet="nt"))
                 for meta, D in nv.rdrp_distance_blocks()],
    }

    rows = []
    for region, bl in blocks.items():
        counts = Counter(m[2] for meta, _ in bl for m in meta)
        targets = sorted(t for t, n in counts.items() if n >= a.min_members)
        print(f"\n[{region}] {sum(counts.values())} references, {len(counts)} genotypes; "
              f"removing each of the {len(targets)} with >= {a.min_members} members")

        for t in targets:
            keep = {(b, i) for b, (meta, _) in enumerate(bl)
                    for i in range(len(meta)) if meta[i][2] != t}
            gone = [(b, i) for b, (meta, _) in enumerate(bl)
                    for i in range(len(meta)) if meta[i][2] == t]
            # a genotype whose removal empties its genogroup leaves nothing to
            # compare against, which is a coverage artefact, not a result
            if len(keep) < 10:
                continue
            table = bounds_from(bl, keep, region)

            for b, i in gone:
                meta, D = bl[b]
                k = [j for j in range(len(meta)) if (b, j) in keep]
                if not k:
                    continue
                d = D[i, k]
                flag = nv.classify(region,
                                   nv._averages(d, [meta[j][2] for j in k]),
                                   nv._averages(d, [meta[j][1] for j in k]), table)[0]
                rows.append(dict(region=region, removed_genotype=t, role="removed",
                                 accession=meta[i][0], genogroup=meta[i][1],
                                 flag=flag, flagged=flag in nv.FLAGGED,
                                 correct=flag in nv.FLAGGED))
            for b, (meta, D) in enumerate(bl):
                for i in range(len(meta)):
                    if (b, i) not in keep:
                        continue
                    k = [j for j in range(len(meta)) if (b, j) in keep and j != i]
                    if not k:
                        continue
                    d = D[i, k]
                    flag = nv.classify(region,
                                       nv._averages(d, [meta[j][2] for j in k]),
                                       nv._averages(d, [meta[j][1] for j in k]), table)[0]
                    rows.append(dict(region=region, removed_genotype=t, role="retained",
                                     accession=meta[i][0], genogroup=meta[i][1],
                                     flag=flag, flagged=flag in nv.FLAGGED,
                                     correct=flag not in nv.FLAGGED))

    print(f"\n{'=' * 72}\nLEAVE-ONE-GENOTYPE-OUT NOVELTY SCREEN\n{'=' * 72}")
    for region in ("orf2", "orf1"):
        rem = [r for r in rows if r["region"] == region and r["role"] == "removed"]
        ret = [r for r in rows if r["region"] == region and r["role"] == "retained"]
        if not rem:
            continue
        sens = sum(r["flagged"] for r in rem) / len(rem)
        spec = sum(not r["flagged"] for r in ret) / len(ret) if ret else float("nan")
        print(f"\n[{region}]  genotypes removed: "
              f"{len({r['removed_genotype'] for r in rem})}")
        print(f"  sensitivity (removed genotype flagged) : {sens:.3f}  "
              f"({sum(r['flagged'] for r in rem)}/{len(rem)} sequences)")
        print(f"  specificity (retained not flagged)     : {spec:.3f}  "
              f"({sum(not r['flagged'] for r in ret)}/{len(ret)} round-sequence pairs)")
        worst = Counter(r["removed_genotype"] for r in rem if not r["flagged"])
        if worst:
            print(f"  genotypes missed entirely or in part: "
                  f"{', '.join(f'{k} x{v}' for k, v in worst.most_common(8))}")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
