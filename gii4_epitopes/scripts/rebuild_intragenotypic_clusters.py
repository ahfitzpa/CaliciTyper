#!/usr/bin/env python3
"""Recover the GII.6, GII.12 and GII.17 intra-genotypic variants.

S4 Table gives the accessions but no cluster column: the GII.6 A/B/C
assignment is in Figure 2A and the GII.17 A-D in S2 Fig, both images.
Assigning accessions to clusters from strain names would be invention, so
the partition is rebuilt from the sequences instead, exactly as the
immunotypes were.

F3 needs the partition, not the labels. The requirement is "pick the two
reference representatives from distinct intra-genotypic clusters", which a
partition answers and a letter does not.

The paper supplies its own control, and this script uses it. From S1 Fig:
GII.12 forms three clusters on nucleotides but "No discrete variant was
detected when amino acids were used". So an amino-acid method that finds
structure in GII.12 is wrong, whatever it reports for the others. GII.6
should give three variants (Fig 2A) and GII.17 four (S2 Fig).

    python scripts/rebuild_intragenotypic_clusters.py
"""
from __future__ import annotations

import csv
import itertools
import subprocess
import sys
from collections import defaultdict

import numpy as np
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from common import INPUTS, RESULTS, write_tsv

S4 = INPUTS / "parra2017_s4_gii6_gii12_gii17.tsv"
FAA = INPUTS / "parra_s4" / "vp1_proteins.faa"
#: variants each genotype is published as having, on amino acids
PUBLISHED = {"GII.6": 3, "GII.17": 4, "GII.12": None}   # None = no discrete variant


def align(seqs, path):
    if path.exists() and path.stat().st_size:
        return [(r.id, str(r.seq).upper()) for r in SeqIO.parse(path, "fasta")]
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.with_suffix(".in.faa")
    raw.write_text("".join(f">{n}\n{s}\n" for n, s in seqs))
    r = subprocess.run(["mafft", "--quiet", "--auto", "--anysymbol", "--thread", "1",
                        str(raw)], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"mafft failed:\n{r.stderr[-400:]}")
    path.write_text(r.stdout)
    raw.unlink()
    return [(x.id, str(x.seq).upper()) for x in SeqIO.parse(path, "fasta")]


def pdist(a, b, min_sites=150, default=1.0):
    ok = [(x, y) for x, y in zip(a, b) if x not in "-X" and y not in "-X"]
    return (sum(x != y for x, y in ok) / len(ok)) if len(ok) >= min_sites else default


def clusters(recs, k):
    names = [n for n, _ in recs]
    m = [[0.0] * (i + 1) for i in range(len(recs))]
    for i in range(len(recs)):
        for j in range(i):
            m[i][j] = pdist(recs[i][1], recs[j][1])
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, m))
    tree.root_at_midpoint()
    front = [tree.root]
    while len(front) < k:
        cand = [c for c in front if c.clades]
        if not cand:
            break
        pick = max(cand, key=lambda c: max(tree.distance(c, t) for t in c.get_terminals()))
        front.remove(pick)
        front.extend(pick.clades)
    return [sorted(t.name for t in c.get_terminals()) for c in front]


def spread(recs):
    v = [pdist(a, b, default=float("nan")) for (_, a), (_, b)
         in itertools.combinations(recs, 2)]
    v = [x for x in v if x == x]
    return (float(np.mean(v)), float(np.max(v))) if v else (float("nan"),) * 2


def main():
    if not FAA.exists():
        sys.exit(f"missing {FAA}; fetch the S4 accessions first "
                 f"(scripts/fetch_s3_sequences.py has the pattern)")
    by_gt = defaultdict(list)
    for r in SeqIO.parse(FAA, "fasta"):
        by_gt[r.id.split("|")[1]].append((r.id, str(r.seq)))

    print("amino-acid divergence per genotype:")
    sp = {}
    for g in sorted(by_gt):
        recs = align(by_gt[g], INPUTS / "parra_s4" / f"{g}.aln.faa")
        sp[g] = spread(recs)
        print(f"  {g:7s} n={len(recs):3d}  mean {sp[g][0] * 100:4.1f}%  "
              f"max {sp[g][1] * 100:4.1f}%")

    ctrl = "GII.12"
    others = [g for g in sp if g != ctrl]
    if ctrl in sp and others:
        ok = all(sp[ctrl][1] < sp[o][1] for o in others)
        print(f"\ncontrol: {ctrl} max divergence below "
              f"{', '.join(others)}: {'yes' if ok else 'NO'}")
        print(f"  the paper reports no discrete {ctrl} variant on amino acids, so a")
        print(f"  method finding structure there is wrong")
        if not ok:
            print("\nControl fails; not writing a partition.")
            return 1

    rows = []
    for g, k in PUBLISHED.items():
        if g not in by_gt or not k:
            continue
        recs = align(by_gt[g], INPUTS / "parra_s4" / f"{g}.aln.faa")
        parts = sorted(clusters(recs, k), key=lambda c: (-len(c), c[0]))
        print(f"\n{g}: {len(parts)} clusters (published {k})")
        for i, c in enumerate(parts, 1):
            print(f"  cluster {i}: n={len(c)}  e.g. {', '.join(x.split('|')[0] for x in c[:4])}")
            for name in c:
                rows.append({"genotype": g, "cluster": f"{g}_c{i}",
                             "accession": name.split("|")[0],
                             "source_label": "",
                             "evidence": "NJ clade of Parra 2017 S4 Table VP1 proteins; "
                                         "letters are in Fig 2A / S2 Fig and are not used"})
    if rows:
        out = RESULTS / "f3_intragenotypic_clusters_reconstructed.tsv"
        write_tsv(out, rows, ["genotype", "cluster", "accession", "source_label", "evidence"])
        print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
