#!/usr/bin/env python3
"""Recover the Parra et al. immunotypes by rebuilding their clustering.

Their per-genotype immunotype assignments (A-L) appear only in Figure 5, an
image. S3 Table lists the exact ORF2 sequences that figure was built from,
so the assignments can be recovered by redoing the clustering on those
sequences rather than read off a picture or guessed.

Method: a neighbour-joining tree over the VP1 amino-acid sequences, cut into
twelve clades. The paper's Figure 5 legend says "Phylogenetic tree showing
the presence of 7 clusters containing two or more genotypes and 5 comprised
of a single genotype", so the immunotypes are clades, not groups formed by
thresholding a distance.

A distance threshold was tried first and does not work: no cut-off from 0.20
to 0.40 yields twelve clusters, and the structure jumps 17 -> 15 -> 8 -> 4
as single linkage chains through intermediate genotypes. The stated ">=20%
aa differences between groups" describes the separation the clades turn out
to have; it is not the rule that forms them. Mean distance within the
published immunotypes runs to about 27%, which is why thresholding at 20%
splits every one of them.

The result is checked against what the paper states in words before being
believed: immunotype G is GII.4 with GII.20; GI.3, GI.7, GI.8 and GI.9 form
one immunotype; there are 12 immunotypes, 7 of two or more genotypes and 5
singletons. If those do not come out, the reconstruction is reported as
failed rather than written anywhere.

    python scripts/rebuild_immunotype_tree.py [--cutoff 0.20]
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from collections import defaultdict

import numpy as np
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from common import INPUTS, RESULTS, write_tsv

FAA = INPUTS / "parra_s3" / "vp1_proteins.faa"
ALN = INPUTS / "parra_s3" / "vp1_proteins.aln.faa"

# what the paper says in words; the reconstruction has to reproduce these
ANCHOR_TOGETHER = [("GII.4", "GII.20"), ("GI.3", "GI.7"), ("GI.3", "GI.8"),
                   ("GI.3", "GI.9"), ("GI.7", "GI.8"), ("GI.7", "GI.9"),
                   ("GI.8", "GI.9")]
EXPECT_TOTAL, EXPECT_MULTI, EXPECT_SINGLE = 12, 7, 5


def align():
    if ALN.exists() and ALN.stat().st_mtime > FAA.stat().st_mtime:
        return ALN
    print("  aligning with mafft --auto --anysymbol --thread 1 ...", flush=True)
    r = subprocess.run(["mafft", "--quiet", "--auto", "--anysymbol", "--thread", "1",
                        str(FAA)], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"mafft failed:\n{r.stderr[-500:]}")
    ALN.write_text(r.stdout)
    return ALN


def pdist(a, b, min_sites=150, default=float("nan")):
    ok = [(x, y) for x, y in zip(a, b) if x not in "-X" and y not in "-X"]
    return (sum(x != y for x, y in ok) / len(ok)) if len(ok) >= min_sites else default


def nj_clades(recs, n_groups):
    """(clades as genotype lists, monophyly count) from a midpoint-rooted NJ tree.

    The tree is split by repeatedly breaking whichever current clade is the
    most internally divergent, which is what reading groups off a drawn tree
    amounts to, until the requested number of groups exists.
    """
    names = [n for n, _ in recs]
    m = [[0.0] * (i + 1) for i in range(len(recs))]
    for i in range(len(recs)):
        for j in range(i):
            m[i][j] = pdist(recs[i][1], recs[j][1], default=1.0)
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, m))
    tree.root_at_midpoint()

    gt = lambda n: n.split("|")[1]
    by_gt = defaultdict(list)
    for lf in tree.get_terminals():
        by_gt[gt(lf.name)].append(lf)
    mono = sum(1 for lv in by_gt.values()
               if len(tree.common_ancestor(lv).get_terminals()) == len(lv))

    front = [tree.root]
    while len(front) < n_groups:
        splittable = [c for c in front if c.clades]
        if not splittable:
            break
        pick = max(splittable,
                   key=lambda c: max(tree.distance(c, t) for t in c.get_terminals()))
        front.remove(pick)
        front.extend(pick.clades)
    groups = [sorted({gt(t.name) for t in c.get_terminals()}) for c in front]
    return sorted(groups, key=lambda g: (-len(g), g[0])), mono, len(by_gt)


# the four immunotype letters the paper names in words; the rest are legible
# only in the figure, so clusters recovered here carry no letter rather than
# an invented one
KNOWN_LETTERS = {tuple(sorted(k)): v for k, v in {
    ("GII.4", "GII.20"): "G",
    ("GI.3", "GI.7", "GI.8", "GI.9"): "C",
    ("GII.6",): "H"}.items()}
LETTER_NOTE = {"J": "GII.17 is named J in the text; its clade here also holds "
                    "GII.13 and GII.21"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--groups", type=int, default=EXPECT_TOTAL)
    a = ap.parse_args()

    if not FAA.exists():
        sys.exit(f"missing {FAA}; run scripts/fetch_s3_sequences.py first")
    aln = align()
    recs = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(aln, "fasta")]
    groups, mono, n_gt = nj_clades(recs, a.groups)
    multi = [g for g in groups if len(g) > 1]
    single = [g for g in groups if len(g) == 1]

    print(f"{len(recs)} proteins, {n_gt} genotypes, {mono} of them monophyletic\n")
    print(f"{len(groups)} clades: {len(multi)} with two or more genotypes, "
          f"{len(single)} singletons")
    print(f"paper states: {EXPECT_TOTAL}, {EXPECT_MULTI}, {EXPECT_SINGLE}\n")
    for i, g in enumerate(groups, 1):
        letter = KNOWN_LETTERS.get(tuple(g), "")
        print(f"  clade {i:2d} {letter or '-':1s}  {', '.join(g)}")

    print("\nchecks against what the paper states in words:")
    ok = True
    for x, y in ANCHOR_TOGETHER:
        together = any(x in g and y in g for g in groups)
        ok &= together
        print(f"  {x:7s} with {y:7s}  {'yes' if together else 'NO -- not reproduced'}")
    shape = (len(groups), len(multi), len(single)) == (EXPECT_TOTAL, EXPECT_MULTI,
                                                       EXPECT_SINGLE)
    print(f"  clade counts {len(groups)}/{len(multi)}/{len(single)}: "
          f"{'match' if shape else 'DO NOT match'}")
    print(f"  every genotype monophyletic: {'yes' if mono == n_gt else f'no ({mono}/{n_gt})'}")

    if not (ok and shape):
        print("\nDoes not reproduce the published structure; nothing written.")
        return 1

    rows = []
    for i, g in enumerate(groups, 1):
        letter = KNOWN_LETTERS.get(tuple(g), "")
        for gt_ in g:
            rows.append({"genotype": gt_, "clade": i, "immunotype": letter,
                         "n_genotypes": len(g), "members": " ".join(g),
                         "letter_source": ("Parra 2017 text" if letter
                                           else "not named in the text; see Fig 5")})
    out = RESULTS / "f4_immunotypes_reconstructed.tsv"
    write_tsv(out, rows, ["genotype", "clade", "immunotype", "n_genotypes",
                          "members", "letter_source"])
    print(f"\nReproduced the published structure. Wrote {out}")
    print("The GROUPING is what this recovers, and it is what a user needs: whether")
    print("two genotypes share an immunotype. Only three letters are named in the")
    print("paper's text, so the rest are left blank rather than invented.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
