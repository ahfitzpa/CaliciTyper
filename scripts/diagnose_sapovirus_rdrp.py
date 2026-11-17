#!/usr/bin/env python3
"""Why do the sapovirus RdRp reference sets miss the published distance ranges?

Three candidate explanations, tested separately:

 A. COMPOSITION -- our sets hold a different set of strains per P-type than the
    paper used, because only records carrying a complete window were admitted.
 B. OUTLIER MEMBERS -- a single misassigned or misaligned strain inflates the
    within-type maximum, rather than the type being genuinely diffuse.
 C. WINDOW -- our RdRp region is not the paper's. The VP1 boundaries are taken
    from annotation; the RdRp 5' boundary is DERIVED (no sapovirus GenBank
    record annotates NS7), so a shifted or over-long window would inflate
    within-type distances and erode separation between types. VP1, whose
    boundary is annotated, reproduces the published ranges; RdRp does not.
    That contrast is what makes this the leading hypothesis.

Test C scans sub-windows of the existing RdRp alignment and asks whether any
of them reproduces the published within-type ranges. If one does, the region
definition is wrong and can be corrected without new sequence data. If none
does, the problem is the strains themselves (A/B) or the published cut-off.

    python scripts/diagnose_sapovirus_rdrp.py
"""
from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import numpy as np
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import config as cg, tree   # noqa: E402

REGIONS = {"RdRp": ("orf1", ROOT / "validation" / "sapovirus_rdrp_taxonomy.csv", "ptype"),
           "VP1": ("orf2", ROOT / "validation" / "sapovirus_vp1_taxonomy.csv", "genotype")}


def published(path, key):
    rows, rng = [], {}
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        if r["note"].startswith(("UNRESOLVED", "RETIRED")):
            continue
        rows.append((r[key], r["accession"].split(".")[0]))
        t = r["pairwise_distance_range"].strip()
        if t and t != "NA":
            lo, _, hi = t.partition("-")
            rng[r[key]] = (float(lo), float(hi or lo))
    return rows, rng


def load(region):
    """{type: [(accession, coded seq)]} over the per-genogroup alignments."""
    out = {}
    d = cg.SAPOVIRUS_REF_DIR[REGIONS[region][0]]
    for orf in sorted(d.glob("*.orf")):
        if orf.stem == "all":
            continue
        types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
        aln = AlignIO.read(orf.with_suffix(".aln"), "clustal")
        codes = tree.encode([str(r.seq).upper() for r in aln])
        for i, rec in enumerate(aln):
            out.setdefault(types[rec.id], []).append((rec.id, codes[i]))
    return out


def pdist(a, b, cols=None):
    if cols is not None:
        a, b = a[cols], b[cols]
    both = (a < 4) & (b < 4)
    n = int(both.sum())
    return (float((a[both] != b[both]).sum()) / n) if n >= 50 else None


def section(title):
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def main():
    for region in ("RdRp", "VP1"):
        _, tax, key = REGIONS[region]
        rows, rng = published(tax, key)
        members = load(region)
        pub_types = {t for t, _ in rows}
        pub_accs = {t: {a for tt, a in rows if tt == t} for t in pub_types}

        section(f"{region}: A. COMPOSITION -- published strains vs strains we hold")
        have_t = set(members)
        print(f"published: {len(pub_types)} types, {len(rows)} strains")
        print(f"ours     : {len(have_t)} types, {sum(len(v) for v in members.values())} strains")
        print(f"types with NO reference at all: {len(pub_types - have_t)}")
        part = [(t, len(members[t]), len(pub_accs[t])) for t in sorted(have_t)
                if len(members[t]) < len(pub_accs[t])]
        print(f"types we hold but with fewer strains than published: {len(part)}")
        for t, got, want in part[:10]:
            print(f"    {t:10s} {got}/{want}")

        section(f"{region}: B. WITHIN-TYPE SPREAD, and whether one member drives it")
        bad = []
        for t in sorted(members):
            if t not in rng or len(members[t]) < 2:
                continue
            names = [n for n, _ in members[t]]
            seqs = [s for _, s in members[t]]
            ds = {}
            for i, j in itertools.combinations(range(len(seqs)), 2):
                v = pdist(seqs[i], seqs[j])
                if v is not None:
                    ds[(i, j)] = v
            if not ds:
                continue
            mx = max(ds.values())
            lo, hi = rng[t]
            if mx <= hi + 0.005:
                continue
            # drop each member in turn: does one member explain the excess?
            best = None
            for k in range(len(seqs)):
                rest = [v for (i, j), v in ds.items() if i != k and j != k]
                if rest and (best is None or max(rest) < best[1]):
                    best = (names[k], max(rest))
            bad.append((t, mx, hi, best))
            fix = "" if best is None else (
                f"  -> without {best[0]}: {best[1]:.3f}"
                + ("  [INSIDE published range]" if best[1] <= hi + 0.005 else ""))
            print(f"  {t:10s} ours max {mx:.3f} vs published {lo}-{hi}, n={len(seqs)}{fix}")
        if not bad:
            print("  every type with a published range is inside it")

        if region != "RdRp":
            continue

        section("RdRp: C. WINDOW -- does a sub-window reproduce the published ranges?")
        # alignments differ in width between genogroup files, so scan within each
        widths = {len(v) for m in members.values() for _, v in m}
        width = min(widths)
        print(f"aligned widths across genogroup files: {sorted(widths)}")
        print(f"current aligned RdRp window: {width} columns "
              f"(region defined as 512 aa / 1536 nt by projection; 5' boundary derived)")
        scored = []
        for frac in (1.0, 0.9, 0.75, 0.6, 0.5, 0.4, 0.3):
            w = int(width * frac)
            for anchor in ("5prime", "centre", "3prime"):
                start = {"5prime": 0, "centre": (width - w) // 2, "3prime": width - w}[anchor]
                cols0 = np.arange(start, start + w)
                inside = tot = 0
                sep_ok = sep_tot = 0
                for t in members:
                    if len(members[t]) < 2:
                        continue
                    seqs = [s for _, s in members[t]]
                    cols = cols0[cols0 < len(seqs[0])]
                    ds = [pdist(a, b, cols) for a, b in itertools.combinations(seqs, 2)]
                    ds = [v for v in ds if v is not None]
                    if not ds:
                        continue
                    if t in rng:
                        tot += 1
                        inside += max(ds) <= rng[t][1] + 0.005
                    # nearest other type, within the same genogroup file
                    others = [pdist(a, b, cols[cols < min(len(a), len(b))])
                              for o in members if o != t
                              and o.split(".")[0] == t.split(".")[0]
                              for a in seqs for _, b in members[o]]
                    others = [v for v in others if v is not None]
                    if others:
                        sep_tot += 1
                        sep_ok += min(others) >= cg.SAPOVIRUS_CUTOFFS["orf1"][0]
                if tot:
                    scored.append((frac, anchor, w, inside, tot, sep_ok, sep_tot))
        print(f"{'window':>8} {'anchor':>8} {'cols':>6}  {'within published':>17}  {'separated at 0.266':>19}")
        for frac, anchor, w, inside, tot, sep_ok, sep_tot in scored:
            flag = "  <-- best" if tot and inside == max(s[3] for s in scored if s[4] == tot) else ""
            print(f"{frac:8.2f} {anchor:>8} {w:6d}  {inside:7d}/{tot:<9d}  {sep_ok:9d}/{sep_tot:<9d}{flag}")
        print("\nIf no window reproduces the published ranges, the region definition is not")
        print("the explanation and the strain composition or the cut-off itself is.")

        section("RdRp: D. ARE THE OUTLIER STRAINS SIMPLY IN THE WRONG P-TYPE?")
        for t, mx, hi, best in bad:
            if best is None:
                continue
            name = best[0]
            seq = next(s for n_, s in members[t] if n_ == name)
            own = [pdist(seq, s2) for n2, s2 in members[t] if n2 != name]
            own = [v for v in own if v is not None]
            near = []
            for o in members:
                # only within the same genogroup: each genogroup has its own
                # alignment, so columns are not homologous across files
                if o == t or o.split(".")[0] != t.split(".")[0]:
                    continue
                ds = [pdist(seq, s2) for _, s2 in members[o]]
                ds = [v for v in ds if v is not None]
                if ds:
                    near.append((min(ds), o))
            near.sort()
            cut = cg.SAPOVIRUS_CUTOFFS["orf1"][0]
            print(f"  {name} (labelled {t}, published range <= {hi})")
            print(f"    distance to its own type : {min(own):.3f}-{max(own):.3f}")
            for v, o in near[:3]:
                mark = "  <-- within the P-type cut-off of another type" if v < cut else ""
                print(f"    nearest {o:10s}        : {v:.3f}{mark}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
