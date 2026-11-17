#!/usr/bin/env python3
"""Is sapovirus VP1 genotype GI.5 correctly assigned P-type GI.P5?

Background. Zhao et al. 2024 name the GI.P genogroup GI.P1-GI.P3 and
GI.P5-GI.P7: there is no GI.P4, because GI.P4 and GI.P1 fell into one cluster
under the 2 x SD criterion and were merged as GI.P1. Their Table 2 gives
DQ366345 (Ehime643) as GI.5 with P-type GI.P5, i.e. the dual type GI.5[P5].

Our own reference sets disagree with the published range for GI.P5: the
within-type maximum is 0.255 against a published 0.052, and both DQ366345 and
LC504396 sit nearer GI.P3 than to each other under a cluster-mean rule. This
asks whether that is a real problem and, if so, which of the possible causes
it is.

Method, following the published framework rather than inventing one:
  * take every GI.P reference we hold, with its published P-type;
  * compute all pairwise RdRp nucleotide distances over the aligned region;
  * apply the 2 x SD criterion per cluster (mean within-cluster distance plus
    two standard deviations) exactly as the putative-novel screen does;
  * ask whether the GI.5-associated RdRp sequences form a cluster distinct
    from every established GI.P type, and in particular whether GI.P1 and
    GI.P5 are separable at all now that GI.P4 has been merged into GI.P1.

No new P-type is proposed from a single unnamed sequence: a lone unnamed
sequence would be provisional (GI.PNA[n]), and a previously named sequence
keeps its published name even as a singleton.

    python scripts/analyse_gi5_ptype.py --out-dir results/gi5_ptype
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import config as cg, tree   # noqa: E402

REFERENCE = "DQ366345"          # Ehime643, the published GI.5 / GI.P5 representative
GI_P_TYPES = ("GI.P1", "GI.P2", "GI.P3", "GI.P5", "GI.P6", "GI.P7")


def pdist(a, b):
    both = (a < 4) & (b < 4)
    n = int(both.sum())
    return (float((a[both] != b[both]).sum()) / n) if n >= 50 else None


def load_group(region_key, stem):
    d = cg.SAPOVIRUS_REF_DIR[region_key]
    orf, aln = d / f"{stem}.orf", d / f"{stem}.aln"
    if not orf.exists() or not aln.exists():
        return [], None
    types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
    a = AlignIO.read(aln, "clustal")
    names = [r.id.split("|")[-1] for r in a]
    codes = tree.encode([str(r.seq).upper() for r in a])
    return [(names[i], types.get(names[i], ""), codes[i]) for i in range(len(names))], codes


def two_sd_bound(vals):
    """mean + 2 SD, the criterion both classifications use."""
    if len(vals) < 2:
        return None, None, None
    m, sd = float(np.mean(vals)), float(np.std(vals, ddof=1))
    return m, sd, m + 2 * sd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default="results/gi5_ptype")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    members, _ = load_group("orf1", "GI")
    if not members:
        raise SystemExit("no GI RdRp reference alignment found")
    print(f"GI.P RdRp references held: {len(members)}")
    by_type = {}
    for n, t, _ in members:
        by_type.setdefault(t, []).append(n)
    for t in sorted(by_type):
        mark = "  <- published GI.5 representative" if REFERENCE in by_type[t] else ""
        print(f"  {t:8s} {len(by_type[t]):2d}  {', '.join(sorted(by_type[t]))}{mark}")
    if not any(REFERENCE in v for v in by_type.values()):
        print(f"\nWARNING: {REFERENCE} (Ehime643) is not in our GI RdRp set at all; "
              f"the published GI.P5 representative is missing.")

    idx = {n: i for i, (n, _, _) in enumerate(members)}
    code = {n: c for n, _, c in members}
    D = {}
    for (n1, _, c1), (n2, _, c2) in itertools.combinations(members, 2):
        v = pdist(c1, c2)
        if v is not None:
            D[(n1, n2)] = D[(n2, n1)] = v

    # ---- 2 x SD bound per P-type
    print("\n2 x SD cluster statistics per P-type (RdRp nucleotide p-distance)")
    print(f"{'P-type':9s} {'n':>2s} {'mean':>7s} {'SD':>7s} {'bound':>7s}  {'max within':>10s}")
    stats = {}
    for t in sorted(by_type):
        ms = by_type[t]
        vals = [D[(x, y)] for x, y in itertools.combinations(ms, 2) if (x, y) in D]
        m, sd, b = two_sd_bound(vals)
        stats[t] = dict(n=len(ms), mean=m, sd=sd, bound=b,
                        max_within=max(vals) if vals else None)
        f = lambda v: f"{v:7.4f}" if isinstance(v, float) else f"{'-':>7s}"
        print(f"{t:9s} {len(ms):2d} {f(m)} {f(sd)} {f(b)}  "
              f"{f(max(vals)) if vals else '        -':>10s}")

    # ---- where does each GI.5-associated sequence sit?
    print("\nmean distance from each sequence to every P-type cluster "
          "(excluding itself), and whether it is inside that cluster's 2 x SD bound")
    rows = []
    for n, t, c in members:
        means = {}
        for ot, ms in by_type.items():
            vals = [D[(n, m)] for m in ms if m != n and (n, m) in D]
            if vals:
                means[ot] = float(np.mean(vals))
        if not means:
            continue
        inside = [ot for ot, v in means.items()
                  if stats[ot]["bound"] is not None and v <= stats[ot]["bound"]]
        nearest = min(means, key=means.get)
        rows.append(dict(accession=n, published_ptype=t, nearest_cluster=nearest,
                         distance_to_nearest=round(means[nearest], 4),
                         distance_to_own=round(means.get(t, float("nan")), 4) if t in means else "",
                         inside_bound_of="; ".join(sorted(inside)) if inside else "none",
                         agrees_with_published=(t in inside) if inside else False,
                         is_published_reference=(n == REFERENCE)))
    rows.sort(key=lambda r: (r["published_ptype"], r["accession"]))
    print(f"{'accession':12s} {'published':9s} {'nearest':9s} {'d(near)':>8s} {'d(own)':>8s}  inside 2xSD bound of")
    for r in rows:
        star = " *" if r["is_published_reference"] else "  "
        print(f"{r['accession']:12s}{star}{r['published_ptype']:9s} {r['nearest_cluster']:9s} "
              f"{r['distance_to_nearest']:8.4f} {str(r['distance_to_own']):>8s}  {r['inside_bound_of']}")

    # ---- can GI.P1 and GI.P5 be told apart at all?
    print("\nbetween-cluster minimum and mean distances among the GI.P types")
    print(f"{'pair':22s} {'min':>8s} {'mean':>8s}  separable at the 0.266 P-type cut-off?")
    pairs_out = []
    for t1, t2 in itertools.combinations(sorted(by_type), 2):
        vals = [D[(x, y)] for x in by_type[t1] for y in by_type[t2] if (x, y) in D]
        if not vals:
            continue
        mn, mu = min(vals), float(np.mean(vals))
        sep = mn >= cg.SAPOVIRUS_CUTOFFS["orf1"][0]
        pairs_out.append(dict(a=t1, b=t2, min=round(mn, 4), mean=round(mu, 4), separable=sep))
        print(f"{t1 + ' vs ' + t2:22s} {mn:8.4f} {mu:8.4f}  {'yes' if sep else 'NO'}")

    with open(out / "gi5_sequences.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    json.dump(dict(stats={k: v for k, v in stats.items()}, pairs=pairs_out,
                   reference_present=any(REFERENCE in v for v in by_type.values())),
              open(out / "gi5_clusters.json", "w"), indent=1, default=float)

    # ---- verdict, stated as evidence rather than as a decision
    print("\n" + "=" * 78)
    print("READING OF THE EVIDENCE")
    print("=" * 78)
    gi5 = [r for r in rows if r["published_ptype"] == "GI.P5"]
    if not gi5:
        print("No sequence in our GI RdRp set carries the published P-type GI.P5, so the")
        print("question cannot be settled from the references we hold: this is a sampling")
        print("gap, not evidence against the published GI.5[P5] assignment.")
    else:
        agree = [r for r in gi5 if r["agrees_with_published"]]
        print(f"GI.P5 sequences held: {len(gi5)}; inside the GI.P5 2 x SD bound: {len(agree)}")
        for r in gi5:
            print(f"  {r['accession']}: nearest cluster {r['nearest_cluster']} at "
                  f"{r['distance_to_nearest']}, inside bound of {r['inside_bound_of']}")
        unsep = [p for p in pairs_out if not p["separable"] and "GI.P5" in (p["a"], p["b"])]
        if unsep:
            print("\nGI.P5 is closer than the 0.266 P-type cut-off to: "
                  + ", ".join(sorted({p['a'] if p['b'] == 'GI.P5' else p['b'] for p in unsep})))
            print("so the cut-off alone cannot separate GI.P5 from those types in this set.")
    print("\nNo new P-type is proposed here. A new type requires either a named published")
    print("type or at least two sequences from distinct locations; a lone unnamed sequence")
    print("would be provisional (GI.PNA[n]), and DQ366345/Ehime643 keeps its published")
    print("GI.P5 name as a named reference regardless of how many relatives we hold.")
    print(f"\nwrote {out}/gi5_sequences.csv and {out}/gi5_clusters.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
