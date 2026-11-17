#!/usr/bin/env python3
"""Temporal hold-out evaluation of the putative-novel screen.

The existing leave-one-out design removes a whole genotype and asks whether
the tool flags its members. That is the only way to manufacture "absent from
the reference set" out of known data, and the bounds are recomputed per
scenario so nothing leaks. But a removed genotype evolved inside the surviving
radiation and sits on a shorter branch than a genuinely novel lineage would,
so the design is harder than reality and its sensitivities are lower bounds.

A temporal hold-out tests the same question against real novelty. Reference
strains are split by collection date: everything collected before a cut-off
forms the reference set and the 2 x SD bounds, and everything collected after
it is screened as if newly sequenced. Genotypes that only appear after the
cut-off are genuinely novel with respect to that reference set, and their
novelty was established by the field rather than by us.

Scoring, for each held-out sequence:
    its genotype present before the cut-off   -> expect `known`
    its genotype absent before the cut-off    -> expect a novelty flag
so false positives are over-flagging an established type and false negatives
are missing a type the classification had not yet designated.

Sequences without a collection date are reported separately as "date unknown"
and excluded from the split, never silently assigned to either side.

    python scripts/temporal_holdout_novelty.py --cutoff 2016-01-01 --out-dir results/temporal_holdout
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import putative_novel as nv   # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MONTHS = {m: i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}


def parse_collection_date(text):
    """GenBank /collection_date to a sortable (year, month) or None.

    Accepts the forms GenBank actually uses: 2014, Mar-2014, 12-Mar-2014,
    2014-03, 2014-03-12, and ranges, from which the earliest is taken.
    """
    if not text:
        return None
    t = text.strip().strip('"').split("/")[0].strip()      # a range: take the start
    m = re.match(r"^(\d{4})-(\d{2})", t)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.match(r"^(\d{4})$", t)
    if m:
        return int(m.group(1)), 1
    m = re.match(r"^(?:\d{1,2}-)?([A-Za-z]{3})-(\d{4})$", t)
    if m:
        return int(m.group(2)), MONTHS.get(m.group(1).title(), 1)
    m = re.search(r"(\d{4})", t)
    return (int(m.group(1)), 1) if m else None


def fetch_dates(accessions, cache, batch=100):
    """{accession: (year, month) or None}, from GenBank /collection_date."""
    known = json.load(open(cache)) if Path(cache).exists() else {}
    todo = [a for a in accessions if a not in known]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        q = urllib.parse.urlencode(dict(db="nuccore", id=",".join(chunk),
                                        rettype="gb", retmode="text"))
        req = urllib.request.Request(f"{EUTILS}?{q}",
                                     headers={"User-Agent": "calicityper-benchmark"})
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                text = r.read().decode(errors="replace")
        except Exception as e:
            print(f"  batch {i // batch}: {e}", flush=True)
            time.sleep(10)
            continue
        acc = None
        for line in text.splitlines():
            if line.startswith("VERSION"):
                acc = line.split()[1].split(".")[0]
                known.setdefault(acc, None)
            elif "/collection_date=" in line and acc:
                d = parse_collection_date(line.split("=", 1)[1])
                if d:
                    known[acc] = list(d)
        for a in chunk:
            known.setdefault(a, None)
        print(f"  dates for {min(i + batch, len(todo))}/{len(todo)}", flush=True)
        time.sleep(1)
    json.dump(known, open(cache, "w"))
    return {a: (tuple(v) if v else None) for a, v in known.items()}


def bounds_from(blocks, keep_idx, region):
    """2 x SD bound table built only from the sequences in keep_idx."""
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
    ap.add_argument("--cutoff", default="2016-01-01", help="YYYY-MM-DD collection-date cut-off")
    ap.add_argument("--out-dir", default="results/temporal_holdout")
    ap.add_argument("--measure", default="ml", choices=["ml", "nj", "p"])
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cy, cm = int(a.cutoff[:4]), int(a.cutoff[5:7])
    cut = (cy, cm)

    aa = dict(nv.vp1_aa_alignment()[0])
    _vp1 = nv.vp1_distance_matrix()
    blocks = {"orf2": [(_vp1[0], nv.as_patristic(_vp1[1], a.measure,
                                                 records=[(m[0], aa[m[0]]) for m in _vp1[0]],
                                                 alphabet="aa"))],
              "orf1": [(meta, nv.as_patristic(D, a.measure,
                                              records=nv.rdrp_alignment(meta[0][1]), alphabet="nt"))
                       for meta, D in nv.rdrp_distance_blocks()]}

    accs = sorted({m[0].split(".")[0] for bl in blocks.values() for meta, _ in bl for m in meta})
    print(f"reference strains: {len(accs)}; fetching collection dates from GenBank")
    dates = fetch_dates(accs, out / "collection_dates.json")
    dated = {a_: d for a_, d in dates.items() if d}
    print(f"  with a collection date: {len(dated)}; date unknown: {len(accs) - len(dated)}")
    if not dated:
        raise SystemExit("no collection dates retrieved; cannot split")
    yrs = Counter(d[0] for d in dated.values())
    print("  collection years:", ", ".join(f"{y}:{n}" for y, n in sorted(yrs.items())[:14]), "...")

    rows, tally = [], defaultdict(Counter)
    for region, bl in blocks.items():
        before, after, unknown = set(), set(), set()
        for b, (meta, _) in enumerate(bl):
            for i, (acc, g, t) in enumerate(meta):
                d = dates.get(acc.split(".")[0])
                (before if (d and d < cut) else after if d else unknown).add((b, i))
        print(f"\n[{region}] before {a.cutoff}: {len(before)}, after: {len(after)}, "
              f"date unknown: {len(unknown)}")
        if len(before) < 10 or not after:
            print("  too few on one side of the cut-off to evaluate")
            continue
        table = bounds_from(bl, before, region)
        known_types = {meta[i][2] for b, (meta, _) in enumerate(bl) for i in range(len(meta))
                       if (b, i) in before}
        for b, (meta, D) in enumerate(bl):
            for i, (acc, g, t) in enumerate(meta):
                if (b, i) not in after:
                    continue
                keep = [j for j in range(len(meta)) if (b, j) in before]
                if not keep:
                    continue
                d = D[i, keep]
                flag = nv.classify(region,
                                   nv._averages(d, [meta[j][2] for j in keep]),
                                   nv._averages(d, [meta[j][1] for j in keep]), table)[0]
                novel_truth = t not in known_types
                flagged = flag in nv.FLAGGED
                outcome = ("true positive" if novel_truth and flagged else
                           "false negative" if novel_truth else
                           "false positive" if flagged else "true negative")
                tally[region][outcome] += 1
                rows.append(dict(region=region, accession=acc, genogroup=g, genotype=t,
                                 collection=".".join(str(x) for x in dates[acc.split(".")[0]]),
                                 genotype_in_reference_set=not novel_truth,
                                 flag=flag, outcome=outcome))

    print("\n" + "=" * 78)
    print(f"TEMPORAL HOLD-OUT, cut-off {a.cutoff}, distance measure {a.measure}")
    print("=" * 78)
    for region, c in tally.items():
        tp, fn = c["true positive"], c["false negative"]
        fp, tn = c["false positive"], c["true negative"]
        n = sum(c.values())
        sens = tp / (tp + fn) if (tp + fn) else float("nan")
        spec = tn / (tn + fp) if (tn + fp) else float("nan")
        print(f"\n[{region}] {n} held-out sequences")
        print(f"  genotype already in the reference set : {tn + fp}  "
              f"(correctly known {tn}, falsely flagged novel {fp})")
        print(f"  genotype designated only afterwards   : {tp + fn}  "
              f"(correctly flagged {tp}, missed {fn})")
        print(f"  novelty sensitivity {sens:.3f}  specificity {spec:.3f}")
        if tp + fn == 0:
            print("  NOTE: no genotype in the hold-out is absent from the reference set, so")
            print("  sensitivity is undefined at this cut-off; try an earlier one.")
    with open(out / "temporal_holdout.csv", "w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print(f"\nwrote {out}/temporal_holdout.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
