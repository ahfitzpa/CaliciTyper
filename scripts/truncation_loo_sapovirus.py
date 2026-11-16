#!/usr/bin/env python3
"""Sapovirus truncation leave-one-out: how much of a region a call needs.

The sapovirus calls are nucleotide p-distance cut-offs (Zhao et al. 2024:
VP1 <0.161 genotype, RdRp <0.266 P-type) derived on COMPLETE regions. Whether
they still hold on a fragment is an empirical question, and until it is
answered the typing path requires SAPOVIRUS_MIN_REGION_COVERAGE (provisionally
0.9). This harness measures it, so the floor is read off a curve rather than
asserted.

For every reference in the sapovirus sets, its own region sequence is cut to a
series of lengths at three positions (5' end, centre, 3' end) -- sapovirus has
no operational typing window of the RIVM kind to anchor amplicons on -- and
each fragment is typed against all OTHER references of its genogroup, its own
accession removed from the alignment and the tree. The call is made exactly as
calicityper type makes it (typing._apply_sapovirus_result), and scored against
the record's published type:

    correct    called the published type
    wrong      called a different type            (a precision failure)
    untypeable no type within the cut-off         (a recall loss)

A reference whose type has no other member is a SINGLETON: once it is left
out, its true type is absent, so the right answer is "untypeable" and any
label is wrong. Singletons therefore measure the cut-off's specificity -- how
often a novel genotype is given a known label -- and are scored that way.

Recorded WITHOUT a coverage threshold; analyse() sweeps the floor over the
table, so one pass answers "what should the floor be?".

    python scripts/truncation_loo_sapovirus.py --out results/truncation_loo_sapovirus.csv --threads 16
    python scripts/truncation_loo_sapovirus.py --analyse-only results/truncation_loo_sapovirus.csv
"""
from __future__ import annotations

import argparse
import collections
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import align, config as cg, typing as ty   # noqa: E402

LENGTHS = (150, 200, 300, 400, 500, 600, 800, 1000, 1200, 0)   # nt; 0 = whole region
POSITIONS = ("5p", "mid", "3p")
REGIONS = {"orf1": "RdRp", "orf2": "VP1"}
SWEEP = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
A = None


def references(region):
    """[(group label, accession, type, ungapped region seq)] for one region."""
    d = cg.SAPOVIRUS_REF_DIR[region]
    out = []
    for orf in sorted(d.glob("*.orf")):
        if orf.stem == "all":
            continue
        group = cg.SAPOVIRUS_PREFIX + orf.stem
        for r in SeqIO.parse(orf, "fasta"):
            _, typ, acc = r.id.split("|")
            out.append((group, acc, typ, str(r.seq).upper().replace("-", "")))
    return out


def fragments(seq, lengths):
    """[(length, position, fragment)]; the whole region once, as length 0."""
    out = []
    for n in lengths:
        if n == 0 or n >= len(seq):
            if n == 0:
                out.append((0, "full", seq))
            continue
        starts = {"5p": 0, "mid": (len(seq) - n) // 2, "3p": len(seq) - n}
        out += [(n, pos, seq[s:s + n]) for pos, s in starts.items()]
    return out


def work(job):
    """One reference: align all its fragments in one MAFFT call, type each."""
    region, group, acc, truth, singleton, frags = job
    ref, alns = align.align_queries(region, group, [(f"{n}_{p}", s) for n, p, s in frags],
                                    exclude=[acc], threads=1)
    rows = []
    for (n, pos, frag), qa in zip(frags, alns):
        res = ty.type_aligned(ref, qa, exclude=[acc], nrep=A.bootstraps, seed=1,
                              min_support=70.0, max_refs_per_type=A.max_refs_per_type)
        row = dict(query=acc, genogroup=group, note="")
        ty.apply_result(row, region, res, min_support=70.0)
        call = row.get(f"{region}_type", "")
        span = int(res.get("region_sites") or 0)
        sites = int(res.get("n_sites") or 0)
        if singleton:
            outcome = "correct" if call in ("untypeable", "unassigned") else "wrong"
        elif call == truth:
            outcome = "correct"
        elif call in ("untypeable", "unassigned"):
            outcome = "untypeable"
        else:
            outcome = "wrong"
        rows.append(dict(region=REGIONS[region], genogroup=group[len(cg.SAPOVIRUS_PREFIX):],
                         accession=acc, truth=truth, singleton=singleton, length=n,
                         position=pos, frag_nt=len(frag), sites=sites, region_sites=span,
                         coverage=round(sites / span, 4) if span else 0.0,
                         distance=row.get(f"{region}_distance", ""),
                         nearest_type=row.get(f"{region}_nearest_type", ""),
                         call=call, bootstrap_best=res.get("best_type", ""),
                         support=row.get(f"{region}_support", ""), outcome=outcome))
    return rows


def init(a):
    global A
    A = a


def build_jobs(regions, lengths, limit=0):
    jobs = []
    for region in regions:
        refs = references(region)
        per_type = collections.Counter((g, t) for g, _, t, _ in refs)
        per_group = collections.Counter(g for g, _, _, _ in refs)
        alone = [a for g, a, _, _ in refs if per_group[g] < 2]
        if alone:
            print(f"{REGIONS[region]}: skipping {len(alone)} reference(s) alone in their genogroup "
                  f"(nothing left to type against once left out): {', '.join(alone)}")
        for group, acc, typ, seq in refs:
            if per_group[group] < 2:
                continue
            jobs.append((region, group, acc, typ, per_type[(group, typ)] == 1,
                         fragments(seq, lengths)))
    return jobs[:limit] if limit else jobs


def analyse(df):
    """Precision/recall by coverage floor, per region; singletons separately."""
    out = []
    for region, d in df.groupby("region"):
        known, single = d[~d.singleton], d[d.singleton]
        print(f"\n=== {region}: {known.accession.nunique()} references with a same-type partner, "
              f"{single.accession.nunique()} singletons ===")
        full = known[known.length == 0]
        print(f"full-length LOO: {(full.outcome == 'correct').mean():.1%} correct, "
              f"{(full.outcome == 'wrong').mean():.1%} wrong, "
              f"{(full.outcome == 'untypeable').mean():.1%} untypeable (n={len(full)})")
        sf = single[single.length == 0]
        if len(sf):
            print(f"singletons, full length: {(sf.outcome == 'correct').mean():.1%} correctly "
                  f"untypeable, {(sf.outcome == 'wrong').mean():.1%} given a wrong known label "
                  f"(n={len(sf)})")
        print("floor  accepted  precision  recall   singleton-false-label")
        for floor in SWEEP:
            acc = known[known.coverage >= floor]
            called = acc[acc.outcome != "untypeable"]
            prec = (called.outcome == "correct").mean() if len(called) else float("nan")
            rec = (acc.outcome == "correct").sum() / len(known) if len(known) else float("nan")
            sacc = single[single.coverage >= floor]
            sfalse = (sacc.outcome == "wrong").mean() if len(sacc) else float("nan")
            print(f"{floor:5.2f}  {len(acc):8d}  {prec:9.3f}  {rec:6.3f}   {sfalse:8.3f}")
            out.append(dict(region=region, floor=floor, accepted=len(acc), precision=prec,
                            recall=rec, singleton_false_label=sfalse))
        print("by fragment length (all positions):")
        t = known.groupby("length").outcome.value_counts(normalize=True).unstack(fill_value=0)
        t["mean_cov"] = known.groupby("length").coverage.mean()
        print(t.round(3).to_string())
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="results/truncation_loo_sapovirus.csv")
    ap.add_argument("--analyse-only", help="skip typing, analyse an existing CSV")
    ap.add_argument("--bootstraps", type=int, default=20,
                    help="bootstrap replicates; reported only, the call is by distance [20]")
    ap.add_argument("--max-refs-per-type", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--regions", default="orf1,orf2")
    ap.add_argument("--lengths", default=",".join(str(x) for x in LENGTHS))
    ap.add_argument("--limit", type=int, default=0, help="first N references only, for a smoke test")
    a = ap.parse_args()

    if a.analyse_only:
        df = pd.read_csv(a.analyse_only)
    else:
        cg.TYPE_SAPOVIRUS = True
        jobs = build_jobs(a.regions.split(","), [int(x) for x in a.lengths.split(",")], a.limit)
        n = sum(len(j[-1]) for j in jobs)
        print(f"{len(jobs)} references, {n} fragments ({a.bootstraps} bootstraps each) "
              f"on {a.threads} threads", flush=True)
        rows = []
        with Pool(a.threads, initializer=init, initargs=(a,)) as pool:
            for k, part in enumerate(pool.imap_unordered(work, jobs, chunksize=1), 1):
                rows += part
                if k % 20 == 0:
                    print(f"  {k}/{len(jobs)} references done", flush=True)
        df = pd.DataFrame(rows)
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(a.out, index=False)
        print(f"wrote {a.out} ({len(df)} rows)")
    summary = analyse(df)
    sp = Path(a.out if not a.analyse_only else a.analyse_only).with_suffix(".floors.csv")
    summary.to_csv(sp, index=False)
    print(f"\nwrote {sp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
