#!/usr/bin/env python3
"""Leave-one-out validation on the bundled reference sets (no internet needed).
Each reference ORF sequence is typed against all OTHER references of its genogroup; the answer is its own label.
Only types with >= 2 references are evaluated (a type with one reference cannot be recovered when it is left out).
Usage: python tests/loo_validation.py [--bootstraps 100] [--threads 2] [--max-refs-per-type N] [--out results/loo.csv]"""
import argparse, sys, collections
from multiprocessing import Pool
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pandas as pd
from Bio import SeqIO
from norotyper import typing as ty, config as cg

A = None
def work(job):
    region, group, acc, truth, seq = job
    r = ty.type_region(seq, group, region, exclude=[acc], nrep=A.bootstraps, min_support=A.min_support,
                       max_refs_per_type=A.max_refs_per_type, backend=A.backend)
    return dict(region=region, group=group, accession=acc, truth=truth, call=r["type"] or "unassigned",
                support=r["support"], best=r.get("best_type", ""), n_sites=r["n_sites"])

def init(a):
    global A; A = a

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstraps", type=int, default=100); ap.add_argument("--min-support", type=float, default=70)
    ap.add_argument("--threads", type=int, default=2); ap.add_argument("--max-refs-per-type", type=int, default=2)
    ap.add_argument("--groups", default="GI,GII"); ap.add_argument("--regions", default="orf1,orf2"); ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="results/loo.csv"); ap.add_argument("--backend", default="mafft", choices=["mafft","clustalw"])
    a = ap.parse_args()
    jobs = []
    for region in a.regions.split(","):
        for g in a.groups.split(","):
            orf = cg.REF_DIR[region] / f"{g}.orf"
            recs = list(SeqIO.parse(orf, "fasta"))
            cnt = collections.Counter(r.id.split("|")[1] for r in recs)
            for r in recs:
                t = r.id.split("|")[1]
                if cnt[t] >= 2:
                    jobs.append((region, g, r.id.split("|")[2], t, str(r.seq)))
    if a.limit: jobs = jobs[::max(1, len(jobs)//a.limit)][:a.limit]
    with Pool(a.threads, initializer=init, initargs=(a,)) as p:
        rows = p.map(work, jobs, chunksize=1)
    df = pd.DataFrame(rows); Path(a.out).parent.mkdir(exist_ok=True, parents=True); df.to_csv(a.out, index=False)
    df["assigned"] = df.call != "unassigned"; df["correct"] = df.call == df.truth
    for region, g in df.groupby("region"):
        print(region, "n=%d assigned=%.1f%% correct(of all)=%.1f%% wrong=%d" % (
            len(g), 100*g.assigned.mean(), 100*g.correct.mean(), int((g.assigned & ~g.correct).sum())))
    print(df[df.assigned & ~df.correct][["region","accession","truth","call","support"]].to_string(index=False))
