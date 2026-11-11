#!/usr/bin/env python3
"""Wall-clock timings for `calicityper type` on representative inputs.

Reports what a user actually wants to know before starting a run: how long one
genome takes, how long a plate takes, and how that scales with --threads.

Genogroup matters more than library size: the bootstrap cost grows with the
square of the genogroup's reference count, so a GII-heavy set is slower per
sequence than a GI one. Inputs are therefore drawn from a real, GII-dominated
set when one is given, and from the bundled references otherwise.

    python scripts/benchmark_runtime.py --out results/runtime.tsv \\
        [--source /path/to/rivm_all.fasta] [--threads 1 4 16] [--sizes 1 4 96]

Run it on a compute node, not a shared login node, or the numbers are noise.
"""
from __future__ import annotations

import argparse
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def take(source, n, out):
    from Bio import SeqIO
    recs = []
    for r in SeqIO.parse(source, "fasta"):
        recs.append(r)
        if len(recs) >= n:
            break
    SeqIO.write(recs, out, "fasta")
    return len(recs)


def genogroup_mix(tsv):
    import csv as _csv
    with open(tsv) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    mix = {}
    for r in rows:
        g = r.get("genogroup") or "unassigned"
        mix[g] = mix.get(g, 0) + 1
    return ", ".join(f"{k} {v}" for k, v in sorted(mix.items(), key=lambda kv: -kv[1])[:3])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(ROOT / "calicityper/data/database/refseq_2.fasta"),
                    help="FASTA to draw inputs from [bundled reference genomes]")
    ap.add_argument("--sizes", type=int, nargs="+", default=[1, 4, 96])
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 16])
    ap.add_argument("--out", default=str(ROOT / "results/runtime.tsv"))
    ap.add_argument("--workdir", default="/tmp/calicityper_runtime")
    a = ap.parse_args()

    wd = Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"host {platform.node()}  source {Path(a.source).name}")
    for n in a.sizes:
        fa = wd / f"in_{n}.fasta"
        got = take(a.source, n, fa)
        if got < n:
            print(f"  only {got} sequences available, skipping size {n}")
            continue
        for t in a.threads:
            out = wd / f"out_{n}_{t}.tsv"
            cmd = ["calicityper", "type", str(fa), "-o", str(out), "--threads", str(t)]
            start = time.monotonic()
            p = subprocess.run(cmd, capture_output=True, text=True)
            wall = time.monotonic() - start
            if p.returncode != 0:
                print(f"  n={n} threads={t}: FAILED\n{p.stderr[-400:]}")
                continue
            mix = genogroup_mix(out)
            rows.append(dict(sequences=n, threads=t, wall_seconds=round(wall, 1),
                             seconds_per_sequence=round(wall / n, 2), genogroups=mix))
            print(f"  n={n:4d} threads={t:2d}  {wall:8.1f}s  "
                  f"{wall / n:6.2f} s/seq  ({mix})")

    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    cols = ["sequences", "threads", "wall_seconds", "seconds_per_sequence", "genogroups"]
    with open(outp, "w") as fh:
        fh.write("# host=%s cpus=%s source=%s\n" % (platform.node(),
                                                    len(__import__("os").sched_getaffinity(0)),
                                                    Path(a.source).name))
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
