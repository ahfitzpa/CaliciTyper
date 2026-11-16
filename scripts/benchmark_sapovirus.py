#!/usr/bin/env python3
"""Sapovirus genotyping benchmark (Zhao et al. 2024 distance cut-offs).

Two independent checks; the third leg, full-length leave-one-out over the
reference strains, is the length-0 rows of scripts/truncation_loo_sapovirus.py
and is summarised here if that table is given with --loo.

1. Published distance ranges. Zhao et al. publish, per genotype, the range of
   pairwise nucleotide p-distances among its strains. Recomputing those from
   our aligned reference sets tests that the sets (extraction, boundaries,
   alignment, the six dual-listed strains' resolution) reproduce the paper,
   which norovirus cannot offer. Also reported per genotype: the smallest
   distance to any OTHER genotype of the genogroup, against the cut-off -- a
   genotype closer than the cut-off to another cannot be told apart by it.

2. RIVM comparison. Every sequence RIVM's tool calls sapovirus is typed with
   `calicityper type --type-sapovirus`, and the VP1 genotype compared with
   RIVM's capsid type (RIVM gives no sapovirus P-type). RIVM uses the older
   nomenclature in places where Zhao et al. split lineages into new
   genogroups, so disagreements are tabulated as (RIVM -> ours) pairs for
   reading, not just counted. Sequences that are themselves reference strains
   are reported separately: they match themselves at distance 0.

Coverage and distance are recorded for every call, so acceptance can be
re-applied with whatever floor the truncation curve gives, without re-running
(--floor, default the current config value).

    python scripts/benchmark_sapovirus.py --rivm-fasta <rivm_all.fasta> --out-dir results/benchmark_sapovirus
"""
from __future__ import annotations

import argparse
import collections
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import config as cg, tree   # noqa: E402

REGIONS = {"orf1": ("RdRp", ROOT / "validation" / "sapovirus_rdrp_taxonomy.csv", "ptype"),
           "orf2": ("VP1", ROOT / "validation" / "sapovirus_vp1_taxonomy.csv", "genotype")}


def published_ranges(path, key):
    out = {}
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        rng = r["pairwise_distance_range"].strip()
        if rng and rng != "NA":
            lo, _, hi = rng.partition("-")
            out[r[key]] = (float(lo), float(hi or lo))
    return out


def pdist(a, b):
    both = (a < 4) & (b < 4)
    n = int(both.sum())
    return (float((a[both] != b[both]).sum()) / n, n) if n >= 50 else (None, n)


def distance_ranges():
    rows = []
    for region, (label, tax, key) in REGIONS.items():
        pub = published_ranges(tax, key)
        gt_cut, _ = cg.SAPOVIRUS_CUTOFFS[region]
        d = cg.SAPOVIRUS_REF_DIR[region]
        for orf in sorted(d.glob("*.orf")):
            if orf.stem == "all":
                continue
            types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
            aln = AlignIO.read(orf.with_suffix(".aln"), "clustal")
            names = [r.id for r in aln]
            codes = tree.encode([str(r.seq).upper() for r in aln])
            within, between = collections.defaultdict(list), collections.defaultdict(list)
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    dist, _ = pdist(codes[i], codes[j])
                    if dist is None:
                        continue
                    ti, tj = types[names[i]], types[names[j]]
                    if ti == tj:
                        within[ti].append(dist)
                    else:
                        between[ti].append(dist)
                        between[tj].append(dist)
            for t in sorted(set(types.values())):
                w = within.get(t, [])
                lo, hi = pub.get(t, (None, None))
                ours_hi = max(w) if w else None
                rows.append(dict(
                    region=label, genogroup=orf.stem, type=t,
                    n=sum(1 for v in types.values() if v == t),
                    ours_min=round(min(w), 4) if w else "", ours_max=round(ours_hi, 4) if w else "",
                    published=f"{lo}-{hi}" if lo is not None else "NA",
                    max_within_published=(ours_hi <= hi + 0.005) if (w and hi is not None) else "",
                    nearest_other_type=round(min(between[t]), 4) if between.get(t) else "",
                    separable_at_cutoff=(min(between[t]) >= gt_cut) if between.get(t) else "",
                    cutoff=gt_cut))
    return pd.DataFrame(rows)


def rivm_compare(rivm_fasta, out_dir, threads, floor):
    rivm = {r["name"]: r for r in csv.DictReader(open(ROOT / "validation" / "results_rivm.csv"))}
    ids = {k for k, r in rivm.items() if "Sapovirus" in r["BLAST"]}
    recs = [r for r in SeqIO.parse(rivm_fasta, "fasta") if r.id in ids]
    q = out_dir / "rivm_sapovirus.fasta"
    SeqIO.write(recs, q, "fasta")
    typed = out_dir / "rivm_sapovirus_typed.tsv"
    print(f"typing {len(recs)} RIVM sapovirus sequences ...", flush=True)
    subprocess.run([sys.executable, "-m", "calicityper.cli", "type", str(q), "-o", str(typed),
                    "--type-sapovirus", "--threads", str(threads)], check=True)
    refs = set()
    for region in REGIONS:
        for orf in cg.SAPOVIRUS_REF_DIR[region].glob("*.orf"):
            refs |= {r.id.split("|")[2].split(".")[0] for r in SeqIO.parse(orf, "fasta")}
    df = pd.read_csv(typed, sep="\t", dtype=str).fillna("")
    out = []
    for _, r in df.iterrows():
        t = rivm[r["query"]]
        out.append(dict(
            query=r["query"], reference_strain=r["query"].split(".")[0] in refs,
            rivm_genogroup=t["BLAST"].replace("Caliciviridae Sapovirus ", ""),
            our_genogroup=r["genogroup"].replace(cg.SAPOVIRUS_PREFIX, ""),
            rivm_vp1=t["capsid type"], our_vp1=r.get("orf2_type", ""),
            vp1_distance=r.get("orf2_distance", ""), vp1_sites=r.get("orf2_sites", ""),
            vp1_coverage=r.get("orf2_coverage", ""),
            vp1_accepted=r.get("orf2_accepted", ""), our_rdrp=r.get("orf1_type", ""),
            rdrp_distance=r.get("orf1_distance", ""), rdrp_accepted=r.get("orf1_accepted", ""),
            note=r.get("note", "")))
    cmp_ = pd.DataFrame(out)
    # acceptance at the floor asked for, re-derived from the recorded columns:
    # within the genotype cut-off AND covering >= floor of the region
    gt_cut = cg.SAPOVIRUS_CUTOFFS["orf2"][0]
    d = pd.to_numeric(cmp_.vp1_distance, errors="coerce")
    cv = pd.to_numeric(cmp_.vp1_coverage, errors="coerce")
    cmp_["vp1_accepted_at_floor"] = (d < gt_cut) & (cv >= floor)
    cmp_.to_csv(out_dir / "rivm_comparison.csv", index=False)
    return cmp_


def summarise_rivm(c):
    print("\n=== RIVM comparison (VP1) ===")
    print(f"sequences: {len(c)}; genogroup agrees with RIVM: "
          f"{(c.rivm_genogroup == c.our_genogroup).sum()}")
    for name, sub in (("all", c), ("not reference strains", c[~c.reference_strain])):
        s = sub[sub.rivm_vp1 != ""]
        called = s[~s.our_vp1.isin(["", "unassigned", "untypeable"])]
        agree = called[called.our_vp1 == called.rivm_vp1]
        print(f"\n[{name}] RIVM VP1-typed {len(s)}: called {len(called)}, agree {len(agree)}, "
              f"differ {len(called) - len(agree)}, untypeable "
              f"{(s.our_vp1 == 'untypeable').sum()}, no call {s.our_vp1.isin(['', 'unassigned']).sum()}")
        acc = called[called.vp1_accepted_at_floor]
        print(f"  accepted at the coverage floor: {len(acc)}, of which agree "
              f"{(acc.our_vp1 == acc.rivm_vp1).sum()}")
        pairs = collections.Counter(zip(called.rivm_vp1, called.our_vp1))
        diffs = [(k, n) for k, n in pairs.most_common() if k[0] != k[1]]
        if diffs:
            print("  differences (RIVM -> ours: n):")
            for (a, b), n in diffs:
                print(f"    {a:10s} -> {b:12s} {n}")
    gg = collections.Counter((a, b) for a, b in zip(c.rivm_genogroup, c.our_genogroup) if a != b)
    if gg:
        print("\ngenogroup differences (RIVM -> ours: n):")
        for (a, b), n in gg.most_common(20):
            print(f"    {a:16s} -> {b:12s} {n}")


def summarise_loo(path):
    d = pd.read_csv(path)
    f = d[d.length == 0]
    print("\n=== full-length leave-one-out (from the truncation run) ===")
    for region, s in f.groupby("region"):
        k, sg = s[~s.singleton], s[s.singleton]
        print(f"{region}: {len(k)} with a same-type partner -> "
              f"{(k.outcome == 'correct').mean():.1%} correct, {(k.outcome == 'wrong').mean():.1%} wrong, "
              f"{(k.outcome == 'untypeable').mean():.1%} untypeable; {len(sg)} singletons -> "
              f"{(sg.outcome == 'correct').mean():.1%} correctly untypeable")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rivm-fasta", required=True)
    ap.add_argument("--out-dir", default="results/benchmark_sapovirus")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--floor", type=float, default=cg.SAPOVIRUS_MIN_REGION_COVERAGE)
    ap.add_argument("--loo", help="truncation_loo_sapovirus.csv, to summarise its full-length rows")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("=== published distance ranges ===", flush=True)
    dr = distance_ranges()
    dr.to_csv(out / "distance_ranges.csv", index=False)
    for region, s in dr.groupby("region"):
        cmp_ = s[s.max_within_published != ""]
        sep = s[s.separable_at_cutoff != ""]
        print(f"{region}: {len(s)} types; max within-type distance inside the published range for "
              f"{int(cmp_.max_within_published.astype(bool).sum())}/{len(cmp_)}; separated from every "
              f"other type of its genogroup by >= the cut-off for "
              f"{int(sep.separable_at_cutoff.astype(bool).sum())}/{len(sep)}")
        bad = cmp_[~cmp_.max_within_published.astype(bool)]
        if len(bad):
            print("  outside the published range:")
            print(bad[["genogroup", "type", "n", "ours_max", "published"]].to_string(index=False))
        close = sep[~sep.separable_at_cutoff.astype(bool)]
        if len(close):
            print("  closer than the cut-off to another type:")
            print(close[["genogroup", "type", "nearest_other_type", "cutoff"]].to_string(index=False))

    print(f"\ncoverage floor for acceptance: {a.floor}")
    c = rivm_compare(a.rivm_fasta, out, a.threads, a.floor)
    summarise_rivm(c)
    if a.loo and Path(a.loo).exists():
        summarise_loo(a.loo)
    print(f"\nwrote {out}/distance_ranges.csv, {out}/rivm_comparison.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
