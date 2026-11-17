#!/usr/bin/env python3
"""Accuracy against the simulated amplicon benchmark of Fitzpatrick et al. 2022.

Zenodo 10.5281/zenodo.7271067 (CC BY 4.0) holds ten simulated norovirus
capsid amplicon datasets whose composition is known by construction: each
sequence carries its expected taxonomy in the FASTA header. Unlike a
comparison against another typing tool, which can only measure agreement,
this measures accuracy against a ground truth that does not depend on any
tool's opinion.

Scope, which must be stated wherever these numbers are used:
  * capsid (VP1) amplicons only, so genotype but no P-type;
  * norovirus only, so nothing about sapovirus;
  * 330-344 nt amplicons generated with the GI (GISKF/GISKR) and GII
    (G2SKF/G2SKR) primer pairs, i.e. the operational typing window rather
    than whole genomes.

Header forms vary across the record, so the expected label is taken as the
part of the final taxonomy field before the first underscore:
    ...;Norovirus;GII;GII.2;GII.2_MZ708533.1.001.expected  -> GII.2
    ...;Norovirus;GII;GII.4_MF140665.1.001.expected        -> GII.4
    ...;Norovirus;GII_FJ692501.1.002.expected              -> GII  (genogroup
                                                              only; excluded
                                                              from genotype
                                                              scoring)

    python scripts/benchmark_zenodo_amplicons.py --data-dir <zenodo dir> --out-dir results/benchmark_zenodo
"""
from __future__ import annotations

import argparse
import collections
import csv
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

GENOTYPE_RE = re.compile(r"^G[IVX]+\.\d+[A-Za-z]*$")


GENOGROUP_RE = re.compile(r"^G[IVX]+$")


def expected_label(header):
    """(genogroup or None, genotype or None, source accession) from a header.

    Header depth varies across the record, so rather than trusting a fixed
    field index, take the deepest field that parses as a genogroup or a
    genotype once its "_ACCESSION.sim.expected" suffix is stripped.
    """
    acc = ""
    m = re.search(r"_([A-Z]{1,2}\d+\.\d+)\.", header)
    if m:
        acc = m.group(1)
    genogroup = genotype = None
    for field in header.split(";"):
        label = field.split("_", 1)[0]
        if GENOTYPE_RE.match(label):
            genotype = label
            genogroup = label.split(".")[0]
        elif GENOGROUP_RE.match(label):
            genogroup = label
    return genogroup, genotype, acc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", default="results/benchmark_zenodo")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--typed", help="reuse an existing typed TSV instead of retyping")
    a = ap.parse_args()
    data, out = Path(a.data_dir), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- assemble the query set, keeping the simulation each sequence came from
    truth, recs = {}, []
    for fa in sorted(data.glob("expected_sequences.*.fasta")):
        sim = fa.name.split(".")[1]
        for r in SeqIO.parse(fa, "fasta"):
            gg, gt, acc = expected_label(r.id)
            qid = f"sim{sim}_{len(truth):04d}"
            truth[qid] = dict(sim=sim, genogroup=gg, genotype=gt, accession=acc,
                              header=r.id, length=len(r.seq))
            r.id, r.description = qid, ""
            recs.append(r)
    q = out / "zenodo_amplicons.fasta"
    SeqIO.write(recs, q, "fasta")
    print(f"{len(recs)} amplicons from {len(set(v['sim'] for v in truth.values()))} simulations")
    print(f"  with an expected genotype: {sum(v['genotype'] is not None for v in truth.values())}")

    typed = Path(a.typed) if a.typed else out / "zenodo_typed.tsv"
    if not a.typed:
        print("typing ...", flush=True)
        subprocess.run([sys.executable, "-m", "calicityper.cli", "type", str(q),
                        "-o", str(typed), "--threads", str(a.threads)], check=True)

    df = pd.read_csv(typed, sep="\t", dtype=str).fillna("")
    import validate_vs_rivm as V

    rows = []
    for _, r in df.iterrows():
        t = truth.get(r["query"])
        if t is None:
            continue
        called_gg = V.genogroup_label(r.get("genogroup", "")) if hasattr(V, "genogroup_label") \
            else r.get("genogroup", "")
        rows.append(dict(
            query=r["query"], sim=t["sim"], accession=t["accession"], length=t["length"],
            expected_genogroup=t["genogroup"] or "", called_genogroup=r.get("genogroup", ""),
            expected_genotype=t["genotype"] or "", called_genotype=V.genotype(r.get("orf2_type", "")),
            support=r.get("orf2_support", ""), sites=r.get("orf2_sites", ""),
            accepted=r.get("orf2_accepted", ""), note=r.get("note", "")))
    c = pd.DataFrame(rows)
    c.to_csv(out / "zenodo_comparison.csv", index=False)

    # ---- genogroup
    gg = c[c.expected_genogroup != ""]
    gg_ok = (gg.called_genogroup == gg.expected_genogroup).sum()
    print(f"\ngenogroup: {gg_ok}/{len(gg)} correct ({100*gg_ok/len(gg):.1f}%)")
    bad = gg[gg.called_genogroup != gg.expected_genogroup]
    for (e, g), n in collections.Counter(zip(bad.expected_genogroup, bad.called_genogroup)).most_common(6):
        print(f"    expected {e:4s} -> called {g or '(none)'}: {n}")

    # ---- genotype
    s = c[(c.expected_genotype != "")].copy()
    called = s[~s.called_genotype.isin(["", "unassigned"])]
    ok = (called.called_genotype == called.expected_genotype).sum()
    print(f"\ngenotype: {len(s)} amplicons with a known genotype")
    print(f"  typed      : {len(called)} ({100*len(called)/len(s):.1f}%)")
    print(f"  correct    : {ok} ({100*ok/len(s):.1f}% of all, {100*ok/max(1,len(called)):.1f}% of typed)")
    print(f"  incorrect  : {len(called)-ok}")
    print(f"  not typed  : {len(s)-len(called)}")
    acc = called[called.accepted == "accepted"]
    if len(acc):
        aok = (acc.called_genotype == acc.expected_genotype).sum()
        print(f"  of calls the tool accepts: {aok}/{len(acc)} correct ({100*aok/len(acc):.1f}%)")

    wrong = called[called.called_genotype != called.expected_genotype]
    if len(wrong):
        print("\n  discordant calls (expected -> called: n):")
        for (e, g), n in collections.Counter(zip(wrong.expected_genotype,
                                                 wrong.called_genotype)).most_common(12):
            print(f"    {e:10s} -> {g:10s} {n}")

    s = s.assign(_typed=~s.called_genotype.isin(["", "unassigned"]),
                 _ok=s.called_genotype == s.expected_genotype)
    per = (s.groupby("expected_genotype")
             .agg(n=("query", "size"), typed=("_typed", "sum"), correct=("_ok", "sum"))
             .reset_index())
    per["accuracy"] = per.correct / per.n
    per.sort_values("n", ascending=False).to_csv(out / "zenodo_per_genotype.csv", index=False)
    print(f"\nwrote {out}/zenodo_comparison.csv and {out}/zenodo_per_genotype.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
