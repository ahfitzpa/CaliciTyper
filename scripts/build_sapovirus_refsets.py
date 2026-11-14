#!/usr/bin/env python3
"""Build and align the sapovirus reference sets, in the norovirus layout.

Inputs are the extracted regions (scripts/extract_sapovirus_regions.py) and the
verified taxonomy tables (scripts/build_sapovirus_taxonomy.py). Output per
region, under reference_sequences_sapovirus_{VP1,RdRp}:

    {group}.fasta        nucleotides, one file per genogroup/P-group
    {group}.orf          the same, headers "group|type|accession"
    {group}.aln          MAFFT alignment in CLUSTAL format
    all.orf, all.aln     every group together

Only records carrying the COMPLETE canonical window are included. The paper
defines its cutoffs on complete sequences, so a reference truncated to 290 of
512 aa would make distances incomparable with the <0.266 P-type rule. Records
that fail the window, carry excessive ambiguity, or whose group has too few
members to support a distance range are excluded and listed in excluded.tsv --
a group with no complete reference must report "no complete reference" rather
than be typed to its nearest neighbour.

    python scripts/build_sapovirus_refsets.py --regions <dir> [--recovered <fasta>]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
MAX_AMBIGUOUS = 0.02          # fraction of N/ambiguity tolerated in a reference
MIN_MEMBERS = 1               # a group needs at least this many complete records


def load_taxonomy(path, key):
    out = {}
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        if r["note"].startswith("UNRESOLVED"):
            continue
        out[r["accession"].split(".")[0]] = r[key]
    return out


def ambiguous_fraction(seq):
    s = str(seq).upper()
    return sum(1 for c in s if c not in "ACGTU") / max(1, len(s))


def mafft_clustal(in_fasta, out_aln, threads=4):
    cmd = ["mafft", "--quiet", "--auto", "--clustalout",
           "--thread", str(threads), str(in_fasta)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"mafft failed on {in_fasta}:\n{res.stderr[-500:]}")
    Path(out_aln).write_text(res.stdout)


def build(region, fasta, taxonomy, out_dir, recovered=None, threads=4):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seqs = {r.id.split(".")[0]: r for r in SeqIO.parse(fasta, "fasta")}
    if recovered and Path(recovered).exists():
        for r in SeqIO.parse(recovered, "fasta"):
            acc = r.id.split(".")[0]
            seqs.setdefault(acc, r)
            # recovered records carry their P-type in the description
            if len(r.description.split()) > 1:
                taxonomy.setdefault(acc, r.description.split()[1])

    groups, excluded = defaultdict(list), []
    for acc, rec in sorted(seqs.items()):
        label = taxonomy.get(acc)
        if not label:
            excluded.append(dict(region=region, accession=acc, label="",
                                 reason="not in the published table"))
            continue
        amb = ambiguous_fraction(rec.seq)
        if amb > MAX_AMBIGUOUS:
            excluded.append(dict(region=region, accession=acc, label=label,
                                 reason=f"{amb:.1%} ambiguous bases, above "
                                        f"{MAX_AMBIGUOUS:.0%}"))
            continue
        groups[label.split(".")[0]].append((label, acc, str(rec.seq)))

    written = {}
    for group, members in sorted(groups.items()):
        if len(members) < MIN_MEMBERS:
            for label, acc, _ in members:
                excluded.append(dict(region=region, accession=acc, label=label,
                                     reason=f"group {group} has only {len(members)} "
                                            f"complete record(s)"))
            continue
        fa = out_dir / f"{group}.fasta"
        orf = out_dir / f"{group}.orf"
        with open(fa, "w") as fh_fa, open(orf, "w") as fh_orf:
            for label, acc, seq in sorted(members):
                fh_fa.write(f">{acc}\n{seq}\n")
                fh_orf.write(f">{group}|{label}|{acc}\n{seq}\n")
        if len(members) > 1:
            mafft_clustal(orf, out_dir / f"{group}.aln", threads)
        written[group] = members
        print(f"  {region} {group:10} {len(members):3d} records, "
              f"{len({m[0] for m in members})} types"
              + ("" if len(members) > 1 else "  (single record: no alignment)"))

    allp = out_dir / "all.orf"
    with open(allp, "w") as fh:
        for group, members in sorted(written.items()):
            for label, acc, seq in sorted(members):
                fh.write(f">{group}|{label}|{acc}\n{seq}\n")
    if sum(len(m) for m in written.values()) > 1:
        mafft_clustal(allp, out_dir / "all.aln", threads)
    return written, excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", required=True, help="output dir of extract_sapovirus_regions.py")
    ap.add_argument("--recovered", help="recovered_rdrp.fasta from recover_sapovirus_rdrp.py")
    ap.add_argument("--out-dir", default=str(ROOT / "calicityper" / "data"))
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    regions = Path(a.regions)
    out = Path(a.out_dir)
    excluded_all = []
    summary = {}
    for region, src, tax_csv, key, dest in (
        ("VP1", regions / "sapovirus_vp1.fasta",
         ROOT / "validation/sapovirus_vp1_taxonomy.csv", "genotype",
         out / "reference_sequences_sapovirus_VP1"),
        ("RdRp", regions / "sapovirus_rdrp.fasta",
         ROOT / "validation/sapovirus_rdrp_taxonomy.csv", "ptype",
         out / "reference_sequences_sapovirus_RdRp"),
    ):
        if not src.exists():
            print(f"missing {src}, skipping {region}")
            continue
        print(f"{region}:")
        written, excluded = build(region, src, load_taxonomy(tax_csv, key), dest,
                                  a.recovered if region == "RdRp" else None, a.threads)
        excluded_all += excluded
        types = {m[0] for members in written.values() for m in members}
        summary[region] = dict(groups=len(written),
                               records=sum(len(m) for m in written.values()),
                               types=len(types))

    with open(out / "sapovirus_excluded.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["region", "accession", "label", "reason"])
        w.writeheader()
        w.writerows(excluded_all)
    print("\nsummary:")
    for region, s in summary.items():
        print(f"  {region}: {s['records']} records, {s['types']} types, "
              f"{s['groups']} groups")
    print(f"  excluded: {len(excluded_all)} (see sapovirus_excluded.tsv)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
