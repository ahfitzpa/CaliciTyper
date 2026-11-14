#!/usr/bin/env python3
"""Add sapovirus genogroups to the genogroup-assignment database.

Two things stop a sapovirus query ever reaching the reference sets, and both
live in the BLAST database rather than in the sets themselves:

1. None of the sapovirus subjects appear in refseq_annotation.csv, so
   typing.region_overlap() returns 0 for every one of them and plan_query
   declines the query for "<100 nt overlap" before any typing is attempted.

2. The database holds 12 sapovirus genogroups while the reference sets, built
   from Zhao et al. 2024, cover 34. A GXXIX query is therefore labelled as
   whichever of the 12 it happens to hit, and routed to the wrong group's
   reference set -- a wrong answer rather than a declined one.

Worse, two of the existing labels CONTRADICT the scheme the reference sets use:
MF766258 is labelled GVII but is GXXIV in Zhao's table, and NC_017936 is
labelled "Bat_sapovirus" but is GXIV. Left alone, those two would route queries
to a reference set for a different genogroup and return a confident wrong
genotype.

So the sapovirus half of the database is rebuilt FROM the published table: one
representative genome per genogroup, labelled with the table's genogroup, and
the four prototype subjects that predate the table (the GI, GII, GIV and GXII
reference genomes) kept with their own labels. Region coordinates come from the
same canonical projection used to build the reference sets
(extract_sapovirus_regions.py), so the annotation and the sets cannot disagree
about where RdRp and VP1 are.

Nothing norovirus is touched, but refseq_2.fasta and the BLAST database are
shared, so their checksums change: re-run the manifest build afterwards, and
the Zenodo archive has to be re-deposited before release.

    python scripts/build_sapovirus_genogroup_db.py \\
        --genbank ~/scratch/calicityper_runs/sapovirus_build/sapovirus.gb [--apply]
"""
from __future__ import annotations

import argparse
import collections
import csv
import importlib.util as _u
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

_spec = _u.spec_from_file_location("exsap", ROOT / "scripts" / "extract_sapovirus_regions.py")
ex = _u.module_from_spec(_spec)
_spec.loader.exec_module(ex)

FASTA = ROOT / "calicityper" / "data" / "database" / "refseq_2.fasta"
ANNOT = ROOT / "calicityper" / "data" / "database" / "refseq_annotation.csv"
DB = ROOT / "calicityper" / "data" / "database" / "refseqdb_2"
PREFIX = "Caliciviridae_Sapovirus_"
# Subjects that predate Zhao et al.'s table and are not in it. These are the
# genogroup prototypes (Dresden = GI, C12 = GII, Ehime1596 = GIV, WD1237 =
# GXII), so their existing labels stand; they are kept, not re-derived.
KEEP = {"NC_006269": "GI", "NC_006554": "GII", "DQ366346": "GIV", "KX000385": "GXII"}


def published_genogroups(path=ROOT / "validation" / "sapovirus_vp1_taxonomy.csv"):
    """{genogroup: [accession]} from the published VP1 table."""
    out = collections.defaultdict(list)
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        if r["note"].startswith("UNRESOLVED"):
            continue
        out[r["genotype"].split(".")[0]].append(r["accession"].split(".")[0])
    return out


def region_coordinates(rec, anchors, aln):
    """1-based inclusive genome coordinates of RdRp and VP1, or None.

    The polyprotein bounds come from the canonical projection, so these are the
    same boundaries the reference sets were cut on.
    """
    cds = ex.orf1_cds(rec)
    if cds is None:
        return None
    prot = ex.protein(cds, rec)
    best = None
    for acc, aprot, abounds in anchors:
        mapped, ident = ex.map_bounds(aprot, abounds, prot, aln)
        if best is None or ident > best[1]:
            best = (mapped, ident)
        if ident >= 0.60:
            break
    mapped, ident = best
    if not mapped or "rdrp" not in mapped or "vp1" not in mapped:
        return None
    off = int(cds.location.start) + ex.codon_start(cds) - 1      # 0-based nt of aa 0
    out = {}
    for region in ("rdrp", "vp1"):
        s, e = mapped[region]
        out[region] = (off + s * 3 + 1, off + e * 3)             # 1-based inclusive
    out["identity"] = round(ident, 4)
    out["aa"] = {k: mapped[k] for k in ("rdrp", "vp1")}
    return out


def pick(group, accs, recs, anchors, aln):
    """The longest record of `group` whose two regions both project."""
    cands = sorted((len(recs[a].seq), a) for a in accs if a in recs)
    for _, acc in reversed(cands):
        coords = region_coordinates(recs[acc], anchors, aln)
        if coords:
            return acc, coords
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genbank", required=True)
    ap.add_argument("--apply", action="store_true",
                    help="write the database and annotation (otherwise report only)")
    a = ap.parse_args()

    recs = {r.id.split(".")[0]: r for r in SeqIO.parse(a.genbank, "genbank")}
    print(f"{len(recs)} GenBank records available")
    aln = ex.aligner()
    anchors = ex.collect_anchors(list(recs.values())) if hasattr(ex, "collect_anchors") else None
    anchors = ex.project_canonical(list(recs.values()), aln)
    print(f"{len(anchors)} projection anchors")

    chosen, failed = {}, []
    for group, accs in sorted(published_genogroups().items()):
        acc, coords = pick(group, accs, recs, anchors, aln)
        if acc is None:
            failed.append(group)
            continue
        chosen[group] = (acc, coords)
    # the four prototypes the table does not list
    for acc, group in KEEP.items():
        if acc in recs:
            coords = region_coordinates(recs[acc], anchors, aln)
            if coords:
                chosen.setdefault(group, (acc, coords))
            else:
                print(f"  {group}: prototype {acc} did not project", file=sys.stderr)
        else:
            print(f"  {group}: prototype {acc} is not in {a.genbank}", file=sys.stderr)

    print(f"\n{len(chosen)} sapovirus genogroups with a subject and coordinates")
    print("%-9s %-12s %-7s %-15s %-15s" % ("group", "accession", "ident", "rdrp", "vp1"))
    for group, (acc, c) in sorted(chosen.items()):
        print("%-9s %-12s %-7.3f %-15s %-15s" % (
            group, acc, c["identity"], f"{c['rdrp'][0]}-{c['rdrp'][1]}",
            f"{c['vp1'][0]}-{c['vp1'][1]}"))
    if failed:
        print(f"\nno usable subject for: {', '.join(failed)}")

    old = [r for r in SeqIO.parse(FASTA, "fasta")]
    kept = [r for r in old if not r.id.split("|")[0].startswith(PREFIX)
            and "Sapovirus" not in r.id.split("|")[0]]
    print(f"\nreplacing {len(old) - len(kept)} sapovirus subjects with {len(chosen)}; "
          f"{len(kept)} non-sapovirus subjects untouched")
    if not a.apply:
        print("report only; pass --apply to write")
        return 0

    shutil.copy(FASTA, FASTA.with_suffix(".fasta.bak"))
    shutil.copy(ANNOT, ANNOT.with_suffix(".csv.bak"))
    with open(FASTA, "w") as fh:
        for r in kept:
            fh.write(f">{r.description}\n{str(r.seq).upper()}\n")
        for group, (acc, _) in sorted(chosen.items()):
            rec = recs[acc]
            fh.write(f">{PREFIX}{group}|{rec.id} {rec.description}\n{str(rec.seq).upper()}\n")

    rows = list(csv.DictReader(open(ANNOT)))
    have = {r["id"].split(".")[0] for r in rows}
    added = 0
    for group, (acc, c) in sorted(chosen.items()):
        if acc in have:
            continue
        rows.append(dict(id=recs[acc].id, rdrp_start=c["rdrp"][0], rdrp_end=c["rdrp"][1],
                         vp1_start=c["vp1"][0], vp1_end=c["vp1"][1]))
        added += 1
    with open(ANNOT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "rdrp_start", "rdrp_end", "vp1_start", "vp1_end"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {FASTA} and {ANNOT} (+{added} annotation rows); .bak kept")

    cmd = ["makeblastdb", "-in", str(FASTA), "-dbtype", "nucl", "-out", str(DB),
           "-title", "calicityper refseq_2"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout + p.stderr, file=sys.stderr)
        return 1
    print("rebuilt the BLAST database; refseq_2.fasta and the database checksums have "
          "changed, so rebuild the manifest and re-deposit on Zenodo before release")
    return 0


if __name__ == "__main__":
    sys.exit(main())
