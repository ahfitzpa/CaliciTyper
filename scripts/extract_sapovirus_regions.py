#!/usr/bin/env python3
"""Extract the RdRp (NS7) and VP1 regions from sapovirus GenBank records.

Sapovirus encodes the nonstructural proteins AND the VP1 capsid in a single
ORF1 polyprotein, so neither region is a CDS of its own and the norovirus
"take ORF1, take ORF2" approach does not apply. Only 18 of the 196 fetched
records annotate the boundaries at all (mat_peptide), and most of those only
annotate VP1.

Method: take the records that DO annotate both regions as anchors, align each
record's ORF1 polyprotein to an anchor's polyprotein (amino acids, global with
free end gaps), and map the anchor's boundaries through that alignment. The
anchor whose alignment identity is highest is used for each record, so a
divergent animal sapovirus is mapped from the closest available annotation
rather than from a distant human one.

Because 18 records carry their own annotation, the mapping error is measurable:
--check reports, for those records, the difference between mapped and annotated
boundaries. Read that number before trusting the rest.

    python scripts/extract_sapovirus_regions.py \\
        --genbank ~/scratch/calicityper_runs/sapovirus_build/sapovirus.gb \\
        --out-dir ~/scratch/calicityper_runs/sapovirus_build/regions
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper.putative_novel import sanitise_aa          # noqa: E402

# Products that identify each region in a mat_peptide annotation. NS6-7 and
# Pro-Pol are the protease-polymerase: the polymerase is its C-terminal half,
# so a boundary taken from them starts upstream of RdRp proper and is only used
# when nothing better exists.
RDRP_WORDS = ("rna dependent rna poly", "rna-dependent rna poly", "rdrp",
              "ns7", "ns6-7", "pro-pol", "polymerase")
VP1_WORDS = ("vp1", "capsid protein", "major structural")
# Plausibility, from the literature and from the records themselves.
VP1_AA_RANGE = (480, 650)
RDRP_AA_RANGE = (250, 700)


def aligner():
    a = PairwiseAligner(mode="global")
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -11.0, -1.0
    a.end_gap_score = 0.0            # a partial record must not be penalised
    return a


def orf1_cds(rec):
    """The ORF1 polyprotein CDS: the longest CDS, which in sapovirus carries
    NS1-NS7 and VP1. ORF2 (VP2) is much shorter."""
    cdss = [f for f in rec.features if f.type == "CDS"]
    if not cdss:
        return None
    return max(cdss, key=lambda f: len(f))


def mat_peptides(rec, words):
    out = []
    for f in rec.features:
        if f.type != "mat_peptide":
            continue
        prod = " ".join(f.qualifiers.get("product", [])).lower()
        if any(w in prod for w in words):
            out.append(f)
    return out


def codon_start(cds):
    """1-based reading-frame offset. 58 of the 403 CDS here are not 1, and
    ignoring it translates the polyprotein out of frame -- the anchor used
    first gave 126 internal stop codons and the mapping collapsed."""
    try:
        return int(cds.qualifiers.get("codon_start", ["1"])[0])
    except (ValueError, TypeError):
        return 1


def aa_bounds(cds, feat):
    """mat_peptide -> (aa start, aa end) within the CDS translation, 0-based."""
    c0, c1 = int(cds.location.start), int(cds.location.end)
    f0, f1 = int(feat.location.start), int(feat.location.end)
    if not (c0 <= f0 < f1 <= c1):
        return None
    off = c0 + codon_start(cds) - 1
    return (f0 - off) // 3, (f1 - off) // 3


def protein(cds, rec):
    """The submitter's own /translation when present -- every CDS in this set
    has one, and it already accounts for codon_start, partial ends and any
    annotated exception. Translating by hand is the fallback."""
    tr = cds.qualifiers.get("translation")
    if tr and tr[0]:
        return tr[0].rstrip("*")
    nt = cds.extract(rec.seq)[codon_start(cds) - 1:]
    nt = nt[: len(nt) // 3 * 3]
    return str(nt.translate()).rstrip("*")


def collect_anchors(recs):
    """Records annotating VP1 and, ideally, RdRp: (id, protein, bounds dict)."""
    anchors = []
    for r in recs:
        cds = orf1_cds(r)
        if cds is None:
            continue
        vp1 = mat_peptides(r, VP1_WORDS)
        rdrp = mat_peptides(r, RDRP_WORDS)
        if not vp1:
            continue
        b = {}
        v = aa_bounds(cds, vp1[0])
        if v:
            b["vp1"] = v
        if rdrp:
            d = aa_bounds(cds, rdrp[0])
            if d:
                b["rdrp"] = d
        if "vp1" not in b:
            continue
        anchors.append((r.id, protein(cds, r), b))
    # prefer anchors that annotate both regions
    anchors.sort(key=lambda t: ("rdrp" not in t[2], t[0]))
    return anchors


def map_bounds(anchor_prot, anchor_bounds, query_prot, aln):
    """Map each anchor region onto the query through a pairwise alignment.
    Returns ({region: (qstart, qend)}, identity)."""
    a = sanitise_aa(anchor_prot)
    q = sanitise_aa(query_prot)
    al = aln.align(a, q)[0]
    pairs = {}
    n = same = 0
    for (a0, a1), (q0, q1) in zip(*al.aligned):
        for k in range(a1 - a0):
            pairs[a0 + k] = q0 + k
            n += 1
            same += a[a0 + k] == q[q0 + k]
    identity = same / n if n else 0.0

    out = {}
    for region, (s, e) in anchor_bounds.items():
        # nearest mapped anchor column at or after the start, and at or before the end
        starts = [pairs[i] for i in range(s, e) if i in pairs]
        if len(starts) < 50:
            continue
        out[region] = (min(starts), max(starts) + 1)
    return out, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genbank", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-identity", type=float, default=0.35,
                    help="reject a mapping below this anchor/query aa identity")
    a = ap.parse_args()

    recs = list(SeqIO.parse(a.genbank, "genbank"))
    anchors = collect_anchors(recs)
    both = [x for x in anchors if "rdrp" in x[2]]
    print(f"{len(recs)} records; {len(anchors)} annotate VP1, {len(both)} annotate RdRp too")
    if not both:
        sys.exit("no record annotates both regions: cannot anchor the mapping")

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    aln = aligner()
    rows, seqs = [], {"rdrp": [], "vp1": []}
    checks = []

    for r in recs:
        cds = orf1_cds(r)
        if cds is None:
            rows.append(dict(accession=r.id, status="no CDS", anchor="", identity="",
                             rdrp_nt="", vp1_nt="", note="record has no CDS feature"))
            continue
        prot = protein(cds, r)
        best = None
        for acc, aprot, abounds in both:
            if acc == r.id:
                continue
            mapped, ident = map_bounds(aprot, abounds, prot, aln)
            if best is None or ident > best[2]:
                best = (acc, mapped, ident)
        acc, mapped, ident = best
        off = int(cds.location.start) + codon_start(cds) - 1
        note = []
        got = {}
        if ident < a.min_identity:
            note.append(f"anchor identity {ident:.2f} below {a.min_identity}: not extracted")
            mapped = {}
        for region, (s, e) in mapped.items():
            aa_len = e - s
            lo, hi = VP1_AA_RANGE if region == "vp1" else RDRP_AA_RANGE
            if not (lo <= aa_len <= hi):
                note.append(f"{region} {aa_len} aa outside {lo}-{hi}")
                continue
            nt0, nt1 = off + s * 3, off + e * 3
            sub = str(r.seq[nt0:nt1])
            got[region] = (nt0 + 1, nt1)          # 1-based inclusive, for the TSV
            seqs[region].append((r.id, sub))
        # measurable accuracy: this record's own annotation, if it has one
        for region, words in (("vp1", VP1_WORDS), ("rdrp", RDRP_WORDS)):
            mp = mat_peptides(r, words)
            if mp and region in mapped:
                truth = aa_bounds(cds, mp[0])
                if truth:
                    checks.append((r.id, region, mapped[region][0] - truth[0],
                                   mapped[region][1] - truth[1], round(ident, 3)))
        rows.append(dict(accession=r.id, status="ok" if got else "failed",
                         anchor=acc, identity=round(ident, 3),
                         rdrp_nt="-".join(map(str, got["rdrp"])) if "rdrp" in got else "",
                         vp1_nt="-".join(map(str, got["vp1"])) if "vp1" in got else "",
                         note="; ".join(note)))

    for region in ("rdrp", "vp1"):
        p = out / f"sapovirus_{region}.fasta"
        with open(p, "w") as fh:
            for name, s in seqs[region]:
                fh.write(f">{name}\n{s}\n")
        print(f"{p.name}: {len(seqs[region])} sequences")

    tsv = out / "coordinates.tsv"
    cols = ["accession", "status", "anchor", "identity", "rdrp_nt", "vp1_nt", "note"]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{tsv.name}: {sum(r['status'] == 'ok' for r in rows)} ok, "
          f"{sum(r['status'] != 'ok' for r in rows)} not extracted")

    if checks:
        print("\nmapping error against records that annotate the region themselves:")
        print("  accession    region  start_delta  end_delta  identity")
        for acc, region, ds, de, ident in checks:
            print(f"  {acc:12} {region:6} {ds:+11d} {de:+10d}  {ident}")
        ds = [abs(c[2]) for c in checks]
        print(f"  median |start error| = {sorted(ds)[len(ds) // 2]} aa over {len(checks)} checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
