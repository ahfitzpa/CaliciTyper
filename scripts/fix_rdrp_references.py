#!/usr/bin/env python3
"""Correct two GII RdRp references whose typing region was taken from the wrong CDS.

JX846926.1 (GII.P8, Hu/GII.7/CHDC3936/1988/USA) and JX289822.1 (GII.P1,
Hu/GII.1/7EK/Hawaii/1971/USA) are partial genomes whose GenBank records carry
several CDS features. When the reference set was built, ORF1 was annotated as
the record's *first* CDS (1..442 and 1..1555), a 5' fragment, so the "RdRp
typing region" stored for each was unrelated sequence: about 0.60 nt p-distance
from every other reference, including its own type. Consequences: no GII.P8 or
GII.P1 query could reach bootstrap support (the clade test needs both of a
type's references), and the novelty screen flagged the GII.P8 reference strain
AB039780.1 (Saitama U25) as a putative new P-type.

The correct regions are the last 762 nt of each record's ORF1 CDS (the CDS
ending at the ORF1/ORF2 junction): JX846926.1 923..4535 -> typing region
3774..4535, 95.0% identical to the other GII.P8 reference AB039780.1;
JX289822.1 2086..3815 -> 3054..3815, 100% identical to the other GII.P1
reference U07611.2.

This rewrites, for those two accessions only:
  GII_annotation.csv  ORF1 coordinates
  GII.orf             the typing region (762 nt)
  all.orf             the ORF1 CDS
  GII.aln, all.aln    the aligned row, re-added with MAFFT --keeplength so no
                      reference column moves
The .dnd guide trees are legacy files the package does not read; they are left
alone. Needs network access (NCBI efetch) and mafft on PATH.

    python scripts/fix_rdrp_references.py
"""
import csv
import io
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

from Bio import AlignIO, SeqIO

D = Path(__file__).resolve().parent.parent / "calicityper" / "data" / "reference_sequences_RdRp"
FIX = {  # accession: (type, ORF1 CDS start, end), 1-based inclusive, from the GenBank record
    "JX846926.1": ("GII.P8", 923, 4535),
    "JX289822.1": ("GII.P1", 2086, 3815),
}
REGION_NT = 762
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def genome(acc):
    url = f"{EUTILS}?db=nuccore&id={acc}&rettype=fasta&retmode=text&tool=calicityper"
    with urllib.request.urlopen(url, timeout=120) as r:
        return str(SeqIO.read(io.StringIO(r.read().decode()), "fasta").seq).upper()


def rewrite_fasta(path, acc, seq):
    """Replace one record's sequence line in place (files are one line per sequence)."""
    lines = path.read_text().split("\n")
    at = [i for i, x in enumerate(lines) if x.startswith(">") and x[1:].split()[0].split("|")[-1] == acc]
    assert len(at) == 1 and not lines[at[0] + 1].startswith(">"), (path, acc)
    assert at[0] + 2 >= len(lines) or lines[at[0] + 2].startswith(">") or lines[at[0] + 2] == "", \
        f"{path}: {acc} is not on a single line"
    lines[at[0] + 1] = seq
    path.write_text("\n".join(lines))


def rewrite_aln(path, acc, region):
    aln = AlignIO.read(path, "clustal")
    ids = [r.id for r in aln]
    k = [i for i, x in enumerate(ids) if x.split("|")[-1] == acc]
    assert len(k) == 1, (path, acc)
    k = k[0]
    rest = [r for i, r in enumerate(aln) if i != k]
    with tempfile.TemporaryDirectory() as td:
        rf, qf = Path(td) / "ref.fa", Path(td) / "q.fa"
        SeqIO.write(rest, rf, "fasta")
        qf.write_text(f">q\n{region}\n")
        out = subprocess.run(["mafft", "--quiet", "--keeplength", "--addfragments", str(qf), str(rf)],
                             check=True, capture_output=True, text=True).stdout
    row = [r for r in SeqIO.parse(io.StringIO(out), "fasta") if r.id == "q"][0]
    row = str(row.seq).upper()
    assert len(row) == aln.get_alignment_length()
    dropped = len(region) - len(row.replace("-", ""))
    # edit the row's segments in place, block by block, so the file's layout
    # (and every other row) is untouched; the per-block conservation lines are
    # left as they were
    lines = path.read_text().split("\n")
    pos = 0
    for i, x in enumerate(lines):
        if x.split(" ", 1)[0] == ids[k]:
            name_w = len(x) - len(x.split()[-1])
            seg = x.split()[-1]
            lines[i] = x[:name_w] + row[pos:pos + len(seg)]
            pos += len(seg)
    assert pos == len(row), (path, acc, pos, len(row))
    path.write_text("\n".join(lines))
    return dropped


def main():
    ann = D / "GII_annotation.csv"
    rows = list(csv.DictReader(open(ann)))
    for acc, (ptype, start, end) in FIX.items():
        g = genome(acc)
        cds = g[start - 1:end]
        region = cds[-REGION_NT:]
        for r in rows:
            if r["id"] == acc:
                r["rdrp_start"], r["rdrp_end"] = str(start), str(end)
        rewrite_fasta(D / "GII.orf", acc, region)
        rewrite_fasta(D / "all.orf", acc, cds)
        d1 = rewrite_aln(D / "GII.aln", acc, region)
        d2 = rewrite_aln(D / "all.aln", acc, region)
        print(f"{acc} ({ptype}): ORF1 {start}..{end}, typing region {end - REGION_NT + 1}..{end}; "
              f"bases dropped by --keeplength: GII.aln {d1}, all.aln {d2}")
    with open(ann, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
