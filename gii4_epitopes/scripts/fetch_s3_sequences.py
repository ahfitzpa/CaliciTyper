#!/usr/bin/env python3
"""Fetch the Parra et al. S3 Table records from GenBank and extract VP1.

Writes one protein FASTA, cached, so the clustering step can be re-run
without touching the network. Sequences are identified as
`<accession>|<genotype>` using the genotype column of S3 Table, which is the
authors' own label.

VP1 is taken from the CDS annotated as VP1 / capsid / ORF2 where present; a
record with no such CDS is reported and skipped rather than guessed at.

No personal address is sent: NCBI asks automated clients for a contact, and
this uses a generic one unless NCBI_EMAIL says otherwise, matching
fetch_rivm_test_set.py.

    python scripts/fetch_s3_sequences.py
"""
from __future__ import annotations

import csv
import io
import os
import sys
import time
import urllib.parse
import urllib.request

from Bio import SeqIO

from common import INPUTS

SRC = INPUTS / "parra_s3" / "accessions.tsv"
OUT = INPUTS / "parra_s3" / "vp1_proteins.faa"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "calicityper-immunotype-rebuild"
WANT = ("vp1", "capsid", "orf2", "major capsid")


def fetch(ids, email, tries=4):
    q = {"db": "nuccore", "rettype": "gb", "retmode": "text",
         "tool": TOOL, "email": email, "id": ",".join(ids)}
    url = EUTILS + "?" + urllib.parse.urlencode(q)
    for n in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf8", "replace")
        except Exception as e:
            if n == tries:
                raise
            print(f"    {type(e).__name__}: {e}; retry {n}/{tries}", flush=True)
            time.sleep(10 * n)


def vp1_protein(rec):
    """(aa, how) for the VP1 CDS, or (None, reason)."""
    for f in rec.features:
        if f.type != "CDS":
            continue
        labels = " ".join(str(v).lower() for k in ("product", "gene", "note")
                          for v in f.qualifiers.get(k, []))
        if any(w in labels for w in WANT):
            aa = f.qualifiers.get("translation", [None])[0]
            if aa:
                return aa, "annotated CDS translation"
            try:
                return str(f.extract(rec.seq).translate(to_stop=True)), "translated CDS"
            except Exception as e:
                return None, f"CDS present but untranslatable: {e}"
    return None, "no VP1/capsid/ORF2 CDS annotated"


def main():
    if not SRC.exists():
        sys.exit(f"missing {SRC}; run scripts/parse_s3_table.py first")
    rows = list(csv.DictReader(open(SRC), delimiter="\t"))
    email = os.environ.get("NCBI_EMAIL", "calicityper@example.org")

    if OUT.exists():
        have = {r.id.split("|")[0] for r in SeqIO.parse(OUT, "fasta")}
        if have >= {r["accession"] for r in rows}:
            print(f"{OUT.name} already has all {len(have)} records")
            return 0

    gt = {r["accession"]: r["genotype"] for r in rows}
    ids = sorted(gt)
    got, skipped = [], []
    for i in range(0, len(ids), 40):
        batch = ids[i:i + 40]
        print(f"  fetching {i + 1}-{i + len(batch)} of {len(ids)}", flush=True)
        text = fetch(batch, email)
        for rec in SeqIO.parse(io.StringIO(text), "genbank"):
            acc = rec.id.split(".")[0]
            aa, how = vp1_protein(rec)
            if aa and len(aa) >= 400:
                got.append((f"{acc}|{gt.get(acc, '?')}", aa))
            else:
                skipped.append((acc, gt.get(acc, "?"),
                                how if not aa else f"VP1 only {len(aa)} aa"))
        time.sleep(1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(f">{n}\n{s}\n" for n, s in got))
    print(f"\n{len(got)} VP1 proteins -> {OUT}")
    if skipped:
        print(f"{len(skipped)} skipped:")
        for acc, g, why in skipped[:15]:
            print(f"  {acc:10s} {g:8s} {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
