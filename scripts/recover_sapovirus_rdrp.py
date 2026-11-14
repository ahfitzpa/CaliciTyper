#!/usr/bin/env python3
"""Recover complete-RdRp references for P-types that have none.

16 of the 45 P-types in Zhao et al. 2024's table have no strain carrying the
full canonical RdRp window: 43 of their 45 listed strains are 5' truncated,
clustering at 283-294 aa because they are 3'-end submissions covering the
typing region rather than whole genomes.

Rather than relax the window -- a distance over 290 aa is not comparable with
cutoffs derived on complete sequences -- this searches GenBank's near-complete
sapovirus genomes for records that carry the full window AND match a listed
partial reference closely enough to inherit its P-type.

The P-type still comes from the published table. A candidate is only accepted
when its RdRp contains the listed partial reference at high nucleotide
identity, so the recovered record is the same P-type as the strain the paper
designated, not a nearest neighbour.

    python scripts/recover_sapovirus_rdrp.py \\
        --genbank <existing.gb> --out-dir <dir> [--min-identity 0.95]
"""
from __future__ import annotations

import argparse
import csv
import io
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import importlib.util as _u                                   # noqa: E402

_spec = _u.spec_from_file_location("exsap", ROOT / "scripts" / "extract_sapovirus_regions.py")
ex = _u.module_from_spec(_spec)
_spec.loader.exec_module(ex)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = "Sapovirus[Organism] AND 6500:8000[SLEN]"


def esearch(term, retmax=2000):
    url = (f"{EUTILS}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(term)}"
           f"&retmax={retmax}&retmode=json&tool=calicityper")
    import json
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)["esearchresult"].get("idlist", [])


def efetch_gb(ids, batch=50, pause=0.4):
    """GenBank flatfiles, in batches, by GET (efetch drops batched POSTs)."""
    out = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(chunk)}"
               f"&rettype=gb&retmode=text&tool=calicityper")
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=180) as r:
                    text = r.read().decode("utf-8", "replace")
                out.extend(SeqIO.parse(io.StringIO(text), "genbank"))
                break
            except Exception as e:                            # noqa: BLE001
                print(f"    retry {attempt + 1}: {type(e).__name__}: {e}", flush=True)
                time.sleep(2 ** attempt)
        print(f"  fetched {min(i + batch, len(ids))}/{len(ids)}", flush=True)
        time.sleep(pause)
    return out


def nt_aligner():
    a = PairwiseAligner(mode="local")
    a.match_score, a.mismatch_score = 2.0, -3.0
    a.open_gap_score, a.extend_gap_score = -7.0, -2.0
    return a


def contains_identity(short, long_seq, aln):
    """Identity of `short` where it aligns inside `long_seq`, and its coverage."""
    if not short or not long_seq:
        return 0.0, 0.0
    al = aln.align(short.upper(), long_seq.upper())[0]
    n = same = 0
    for (a0, a1), (b0, b1) in zip(*al.aligned):
        n += a1 - a0
        same += sum(x == y for x, y in zip(short[a0:a1].upper(), long_seq[b0:b1].upper()))
    return (same / n if n else 0.0), (n / len(short))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genbank", required=True, help="the records already fetched")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-identity", type=float, default=0.95)
    ap.add_argument("--min-coverage", type=float, default=0.90)
    ap.add_argument("--limit", type=int, default=0, help="cap candidates, for testing")
    a = ap.parse_args()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    have = list(SeqIO.parse(a.genbank, "genbank"))
    have_ids = {r.id.split(".")[0] for r in have}

    # P-types with no complete-window strain, and the partial references that
    # define them
    tab = defaultdict(list)
    for r in csv.DictReader(l for l in open(ROOT / "validation/sapovirus_rdrp_taxonomy.csv")
                            if not l.startswith("#")):
        if not r["note"].startswith("UNRESOLVED"):
            tab[r["ptype"]].append(r["accession"])
    complete = {r.id.split(".")[0] for r in
                SeqIO.parse(out.parent / "regions" / "sapovirus_rdrp.fasta", "fasta")} \
        if (out.parent / "regions" / "sapovirus_rdrp.fasta").exists() else set()
    missing = {p: accs for p, accs in tab.items()
               if not ({x.split(".")[0] for x in accs} & complete)}
    print(f"{len(missing)} P-types need a complete-window reference")

    # the partial RdRp of each missing P-type's strains, as nucleotides
    aln_aa = ex.aligner()
    anchors = ex.project_canonical(have, aln_aa)
    by_id = {r.id.split(".")[0]: r for r in have}
    probes = {}
    for p, accs in missing.items():
        for acc in accs:
            rec = by_id.get(acc.split(".")[0])
            if rec is None:
                continue
            cds = ex.orf1_cds(rec)
            if cds is None:
                continue
            prot = ex.protein(cds, rec)
            mapped, ident = ex.map_bounds(anchors[0][1], ex.CANONICAL_BOUNDS, prot, aln_aa)
            if "rdrp" not in mapped:
                continue
            off = int(cds.location.start) + ex.codon_start(cds) - 1
            s, e = mapped["rdrp"]
            probes.setdefault(p, []).append((acc, str(rec.seq[off + s * 3: off + e * 3])))
    print(f"partial RdRp probes available for {len(probes)} of them")

    ids = esearch(QUERY)
    print(f"GenBank near-complete sapovirus genomes: {len(ids)}")
    cand_ids = [i for i in ids]
    if a.limit:
        cand_ids = cand_ids[:a.limit]
    cands = efetch_gb(cand_ids)
    cands = [r for r in cands if r.id.split(".")[0] not in have_ids]
    print(f"{len(cands)} candidates after dropping records already held")

    # extract the canonical window from each candidate
    extracted = {}
    for rec in cands:
        cds = ex.orf1_cds(rec)
        if cds is None:
            continue
        prot = ex.protein(cds, rec)
        best = None
        for acc, aprot, abounds in anchors:
            mapped, ident = ex.map_bounds(aprot, abounds, prot, aln_aa)
            if best is None or ident > best[1]:
                best = (mapped, ident)
            if ident >= 0.60:
                break
        mapped, ident = best
        if "rdrp" not in mapped or ident < 0.30:
            continue
        s, e = mapped["rdrp"]
        if not (ex.RDRP_AA_RANGE[0] <= e - s <= ex.RDRP_AA_RANGE[1]):
            continue
        off = int(cds.location.start) + ex.codon_start(cds) - 1
        extracted[rec.id] = (str(rec.seq[off + s * 3: off + e * 3]), ident, rec.description[:60])
    print(f"{len(extracted)} candidates carry the full canonical window")

    naln = nt_aligner()
    rows, kept = [], {}
    for p, plist in sorted(probes.items()):
        best = None
        for probe_acc, probe in plist:
            for cid, (seq, ident, desc) in extracted.items():
                pid, cov = contains_identity(probe, seq, naln)
                if cov < a.min_coverage or pid < a.min_identity:
                    continue
                if best is None or pid > best["identity"]:
                    best = dict(ptype=p, accession=cid, identity=round(pid, 4),
                                coverage=round(cov, 4), matched_partial=probe_acc,
                                description=desc, nt=len(seq))
        if best:
            kept[p] = (best["accession"], extracted[best["accession"]][0])
            rows.append(best)
            print(f"  {p:10} <- {best['accession']:12} id={best['identity']:.3f} "
                  f"cov={best['coverage']:.2f}  ({best['description'][:40]})")
        else:
            rows.append(dict(ptype=p, accession="", identity="", coverage="",
                             matched_partial="", description="no complete-window record "
                             "matched a listed strain", nt=""))
    with open(out / "recovered_ptypes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "ptype", "accession", "identity", "coverage", "matched_partial", "nt", "description"])
        w.writeheader()
        w.writerows(rows)
    with open(out / "recovered_rdrp.fasta", "w") as fh:
        for p, (acc, seq) in sorted(kept.items()):
            fh.write(f">{acc} {p}\n{seq}\n")
    print(f"\nrecovered {len(kept)} of {len(missing)} missing P-types "
          f"-> {out / 'recovered_rdrp.fasta'}")
    print(f"still missing: {sorted(set(missing) - set(kept))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
