"""Typing pipeline: BLAST genogroup -> region check (>=100 nt) -> ClustalW profile alignment -> NJ + bootstrap."""
import os, subprocess, tempfile
from functools import lru_cache
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from . import config as cg
from . import tree

MIN_OVERLAP = 100
MIN_QUERY_FRAC_COLS = 0.5   # keep alignment columns where >=50% of references have a base


@lru_cache(maxsize=None)
def load_reference(region, group):
    """Return (aln_path, {accession: type}) for a region ('orf1'/'orf2') and genogroup."""
    d = cg.REF_DIR[region]
    aln = d / f"{group}.aln"
    if not aln.exists():
        aln, orf = d / "all.aln", d / "all.orf"
    else:
        orf = d / f"{group}.orf"
    types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
    return str(aln), types


def _clustalw_profile(aln, query_fa, out_fa):
    cmd = [cg.exe("clustalw"), f"-PROFILE1={aln}", f"-PROFILE2={query_fa}", "-SEQUENCES",
           f"-OUTFILE={out_fa}", "-OUTPUT=FASTA"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def type_region(seq, group, region, exclude=(), nrep=100, seed=1, min_support=70.0, max_refs_per_type=None):
    """Type one sequence for one region against the reference set of `group`."""
    aln, types = load_reference(region, group)
    excl = set(exclude)
    types = {k: v for k, v in types.items() if k not in excl}
    with tempfile.TemporaryDirectory(prefix="norotyper_") as td:
        q, out = os.path.join(td, "q.fa"), os.path.join(td, "out.fa")
        with open(q, "w") as fh:
            fh.write(f">QUERY\n{str(seq).upper()}\n")
        ref_aln = aln
        if exclude:   # leave-one-out: write the reference alignment without the excluded sequences
            from Bio import AlignIO
            a = AlignIO.read(aln, "clustal")
            kept = [r for r in a if r.id.split("|")[-1] not in excl]
            ref_aln = os.path.join(td, "ref.aln")
            from Bio.Align import MultipleSeqAlignment
            AlignIO.write(MultipleSeqAlignment(kept), ref_aln, "clustal")
        _clustalw_profile(ref_aln, q, out)
        recs = list(SeqIO.parse(out, "fasta"))
    names = [r.id for r in recs]
    # alignment labels are either a bare accession or 'group|type|accession' depending on the file
    types = {n: types[n.split("|")[-1]] for n in names if n != "QUERY" and n.split("|")[-1] in types}
    seqs = [str(r.seq) for r in recs]
    codes = tree.encode(seqs)
    qi = names.index("QUERY")
    refmask = np.array([i != qi for i in range(len(names))])
    frac = (codes[refmask] < 4).mean(0)
    cols = (codes[qi] < 4) & (frac >= MIN_QUERY_FRAC_COLS)
    if cols.sum() < MIN_OVERLAP:
        return dict(type=None, support=0.0, supports={}, n_sites=int(cols.sum()), n_refs=len(types),
                    note="too few aligned sites")
    res = tree.assign(names, types, "QUERY", codes[:, cols], nrep=nrep, seed=seed, min_support=min_support,
                      max_refs_per_type=max_refs_per_type)
    return res


@lru_cache(maxsize=1)
def _annotation():
    a = pd.read_csv(cg.ANNOTATION)
    a["acc"] = a["id"].astype(str)
    return a.set_index("acc")


def blast_genogroup(query_fa):
    """Return per-query dict(genogroup, subject, sstart, send, reverse)."""
    cmd = [cg.exe("blastn"), "-query", query_fa, "-db", str(cg.BLAST_DB), "-evalue", "1e-5", "-task", "blastn",
           "-max_target_seqs", "5", "-max_hsps", "1",
           "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    best = {}
    for line in out.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ev, bs = line.split("\t")
        if q not in best or float(bs) > best[q]["bitscore"]:
            best[q] = dict(subject=s, sstart=int(ss), send=int(se), bitscore=float(bs), pident=float(pid),
                           aln_len=int(ln), evalue=float(ev))
    for q, h in best.items():
        h["genogroup"] = h["subject"].split("|")[0]
        h["reverse"] = h["sstart"] > h["send"]
    return best


def region_overlap(hit, region):
    ann = _annotation()
    acc = hit["subject"].split("|")[-1].split()[0]
    if acc not in ann.index:
        return 0
    row = ann.loc[acc]
    lo, hi = min(hit["sstart"], hit["send"]), max(hit["sstart"], hit["send"])
    a, b = (row.rdrp_start, row.rdrp_end) if region == "orf1" else (row.vp1_start, row.vp1_end)
    if pd.isna(a) or pd.isna(b):
        return 0
    return max(0, min(hi, b) - max(lo, a))


def type_query(rec_id, seq, hit, nrep=100, seed=1, min_support=70.0, max_refs_per_type=None):
    row = dict(query=rec_id, genogroup="unassigned", orf1_type="", orf1_support="", orf2_type="", orf2_support="",
               note="")
    if hit is None:
        row["note"] = "no BLAST hit (not a norovirus / Caliciviridae sequence in reference set)"
        return row
    g = hit["genogroup"]
    row["genogroup"] = g if g in cg.NOROVIRUS_GROUPS else g
    if g not in cg.NOROVIRUS_GROUPS:
        row["note"] = "other Caliciviridae: genogroup only"
        return row
    if hit["reverse"]:
        seq = str(Seq(str(seq)).reverse_complement())
    notes = []
    for region in ("orf1", "orf2"):
        if region_overlap(hit, region) < MIN_OVERLAP:
            notes.append(f"{region}: <{MIN_OVERLAP} nt overlap")
            continue
        if not (cg.REF_DIR[region] / f"{g}.aln").exists() and not (cg.REF_DIR[region] / "all.aln").exists():
            continue
        res = type_region(seq, g, region, nrep=nrep, seed=seed, min_support=min_support,
                          max_refs_per_type=max_refs_per_type)
        row[f"{region}_type"] = res["type"] or "unassigned"
        row[f"{region}_support"] = round(res["support"], 1)
        if res["type"] is None:
            notes.append(f"{region}: best {res.get('best_type', '-')} {res['support']:.0f}% < {min_support:.0f}%")
    row["note"] = "; ".join(notes)
    return row
