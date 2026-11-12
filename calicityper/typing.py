"""Typing pipeline: BLAST genogroup -> region check (>=100 nt) -> ClustalW profile alignment -> NJ + bootstrap."""
import os, subprocess, tempfile
from functools import lru_cache
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from . import config as cg
from . import tree, align

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


def _type_region_clustalw(seq, group, region, exclude=(), nrep=100, seed=1, min_support=70.0, max_refs_per_type=None):
    """Legacy path: ClustalW profile alignment, one query per call (slow)."""
    aln, types = load_reference(region, group)
    excl = set(exclude)
    types = {k: v for k, v in types.items() if k not in excl}
    with tempfile.TemporaryDirectory(prefix="calicityper_") as td:
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


def type_aligned(ref, qaln, exclude=(), nrep=100, seed=1, min_support=70.0, max_refs_per_type=None):
    """Bootstrap-NJ typing of one query already aligned to `ref` (align.RefSet)."""
    ex = set(exclude)
    keep = [i for i, n in enumerate(ref.names) if n not in ex]
    names = [ref.names[i] for i in keep] + ["QUERY"]
    types = {ref.names[i]: ref.types[i] for i in keep}
    codes = np.vstack([ref.codes[keep], tree.encode([qaln])])
    frac = (codes[:-1] < 4).mean(0)
    cols = (codes[-1] < 4) & (frac >= MIN_QUERY_FRAC_COLS)
    if cols.sum() < MIN_OVERLAP:
        return dict(type=None, support=0.0, supports={}, n_sites=int(cols.sum()), n_refs=len(types),
                    note="too few aligned sites")
    return tree.assign(names, types, "QUERY", codes[:, cols], nrep=nrep, seed=seed, min_support=min_support,
                       max_refs_per_type=max_refs_per_type)


def type_region(seq, group, region, exclude=(), nrep=100, seed=1, min_support=70.0, max_refs_per_type=None,
                backend="mafft", threads=1):
    """Type one sequence for one region against the reference set of `group`."""
    if backend == "clustalw":
        return _type_region_clustalw(seq, group, region, exclude, nrep, seed, min_support, max_refs_per_type)
    ref, (qaln,) = align.align_queries(region, group, [("q", seq)], exclude, threads)
    return type_aligned(ref, qaln, exclude, nrep, seed, min_support, max_refs_per_type)


@lru_cache(maxsize=1)
def _annotation():
    a = pd.read_csv(cg.ANNOTATION)
    a["acc"] = a["id"].astype(str)
    return a.set_index("acc")


def blast_genogroup(query_fa, threads=1):
    """Return per-query dict(genogroup, subject, sstart, send, reverse)."""
    cmd = [cg.exe("blastn"), "-query", query_fa, "-db", str(cg.BLAST_DB), "-evalue", "1e-5", "-task", "blastn",
           "-max_target_seqs", "5", "-max_hsps", "1", "-num_threads", str(max(1, threads)),
           "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    best = {}
    for line in out.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ev, bs = line.split("\t")
        if q not in best or float(bs) > best[q]["bitscore"]:
            best[q] = dict(subject=s, qstart=int(qs), qend=int(qe), sstart=int(ss), send=int(se), bitscore=float(bs), pident=float(pid),
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
    row = dict(query=rec_id, genogroup="unassigned", orf1_type="", orf1_support="", orf1_runners="",
               orf2_type="", orf2_support="", orf2_runners="", note="")
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


CROP_PAD = 300   # nt kept on each side of the expected ORF when cropping a long query before alignment


def orf_window(hit, region):
    """Approximate window (0-based, half-open) of `region` on the ORIGINAL query, mapped through the BLAST HSP."""
    ann = _annotation()
    acc = hit["subject"].split("|")[-1].split()[0]
    if acc not in ann.index:
        return None
    row = ann.loc[acc]
    a, b = (row.rdrp_start, row.rdrp_end) if region == "orf1" else (row.vp1_start, row.vp1_end)
    if pd.isna(a) or pd.isna(b):
        return None
    f = (lambda x: hit["qstart"] + (hit["sstart"] - x)) if hit["reverse"] else (lambda x: hit["qstart"] + (x - hit["sstart"]))
    lo, hi = sorted((f(a), f(b)))
    return max(0, int(lo) - CROP_PAD - 1), int(hi) + CROP_PAD


NUCLEOTIDES = frozenset("ACGTURYKMSWBDHVN-")     # IUPAC codes plus gap


def invalid_characters(seq):
    """Characters in seq that are not IUPAC nucleotide codes (case-insensitive), sorted."""
    return sorted(set(str(seq).upper()) - NUCLEOTIDES)


def plan_query(rec_id, seq, hit):
    """Decide what to do with one query. Returns (row, jobs) where jobs = [(region, group, oriented_seq)]."""
    row = dict(query=rec_id, genogroup="unassigned", orf1_type="", orf1_support="", orf1_runners="",
               orf2_type="", orf2_support="", orf2_runners="", note="")
    bad = invalid_characters(seq)
    if bad:
        # Typically two FASTA records merged onto one line (a header inside the
        # sequence). MAFFT rejects these, and one such query would otherwise take
        # down the whole alignment chunk it shares with other queries.
        shown = "".join(bad[:10]) + ("..." if len(bad) > 10 else "")
        row["note"] = f"not typed: sequence contains non-nucleotide characters ({shown!r}); check the FASTA is well formed"
        return row, []
    ungapped = len(str(seq).replace("-", "").replace("N", ""))
    if ungapped < cg.MIN_INPUT_NT:
        row["note"] = (f"rejected: {ungapped} nt of sequence is shorter than the shortest "
                       f"RIVM typing window ({cg.MIN_INPUT_NT} nt), so no window can be covered")
        return row, []
    if hit is None:
        row["note"] = "no BLAST hit (not a norovirus / Caliciviridae sequence in reference set)"
        return row, []
    g = hit["genogroup"]
    row["genogroup"] = g
    if g not in cg.NOROVIRUS_GROUPS:
        row["note"] = "other Caliciviridae: genogroup only"
        return row, []
    if g not in cg.TYPEABLE_GROUPS:
        row["note"] = (f"norovirus {g}: genogroup only, genotyping is validated for "
                       f"{'/'.join(cg.TYPEABLE_GROUPS)} only")
        return row, []
    jobs, notes = [], []
    for region in ("orf1", "orf2"):
        if region_overlap(hit, region) < MIN_OVERLAP:
            notes.append(f"{region}: <{MIN_OVERLAP} nt overlap")
            continue
        if not (cg.REF_DIR[region] / f"{g}.aln").exists() and not (cg.REF_DIR[region] / "all.aln").exists():
            continue
        sub = str(seq)
        w = orf_window(hit, region)
        if w:                                   # long query (e.g. whole genome): align only the ORF neighbourhood
            sub = sub[w[0]:w[1]]
        if hit["reverse"]:
            sub = str(Seq(sub).reverse_complement())
        jobs.append((region, g, sub))
    row["note"] = "; ".join(notes)
    return row, jobs


MAX_RUNNERS = 3          # how many ranked types to report per region


def format_runners(supports, n=MAX_RUNNERS):
    """'GII.P31:97.0|GII.P16:2.0' -- the n best-supported types, best first.

    assign() already ranks these; they were computed and then discarded.
    Useful for near-threshold calls, where the runner-up is the real question.
    """
    ranked = sorted((supports or {}).items(), key=lambda kv: (-kv[1], kv[0]))[:n]
    return "|".join(f"{t}:{v:.1f}" for t, v in ranked)


def apply_result(row, region, res, min_support):
    row[f"{region}_type"] = res["type"] or "unassigned"
    row[f"{region}_support"] = round(res["support"], 1)
    row[f"{region}_runners"] = format_runners(res.get("supports"))
    if res["type"] is None:
        extra = f"{region}: best {res.get('best_type', '-')} {res['support']:.0f}% < {min_support:.0f}%"
        row["note"] = "; ".join(x for x in (row["note"], extra) if x)


# ---------------------------------------------- RIVM operational typing windows

def rivm_window_columns(region, group, genotype=None):
    """Alignment columns of the RIVM window for `region`, as a numpy index array.

    ORF2: the first N reference-covered columns of VP1 (region C, measured from
    the 5' end of VP1). ORF1: the last N columns of the alignment, which already
    spans the 3' RdRp typing region.

    Columns are chosen by reference occupancy so that an insertion present in
    only a few references cannot shift the window -- the reason windows are
    column ranges here and not genomic coordinates.
    """
    import numpy as np
    from . import align
    lengths = cg.RIVM_WINDOW_NT[region]
    n = lengths.get(genotype) or lengths.get(group) or lengths["_default"]
    ref = align.refset(region, group)
    covered = np.flatnonzero(ref.occupancy >= 0.5)
    if len(covered) == 0:
        return np.array([], dtype=int), n
    return (covered[:n] if region == "orf2" else covered[-n:]), n


def rivm_window_coverage(region, group, aligned_query, genotype=None):
    """(covered nt, window nt, fraction) for one query over its RIVM window."""
    import numpy as np
    cols, n = rivm_window_columns(region, group, genotype)
    if len(cols) == 0:
        return 0, n, 0.0
    q = np.frombuffer(aligned_query.upper().encode(), dtype="S1")[cols]
    covered = int(np.isin(q, [b"A", b"C", b"G", b"T", b"U"]).sum())
    return covered, len(cols), round(covered / len(cols), 4)
