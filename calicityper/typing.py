"""Typing pipeline: BLAST genogroup -> region check (>=100 nt) -> ClustalW profile alignment -> NJ + bootstrap."""
import csv, math, os, subprocess, tempfile
from functools import lru_cache
import numpy as np
from . import config as cg
from . import tree, align

MIN_OVERLAP = 100
MIN_QUERY_FRAC_COLS = 0.5   # keep alignment columns where >=50% of references have a base

_RC_TABLE = str.maketrans(
    "ACGTURYKMSWBDHVNacgturykmswbdhvn-",
    "TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn-",
)


def _revcomp(seq):
    return str(seq).translate(_RC_TABLE)[::-1]


def read_fasta(path: str | os.PathLike[str]) -> list[tuple[str, str]]:
    """Fast stdlib FASTA reader returning [(record_id, seq_str)] matching SeqIO.parse."""
    with open(path, "rb") as raw_fh:
        if raw_fh.read(2) == b"\x1f\x8b":
            raise ValueError(
                f"{path} appears to be gzip-compressed (.gz); decompress it first "
                "(e.g. `gunzip -k <file>`)."
            )
    recs: list[tuple[str, str]] = []
    rid: str | None = None
    buf: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    recs.append((rid, "".join(buf)))
                parts = line[1:].strip().split(None, 1)
                rid = parts[0] if parts else ""
                buf = []
            elif rid is not None:
                buf.append(line)
    if rid is not None:
        recs.append((rid, "".join(buf)))
    return recs


@lru_cache(maxsize=None)
def load_reference(region, group):
    """Return (aln_path, {accession: type}) for a region ('orf1'/'orf2') and genogroup."""
    aln, orf = cg.ref_paths(region, group)
    types = {rid.split("|")[2]: rid.split("|")[1] for rid, _ in read_fasta(orf)}
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
        recs = read_fasta(out)
    names = [rid for rid, _ in recs]
    # alignment labels are either a bare accession or 'group|type|accession' depending on the file
    types = {n: types[n.split("|")[-1]] for n in names if n != "QUERY" and n.split("|")[-1] in types}
    seqs = [sq for _, sq in recs]
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
    acc, ident, n_ident = tree.nearest_reference(names, "QUERY", codes[:, cols])
    best_t = types.get(acc, "")
    sec_t, sec_d, sec_n = tree.second_nearest_type(names, types, "QUERY", codes[:, cols], best_t)
    res.update(nearest=acc.split("|")[-1], nearest_identity=ident, nearest_sites=n_ident,
               nearest_type=best_t, second_type=sec_t, second_distance=sec_d, second_sites=sec_n)
    return res


def type_aligned(ref, qaln, exclude=(), nrep=100, seed=1, min_support=70.0,
                 max_refs_per_type=None, fast_boot=False):
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
                    region_sites=int((frac >= MIN_QUERY_FRAC_COLS).sum()), note="too few aligned sites")
    res = tree.assign(names, types, "QUERY", codes[:, cols], nrep=nrep, seed=seed, min_support=min_support,
                      max_refs_per_type=max_refs_per_type, fast_boot=fast_boot)
    res["region_sites"] = int((frac >= MIN_QUERY_FRAC_COLS).sum())   # columns the region spans
    acc, ident, n_ident = tree.nearest_reference(names, "QUERY", codes[:, cols])
    best_t = types.get(acc, "")
    sec_t, sec_d, sec_n = tree.second_nearest_type(names, types, "QUERY", codes[:, cols], best_t)
    res.update(nearest=acc, nearest_identity=ident, nearest_sites=n_ident,
               nearest_type=best_t, second_type=sec_t, second_distance=sec_d, second_sites=sec_n)
    return res


def type_region(seq, group, region, exclude=(), nrep=100, seed=1, min_support=70.0, max_refs_per_type=None,
                backend="mafft", threads=1, fast_boot=False):
    """Type one sequence for one region against the reference set of `group`."""
    if backend == "clustalw":
        return _type_region_clustalw(seq, group, region, exclude, nrep, seed, min_support, max_refs_per_type)
    ref, (qaln,) = align.align_queries(region, group, [("q", seq)], exclude, threads)
    return type_aligned(ref, qaln, exclude, nrep, seed, min_support, max_refs_per_type, fast_boot=fast_boot)


def _opt_coord(val):
    s = (val or "").strip()
    if not s or s.lower() in ("na", "nan", "none"):
        return None
    v = float(s)
    return None if math.isnan(v) else v


@lru_cache(maxsize=1)
def _annotation():
    """Return {acc: (rdrp_start, rdrp_end, vp1_start, vp1_end)} parsed via stdlib csv."""
    out = {}
    with open(cg.ANNOTATION, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            acc = str(row.get("id", "")).strip()
            if not acc:
                continue
            out[acc] = (
                _opt_coord(row.get("rdrp_start")),
                _opt_coord(row.get("rdrp_end")),
                _opt_coord(row.get("vp1_start")),
                _opt_coord(row.get("vp1_end")),
            )
    return out


def blast_genogroup(query_fa, threads=1):
    """Return per-query dict(genogroup, subject, sstart, send, reverse)."""
    cmd = [cg.exe("blastn"), "-query", query_fa, "-db", str(cg.BLAST_DB), "-evalue", "1e-5", "-task", "blastn",
           "-max_target_seqs", "5", "-max_hsps", "1", "-dbsize", str(cg.BLAST_DBSIZE),
           "-num_threads", str(max(1, threads)),
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
    coords = ann.get(acc)
    if coords is None:
        return 0
    lo, hi = min(hit["sstart"], hit["send"]), max(hit["sstart"], hit["send"])
    a, b = (coords[0], coords[1]) if region == "orf1" else (coords[2], coords[3])
    if a is None or b is None:
        return 0
    return max(0, min(hi, b) - max(lo, a))


def type_query(rec_id, seq, hit, nrep=100, seed=1, min_support=70.0, max_refs_per_type=None, fast_boot=False):
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
        seq = _revcomp(seq)
    notes = []
    for region in ("orf1", "orf2"):
        if region_overlap(hit, region) < MIN_OVERLAP:
            notes.append(f"{region}: <{MIN_OVERLAP} nt overlap")
            continue
        if cg.ref_paths(region, g) is None:
            continue
        res = type_region(seq, g, region, nrep=nrep, seed=seed, min_support=min_support,
                          max_refs_per_type=max_refs_per_type, fast_boot=fast_boot)
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
    coords = ann.get(acc)
    if coords is None:
        return None
    a, b = (coords[0], coords[1]) if region == "orf1" else (coords[2], coords[3])
    if a is None or b is None:
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
    if hit["aln_len"] < cg.MIN_GENOGROUP_HIT_NT:
        row["note"] = (f"no genogroup: best BLAST hit ({hit['genogroup']}) aligns only "
                       f"{hit['aln_len']} nt, below {cg.MIN_GENOGROUP_HIT_NT} nt, so it is a chance match")
        return row, []
    g = hit["genogroup"]
    row["genogroup"] = g
    if cg.is_sapovirus(g):
        if not cg.TYPE_SAPOVIRUS:
            row["note"] = ("sapovirus: genogroup only (sapovirus genotyping is not yet benchmarked; "
                           "--type-sapovirus enables it)")
            return row, []
    elif g not in cg.NOROVIRUS_GROUPS:
        row["note"] = "other Caliciviridae: genogroup only"
        return row, []
    elif g not in cg.TYPEABLE_GROUPS:
        row["note"] = (f"norovirus {g}: genogroup only, genotyping is validated for "
                       f"{'/'.join(cg.TYPEABLE_GROUPS)} only")
        return row, []
    jobs, notes = [], []
    for region in ("orf1", "orf2"):
        if region_overlap(hit, region) < MIN_OVERLAP:
            notes.append(f"{region}: <{MIN_OVERLAP} nt overlap")
            continue
        if cg.ref_paths(region, g) is None:
            if cg.is_sapovirus(g):
                label = "RdRp" if region == "orf1" else "VP1"
                notes.append(f"{region}: no sapovirus {g[len(cg.SAPOVIRUS_PREFIX):]} {label} reference set")
            continue
        sub = str(seq)
        w = orf_window(hit, region)
        if w:                                   # long query (e.g. whole genome): align only the ORF neighbourhood
            sub = sub[w[0]:w[1]]
        if hit["reverse"]:
            sub = _revcomp(sub)
        jobs.append((region, g, sub))
    row["note"] = "; ".join(notes)
    return row, jobs


def dual_type(genotype, ptype):
    """The dual-typing designation: VP1 genotype first, P-type in brackets.

    Follows Zhao et al. 2024 for sapovirus and the established norovirus
    convention: the P-group is written out only when it differs from the VP1
    genogroup, so GI.5 with GI.P5 is "GI.5[P5]" while GIV.1 with GII.P10 --
    whose RdRp sits in a different group from its capsid -- is
    "GIV.1[GII.P10]". That asymmetry is the point of dual typing: it makes an
    incongruent pairing visible in the name itself.
    """
    g, p = (genotype or "").strip(), (ptype or "").strip()
    if not g or g in ("unassigned", "untypeable") or not p or p in ("unassigned", "untypeable"):
        return ""
    # The GII.4 variant is reported in its own right (orf2_type, gii4_window_cov);
    # the dual-typing convention is the genotype stem, GII.4[P16] not
    # GII.4_Sydney[P16].
    g = g.split("_", 1)[0]
    gg = g.split(".")[0]                      # VP1 genogroup, e.g. GI from GI.5
    pg, _, rest = p.partition(".P")           # P-group and P-type number
    if rest and pg == gg:
        return f"{g}[P{rest}]"                # same group: short form
    return f"{g}[{p}]"                        # different group: write it out


MAX_RUNNERS = 3          # how many ranked types to report per region


def format_runners(supports, n=MAX_RUNNERS):
    """'GII.P31:97.0|GII.P16:2.0' -- the n best-supported types, best first.

    assign() already ranks these; they were computed and then discarded.
    Useful for near-threshold calls, where the runner-up is the real question.
    """
    ranked = sorted((supports or {}).items(), key=lambda kv: (-kv[1], kv[0]))[:n]
    return "|".join(f"{t}:{v:.1f}" for t, v in ranked)


def apply_result(row, region, res, min_support):
    """Write one region's result onto the row.

    The best-supported type is ALWAYS reported, and whether it may be relied on
    is carried in its own `_accepted` column rather than by blanking the call.
    Withholding the type hid the one thing a reader needs in order to judge a
    borderline case -- what the tool actually found -- and forced anyone
    wanting it to re-run with a lower threshold. Acceptance has two conditions,
    both measured in this harness: at least `min_support`% bootstrap support,
    and at least cg.MIN_INFORMATIVE_SITES aligned sites, the point at which
    leave-one-out F1 reaches 97% (see validation/benchmarking materials, curated separately).
    """
    if cg.is_sapovirus(row.get("genogroup", "")):
        return _apply_sapovirus_result(row, region, res)
    best = res.get("best_type") or res.get("type") or ""
    sites = int(res.get("n_sites") or 0)
    row[f"{region}_type"] = best or "unassigned"
    row[f"{region}_support"] = round(res["support"], 1)
    row[f"{region}_runners"] = format_runners(res.get("supports"))
    row[f"{region}_sites"] = sites
    row[f"{region}_nearest"] = res.get("nearest", "")
    row[f"{region}_nearest_identity"] = res.get("nearest_identity", "")
    row[f"{region}_nearest_type"] = res.get("nearest_type", "")
    why = []
    if not best:
        why.append("no type could be computed")
    if res["support"] < min_support:
        why.append(f"support {res['support']:.0f}% < {min_support:.0f}%")
    min_s = cg.min_sites_for(region)
    if sites < min_s:
        why.append(f"{sites} aligned sites < {min_s}")
    row[f"{region}_accepted"] = "accepted" if not why else "rejected: " + "; ".join(why)
    if why:
        row["note"] = "; ".join(x for x in (row["note"], f"{region}: " + "; ".join(why)) if x)
    if region == "orf2" and best.startswith("GII.4_"):
        ident_val = res.get("nearest_identity")
        if res["support"] < min_support:
            warn = (
                "orf2: variant_reference_missing_risk (GII.4 variant unconfirmed against "
                "pinned reference manifest; inspect calicityper antigenic and phylogeny)"
            )
            row["note"] = "; ".join(x for x in (row["note"], warn) if x)
        elif ident_val not in (None, "") and float(ident_val) < 94.0:
            warn = (
                f"orf2: variant_reference_missing_risk (divergent GII.4 nearest "
                f"{res.get('nearest_type', best)} {float(ident_val):.1f}% < 94.0%; "
                f"variant may be absent from pinned reference manifest -- inspect "
                f"calicityper antigenic and phylogeny)"
            )
            row["note"] = "; ".join(x for x in (row["note"], warn) if x)


def wilson_ci(k, n, z=1.959963984540054):
    """95% Wilson score confidence interval (lo, hi) for a binomial proportion k/n."""
    if not n or n <= 0:
        return None
    p = max(0.0, min(1.0, float(k) / float(n)))
    z2 = z * z
    denom = 1.0 + z2 / n
    centre = (p + z2 / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _apply_sapovirus_result(row, region, res):
    """Sapovirus call by the published nucleotide distance cut-offs (Zhao et al. 2024).

    The nearest reference's type is the call only if the query lies within the
    genotype cut-off of it; otherwise the query is untypeable (a possible novel
    genotype if still inside the genogroup cut-off). Bootstrap support for the
    called type, a 95% Wilson CI on the observed p-distance, and the distance
    to the nearest distinct alternative type (`second_type`/`distance_margin`)
    are reported beside the call to quantify uncertainty on short fragments or
    overlapping P-type clusters. Acceptance also needs the query to cover at
    least SAPOVIRUS_MIN_REGION_COVERAGE of the region (default 0.0; see config).
    """
    gt_cut, gg_cut = cg.SAPOVIRUS_CUTOFFS[region]
    label = "RdRp" if region == "orf1" else "VP1"
    nearest_type = res.get("nearest_type") or ""
    ident = res.get("nearest_identity")
    sites = int(res.get("n_sites") or 0)
    ident_sites = int(res.get("nearest_sites") or sites)
    span = int(res.get("region_sites") or 0)
    cov = sites / span if span else 0.0
    sec_type = res.get("second_type") or ""
    sec_dist = res.get("second_distance")
    row[f"{region}_runners"] = format_runners(res.get("supports"))
    row[f"{region}_sites"] = sites
    row[f"{region}_coverage"] = round(cov, 3) if span else ""
    row[f"{region}_nearest"] = res.get("nearest", "")
    row[f"{region}_nearest_identity"] = ident if ident not in (None, "") else ""
    row[f"{region}_nearest_type"] = nearest_type
    row[f"{region}_distance_ci"] = ""
    row[f"{region}_second_type"] = sec_type
    row[f"{region}_second_distance"] = sec_dist if sec_dist is not None else ""
    row[f"{region}_distance_margin"] = ""
    why = []
    if ident in (None, "") or not nearest_type:
        row[f"{region}_type"], row[f"{region}_support"], row[f"{region}_distance"] = "unassigned", "", ""
        why.append("no reference overlaps the query enough to measure a distance")
    else:
        d = round(1.0 - float(ident) / 100.0, 4)
        row[f"{region}_distance"] = d
        ci = wilson_ci(round(d * ident_sites), ident_sites) if ident_sites > 0 else None
        if ci is not None:
            row[f"{region}_distance_ci"] = f"{ci[0]:.4f}-{ci[1]:.4f}"
        if sec_dist is not None:
            row[f"{region}_distance_margin"] = round(float(sec_dist) - d, 4)
        if d < gt_cut:
            row[f"{region}_type"] = nearest_type
            row[f"{region}_support"] = round(float((res.get("all_supports") or {}).get(nearest_type, 0.0)), 1)
            if ci is not None and ci[1] >= gt_cut:
                row["note"] = "; ".join(x for x in (
                    row["note"],
                    f"{region}: 95% Wilson CI ({ci[0]:.3f}-{ci[1]:.3f}) straddles {gt_cut} {label} cut-off"
                ) if x)
            if sec_dist is not None and float(sec_dist) < gt_cut and sec_type:
                row["note"] = "; ".join(x for x in (
                    row["note"],
                    f"{region}: second type {sec_type} (p={float(sec_dist):.3f}) also within {gt_cut} cut-off"
                ) if x)
        else:
            row[f"{region}_type"], row[f"{region}_support"] = "untypeable", ""
            if d < gg_cut:
                why.append(f"{label} p-distance {d:.3f} to the nearest reference ({nearest_type}) is over "
                           f"the {gt_cut} genotype cut-off but within the {gg_cut} genogroup cut-off: "
                           f"possible novel genotype")
            else:
                why.append(f"{label} p-distance {d:.3f} is over the {gg_cut} genogroup cut-off: "
                           f"the BLAST genogroup is not confirmed")
    if span and cov < cg.SAPOVIRUS_MIN_REGION_COVERAGE:
        why.append(f"covers {sites}/{span} {label} columns ({cov:.0%} < "
                   f"{cg.SAPOVIRUS_MIN_REGION_COVERAGE:.0%}; the cut-offs were derived on complete regions)")
    called, boot = row[f"{region}_type"], res.get("best_type") or ""
    if called not in ("untypeable", "unassigned") and boot and boot != called:
        row["note"] = "; ".join(x for x in (row["note"], f"{region}: bootstrap favours {boot} "
                                            f"over the distance call {called}") if x)
    row[f"{region}_accepted"] = "accepted" if not why else "rejected: " + "; ".join(why)
    if why:
        row["note"] = "; ".join(x for x in (row["note"], f"{region}: " + "; ".join(why)) if x)


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


def gii4_variant_window_columns(group="GII"):
    """Alignment columns of the GII.4 variant window (P2 to the VP1 3' end)."""
    import numpy as np
    from . import align
    lo_aa, hi_aa = cg.GII4_VARIANT_WINDOW_AA
    ref = align.refset("orf2", group)
    covered = np.flatnonzero(ref.occupancy >= 0.5)
    hi = len(covered) // 3 if hi_aa is None else hi_aa
    return covered[(lo_aa - 1) * 3: hi * 3]


def gii4_variant_window_coverage(aligned_query, group="GII"):
    """(covered nt, window nt, fraction) over the GII.4 variant window."""
    import numpy as np
    cols = gii4_variant_window_columns(group)
    if len(cols) == 0:
        return 0, 0, 0.0
    q = np.frombuffer(aligned_query.upper().encode(), dtype="S1")[cols]
    covered = int(np.isin(q, [b"A", b"C", b"G", b"T", b"U"]).sum())
    return covered, len(cols), round(covered / len(cols), 4)
