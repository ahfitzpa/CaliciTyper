"""Screen calls for putative novel ORF1 (RdRp) and ORF2 (VP1) types (2x SD criterion).

Norovirus classification (Chhabra et al. 2019, J Gen Virol 100:1393) assigns
capsid genotypes and genogroups on uncorrected p-distances of the **complete
VP1 amino acid** sequence, and P-types on **partial RdRp nucleotide**
sequence. A new genotype or genogroup is recognised when it forms a separate
cluster under the 2x SD criterion: "the average distance between all sequences
within a new genogroup or genotype and its nearest cluster(s) should not
overlap within two standard deviations".

This module applies that rule as a screen, per cluster:

  each reference genotype (and genogroup) C has a bound
        bound_C = mean within-C distance + 2 x SD of within-C distances
  a query is placed by its average distance to the members of C
  it is `known` if that average is within the bound of some genotype;
  otherwise `putative_novel_type`, or `putative_novel_group` if it also
  falls outside every genogroup's bound.

Clusters with fewer than three references have at most one within-cluster
distance, so no SD of their own; they borrow the pooled SD of all clusters at
that level (far outliers, beyond Q3 + 3 IQR, dropped), and a singleton also
borrows the pooled mean. `scripts/build_putative_novel_bounds.py` writes the
per-cluster table this reads, recording which clusters borrowed.

A flagged query is one sequence, and the paper names a cluster from a single
sequence (or a single location) "not assigned" (e.g. GII.NA1). So a flag means
"tentative, worth a closer look", never a new designation. `not_assessable`
means the query lacks what the rule is defined on: a complete VP1 open reading
frame, or the RdRp typing region.

Sapovirus has published cut-offs (Zhao et al. 2024, BMC Genomics 25:1131) but
no bundled VP1/RdRp reference sets yet, so sapovirus calls are not assessed.
"""
from __future__ import annotations

import csv
import hashlib
import os
from functools import lru_cache
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

from . import config as cg

BOUNDS_FILE = cg.DATA / "putative_novel_bounds.tsv"

# What each region's criterion is defined on.
REGION = {
    "orf2": dict(label="VP1", alphabet="aa", basis="complete VP1 amino acid sequence"),
    "orf1": dict(label="RdRp", alphabet="nt", basis="partial RdRp nucleotide sequence (typing region)"),
}
# How a distance between two sequences is measured.
#   ml  patristic distance on an IQ-TREE maximum likelihood tree -- the paper's
#       measure (it uses PhyML), and the default
#   nj  patristic distance on a neighbour-joining tree, for when IQ-TREE is not
#       installed, and as a comparison
#   p   the observed pairwise distance, with no tree at all
# "patristic" is accepted as a name for "ml" so older command lines keep working.
DISTANCE = os.environ.get("CALICITYPER_NOVEL_DISTANCE", "ml").lower()
if DISTANCE == "patristic":
    DISTANCE = "ml"
MEASURES = ("ml", "nj", "p")

VP1_LEN_TOL = 0.10        # complete VP1 must be within 10% of the reference length range
RDRP_MIN_COVER = 0.95     # fraction of the RdRp typing region the query must cover
RDRP_OCCUPANCY = 0.90     # a column belongs to the typing region if >= 90% of references cover it

KNOWN, NOVEL_TYPE, NOVEL_GROUP, NA = ("known", "putative_novel_type",
                                  "putative_novel_group", "not_assessable")
FLAGGED = (NOVEL_TYPE, NOVEL_GROUP)
NOT_NOROVIRUS = ("not a norovirus: sapovirus cut-offs are published (Zhao et al. 2024) "
                 "but no sapovirus VP1/RdRp reference sets are bundled yet")


# ------------------------------------------------------------------ distances

@lru_cache(maxsize=None)
def _aligner(alphabet):
    a = PairwiseAligner(mode="global")
    if alphabet == "aa":
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -10.0, -0.5
    else:
        a.match_score, a.mismatch_score = 2.0, -3.0
        a.open_gap_score, a.extend_gap_score = -5.0, -2.0
    # a partial sequence must not pay for what it does not cover
    a.end_gap_score = 0.0
    return a


AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AMBIGUOUS_AA = frozenset("XBZJ")     # X plus Asx / Glx / Xle


def sanitise_aa(seq):
    """Map anything outside the 20 standard residues to 'X'.

    Translating an ambiguous codon yields an ambiguity letter, and BLOSUM62 has
    no row for all of them: 'MTA' (Ile or Leu) becomes 'J', which makes the
    pairwise aligner raise "sequence contains letters not in the alphabet" and
    ends the whole run. Ambiguous codons are normal in GenBank records (329 N,
    219 Y, 173 R across the 5,986-record RIVM set), not bad input.
    """
    return "".join(c if c in AA_ALPHABET else "X" for c in seq.upper())


def p_distance(x, y, alphabet):
    """Uncorrected p-distance after pairwise alignment, with pairwise deletion:
    mismatches / positions where neither sequence has a gap (as MEGA does).

    Positions where either sequence is ambiguous are skipped too, so an N-rich
    stretch lowers the number compared rather than inflating the distance --
    the same rule msa_distances() applies to nucleotides.
    """
    if alphabet == "aa":
        x, y = sanitise_aa(x), sanitise_aa(y)
    aln = _aligner(alphabet).align(x, y)[0]
    n = diff = 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for c, d in zip(x[a0:a1], y[b0:b1]):
            if c in AMBIGUOUS_AA or d in AMBIGUOUS_AA:
                continue
            n += 1
            diff += c != d
    return (diff / n if n else float("nan")), n


# ---------------------------------------------------------- reference segments

def genotype(t):
    """Genotype stem: GII.4 variants (GII.4_Sydney, ...) are below genotype level."""
    return t.split("_")[0]


def _translate(nt):
    return str(Seq(nt[: len(nt) // 3 * 3]).translate()).rstrip("*")


@lru_cache(maxsize=None)
def genogroup_of():
    """{accession: genogroup} for both regions' references."""
    out = {}
    for region in ("orf1", "orf2"):
        for r in SeqIO.parse(cg.REF_DIR[region] / "all.orf", "fasta"):
            g, _, acc = r.id.split("|")[:3]
            out[acc] = g
    return out


@lru_cache(maxsize=None)
def vp1_references():
    """[(accession, genogroup, type, VP1 amino acids)] -- complete ORF2, translated."""
    out = []
    for r in SeqIO.parse(cg.REF_DIR["orf2"] / "all.orf", "fasta"):
        g, t, acc = r.id.split("|")[:3]
        out.append((acc, g, t, _translate(str(r.seq).upper())))
    return out


@lru_cache(maxsize=None)
def _members(region, group):
    """Row indices of `group`'s own references in its (possibly shared) refset.
    Genogroups without their own alignment fall back to all.aln, which holds
    every genogroup, so the refset is filtered to the right one."""
    from . import align
    ref = align.refset(region, group)
    gg = genogroup_of()
    return ref, [k for k, n in enumerate(ref.names) if gg.get(n) == group]


def rdrp_region_columns(group):
    """Columns of the RdRp typing region: those >= RDRP_OCCUPANCY of the
    genogroup's references cover. References mix complete ORF1 with partial
    RdRp; this is the region every one of them was typed on."""
    ref, mem = _members("orf1", group)
    occ = (ref.codes[mem] < 4).mean(0) if mem else ref.occupancy
    return np.flatnonzero(occ >= RDRP_OCCUPANCY)


def as_patristic(D, method=None, records=None, alphabet="nt", threads=1):
    """Measure distances between the taxa of `D` (and `records`) under `method`.

    `D` is the observed distance matrix; `records` is [(name, aligned seq)] in
    the same order, needed only for the ML tree, which is inferred from the
    alignment rather than from distances.
    """
    import numpy as _np
    m = (method or DISTANCE)
    m = "ml" if m == "patristic" else m
    D = _np.asarray(D, dtype=float)
    if m == "p":
        return D
    if m == "nj":
        from .patristic import patristic_matrix
        return patristic_matrix(D)
    if m == "ml":
        if records is None:
            raise ValueError("the ml measure needs the alignment, not only distances")
        from .patristic import ml_patristic
        return ml_patristic(records, alphabet, threads)
    raise ValueError(f"unknown distance measure {m!r} (use one of {MEASURES})")


def augment(D, query_row):
    """Square matrix of `D` with the query appended as its last taxon."""
    n = len(D)
    out = np.zeros((n + 1, n + 1))
    out[:n, :n] = D
    q = np.nan_to_num(np.asarray(query_row, dtype=float), nan=float(np.nanmax(D)) if n else 1.0)
    out[n, :n] = out[:n, n] = q
    return out


def _cache_dir():
    """First writable cache location.

    A home directory is not always writable on a cluster -- it may be over
    quota or read-only on compute nodes -- and a cache that cannot be written
    must not take the run down with it. CALICITYPER_CACHE overrides; otherwise
    the standard cache path is tried, then the reference data directory (the
    natural home for something derived from it), then a temporary directory
    that lasts the process.
    """
    import tempfile
    env = os.environ.get("CALICITYPER_CACHE")
    candidates = [Path(env)] if env else []
    candidates += [Path(os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")) / "calicityper",
                   cg.DATA / ".cache"]
    for d in candidates:
        try:
            d.mkdir(parents=True, exist_ok=True)
            probe = d / ".writable"
            probe.touch()
            probe.unlink()
            return d
        except OSError:
            continue
    return Path(tempfile.mkdtemp(prefix="calicityper_cache_"))


@lru_cache(maxsize=1)
def vp1_reference_matrix():
    """([(accession, genogroup, genotype, aa)], observed aa distance matrix).

    The matrix is 13,695 pairwise protein alignments, about 45 seconds, so it
    is cached on disk under a key derived from the reference accessions -- a
    changed reference set gets a new key rather than a stale matrix. The cache
    lives outside the package data directory so that it cannot disturb the
    reference manifest's checksums.
    """
    refs = vp1_references()
    key = hashlib.sha256("|".join(r[0] for r in refs).encode()).hexdigest()[:16]
    path = _cache_dir() / f"vp1_aa_distances_{key}.npy"
    if path.exists():
        return refs, np.load(path)
    n = len(refs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = p_distance(refs[i][3], refs[j][3], "aa")[0]
    np.save(path, D)
    return refs, D


@lru_cache(maxsize=None)
def rdrp_alignment(group):
    """[(accession, aligned nucleotides)] over the RdRp typing region."""
    ref, mem = _members("orf1", group)
    if not mem:
        return []
    cols = rdrp_region_columns(group)
    return [(ref.names[k], "".join(np.asarray(list(ref.seqs[k]))[cols])) for k in mem]


@lru_cache(maxsize=1)
def vp1_aa_alignment():
    """([(accession, aligned amino acids)], path to the alignment file).

    One MAFFT run over the reference VP1 proteins, cached on disk beside the
    distance matrix: the ML tree needs an alignment, and re-aligning 166
    proteins per query would cost more than the tree.
    """
    refs = vp1_references()
    key = hashlib.sha256("|".join(r[0] for r in refs).encode()).hexdigest()[:16]
    path = _cache_dir() / f"vp1_aa_{key}.aln.fasta"
    if not path.exists():
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory(prefix="calicityper_aa_") as td:
            raw = Path(td) / "aa.fasta"
            raw.write_text("".join(f">{a}\n{sq}\n" for a, _, _, sq in refs))
            out = subprocess.run([cg.exe("mafft"), "--quiet", "--auto", "--anysymbol",
                                  "--thread", "2", str(raw)], capture_output=True, text=True)
            if out.returncode != 0:
                raise RuntimeError(f"mafft failed on the VP1 proteins:\n{out.stderr[-800:]}")
            path.write_text(out.stdout)
    recs = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(path, "fasta")]
    return recs, path


def add_to_vp1_alignment(aa):
    """Align one query protein to the cached VP1 reference alignment."""
    import subprocess
    import tempfile
    recs, path = vp1_aa_alignment()
    with tempfile.TemporaryDirectory(prefix="calicityper_add_") as td:
        q = Path(td) / "q.fasta"
        q.write_text(f">QUERY\n{aa}\n")
        out = subprocess.run([cg.exe("mafft"), "--quiet", "--anysymbol", "--keeplength",
                              "--add", str(q), str(path)], capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError(f"mafft failed adding the query protein:\n{out.stderr[-800:]}")
        added = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(
            __import__("io").StringIO(out.stdout), "fasta")]
    return added[:-1], added[-1][1]


@lru_cache(maxsize=None)
def rdrp_reference_matrix(group):
    """(refset, member indices, observed nucleotide distance matrix) for `group`."""
    ref, mem = _members("orf1", group)
    if not mem:
        return ref, mem, np.zeros((0, 0))
    cols = rdrp_region_columns(group)
    C = ref.codes[mem][:, cols]
    D = np.vstack([msa_distances(C[a], C)[0] for a in range(len(mem))])
    D = np.nan_to_num(D, nan=0.0)
    np.fill_diagonal(D, 0.0)
    return ref, mem, D


def msa_distances(qcodes, rcodes):
    """p-distance of one encoded aligned query to each encoded aligned reference,
    with pairwise deletion. Returns (p, n) arrays."""
    both = (qcodes < 4)[None, :] & (rcodes < 4)
    n = both.sum(1)
    diff = ((rcodes != qcodes[None, :]) & both).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, diff / np.maximum(n, 1), np.nan), n


# --------------------------------------------------------------- query segments

def complete_vp1(nt, ref_lengths):
    """The query's complete VP1 as amino acids, or None.

    Complete means an ATG...stop open reading frame, no internal stop, whose
    length is within VP1_LEN_TOL of the reference VP1 length range. `nt` is the
    oriented ORF2 neighbourhood the typing step aligned.
    """
    s = nt.upper().replace("U", "T")
    lo = min(ref_lengths) * (1 - VP1_LEN_TOL)
    hi = max(ref_lengths) * (1 + VP1_LEN_TOL)
    best = None
    for frame in range(3):
        aa = str(Seq(s[frame: frame + (len(s) - frame) // 3 * 3]).translate())
        start = 0
        while True:
            m = aa.find("M", start)
            if m < 0:
                break
            stop = aa.find("*", m)
            if stop < 0:
                break                              # runs off the end: truncated, not complete
            orf = aa[m:stop]
            if lo <= len(orf) <= hi and (best is None or len(orf) > len(best)):
                best = orf
            start = stop + 1
    return best


# ----------------------------------------------------------------------- bounds

def _trimmed(v):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    if len(v) < 4:
        return v
    q1, q3 = np.quantile(v, [0.25, 0.75])
    return v[v <= q3 + 3 * (q3 - q1)]


MIN_OWN = 3              # references a cluster needs for its own mean and SD


def is_tentative(cluster):
    """Not-assigned clusters (GNA1, GII.NA1, GII.PNA4, ...): the 2x SD criterion
    cannot be applied to them, so they never supply their own bound."""
    return "NA" in cluster.split(".")[-1] or cluster.startswith("GNA")


def cluster_stats(pairs, members, own_level=True):
    """Bounds per cluster.

    pairs   : {cluster: [distances between its members]}
    members : {cluster: number of references}
    own_level : False puts every cluster on the pooled bound (used for
                genogroups, whose internal diversity varies too much between
                so few groups for one diverse group not to swallow the rest)

    A cluster uses its own mean + 2 SD only when it has >= MIN_OWN references
    and is not tentative; otherwise the pooled (outlier-trimmed) distribution.
    Returns ({cluster: dict(mean, sd, sd_source, bound, n_pairs)}, (pooled mean, sd)).
    """
    pooled = _trimmed([d for v in pairs.values() for d in v])
    pm, psd = float(pooled.mean()), float(pooled.std(ddof=1))
    out = {}
    for c, v in pairs.items():
        v = np.asarray(v, float)
        v = v[~np.isnan(v)]
        if own_level and members.get(c, 0) >= MIN_OWN and len(v) >= 2 and not is_tentative(c):
            m, sd, src = float(v.mean()), float(v.std(ddof=1)), "own"
        else:
            m, sd = pm, psd
            src = ("pooled_tentative" if is_tentative(c) else
                   "pooled_small" if own_level else "pooled_level")
        out[c] = dict(mean=m, sd=sd, sd_source=src, bound=m + 2 * sd, n_pairs=len(v))
    return out, (pm, psd)


OWN_LEVEL = {"genotype": True, "genogroup": False}


def vp1_distance_matrix():
    refs = vp1_references()
    n = len(refs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = p_distance(refs[i][3], refs[j][3], "aa")[0]
    return [(a, g, genotype(t)) for a, g, t, _ in refs], D


def rdrp_distance_blocks():
    """RdRp distances exist only within one genogroup's alignment."""
    for g in sorted(set(genogroup_of().values())):
        ref, mem = _members("orf1", g)
        if not mem:
            continue
        cols = rdrp_region_columns(g)
        C = ref.codes[mem][:, cols]
        D = np.vstack([msa_distances(C[a], C)[0] for a in range(len(mem))])
        yield [(ref.names[k], g, genotype(ref.types[k])) for k in mem], D


def within_pairs(meta, D, level, skip=None):
    """{cluster: [distances between its members]} at `level` (genotype/genogroup),
    optionally leaving one reference out (for validation)."""
    key = (lambda m: m[2]) if level == "genotype" else (lambda m: m[1])
    out = {}
    for i, mi in enumerate(meta):
        if i == skip:
            continue
        out.setdefault(key(mi), [])
        for j in range(i + 1, len(meta)):
            if j != skip and key(meta[j]) == key(mi):
                out[key(mi)].append(D[i, j])
    return out


def build_bounds():
    """Per-cluster bound table, for both regions, both levels and both distance
    measures. Used by scripts/build_putative_novel_bounds.py; ships as BOUNDS_FILE.

    A bound is only meaningful in the units it was computed in, so the table
    carries a `distance` column and holds a set of rows per measure. Otherwise
    running with CALICITYPER_NOVEL_DISTANCE=p would compare observed distances
    against patristic bounds, which is how an apparently small change of
    setting turns into a wrong answer."""
    raw = {"orf2": [vp1_distance_matrix()], "orf1": list(rdrp_distance_blocks())}
    # the alignment each block was measured on, in the same row order, for the
    # ML tree; VP1's reference alignment is re-ordered to match the matrix
    aa_by_acc = dict(vp1_aa_alignment()[0])
    alns = {"orf2": [[(m[0], aa_by_acc[m[0]]) for m in raw["orf2"][0][0]]],
            "orf1": [rdrp_alignment(meta[0][1]) for meta, _ in raw["orf1"]]}
    rows = []
    for measure in MEASURES:
        rows += _bounds_for(raw, alns, measure)
    return rows


def _bounds_for(raw, alns, measure):
    rows = []
    # the bound and the query's distance must be in the same units, so the
    # transform applied here is the one assess() applies
    blocks = {}
    for region, blks in raw.items():
        alphabet = REGION[region]["alphabet"]
        blocks[region] = [
            (meta, as_patristic(D, measure, records=alns[region][b], alphabet=alphabet))
            for b, (meta, D) in enumerate(blks)]
    for region, blks in blocks.items():
        for level in ("genotype", "genogroup"):
            pairs, size = {}, {}
            for meta, D in blks:
                for c, v in within_pairs(meta, D, level).items():
                    pairs.setdefault(c, []).extend(v)
                for m in meta:
                    c = m[2] if level == "genotype" else m[1]
                    size[c] = size.get(c, 0) + 1
            stats, (pm, psd) = cluster_stats(pairs, size, OWN_LEVEL[level])
            for c in sorted(stats):
                rows.append(dict(region=region, level=level, cluster=c, n_members=size[c],
                                 alphabet=REGION[region]["alphabet"], distance=measure,
                                 **stats[c]))
            rows.append(dict(region=region, level=level, cluster="*pooled*", n_members="",
                             alphabet=REGION[region]["alphabet"], distance=measure,
                             mean=pm, sd=psd, sd_source="trimmed", bound=pm + 2 * psd,
                             n_pairs=sum(len(v) for v in pairs.values())))
    return rows


@lru_cache(maxsize=None)
def bounds(measure=None):
    """{(region, level, cluster): bound} from the shipped table, for one measure.

    A table written before the `distance` column existed holds observed-distance
    bounds, so it is read as such rather than silently treated as patristic."""
    m = measure or DISTANCE
    with open(BOUNDS_FILE) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if rows and "distance" not in rows[0]:
        if m != "p":
            raise RuntimeError(
                f"{BOUNDS_FILE} predates patristic distances and holds observed-distance "
                f"bounds only. Rebuild it with scripts/build_putative_novel_bounds.py, or "
                f"set CALICITYPER_NOVEL_DISTANCE=p to use it as it stands.")
        rows = [dict(r, distance="p") for r in rows]
    out = {(r["region"], r["level"], r["cluster"]): float(r["bound"])
           for r in rows if r["distance"] == m}
    if not out:
        raise RuntimeError(f"no {m!r} bounds in {BOUNDS_FILE}; rebuild it with "
                           f"scripts/build_putative_novel_bounds.py")
    return out


def bound_of(region, level, cluster, table=None):
    t = table or bounds()
    return t.get((region, level, cluster), t[(region, level, "*pooled*")])


# ------------------------------------------------------------------------ screen

def classify(region, type_avg, group_avg, table=None):
    """Apply the rule to a query's average distances.

    type_avg  : {genotype: average distance to its members}
    group_avg : {genogroup: average distance to its members}, or None when the
                query already sits inside a known genotype (so inside its
                genogroup) and the genogroup level need not be computed
    Returns (flag, nearest genotype, its average, its bound).
    """
    t_near = min(type_avg, key=type_avg.get)
    inside = [t for t, d in type_avg.items() if d <= bound_of(region, "genotype", t, table)]
    if inside:
        best = min(inside, key=type_avg.get)
        return KNOWN, best, type_avg[best], bound_of(region, "genotype", best, table)
    flag = NOVEL_TYPE
    if group_avg is not None and not any(
            d <= bound_of(region, "genogroup", g, table) for g, d in group_avg.items()):
        flag = NOVEL_GROUP
    return flag, t_near, type_avg[t_near], bound_of(region, "genotype", t_near, table)


def _na(region, why):
    return dict(flag=NA, detail=f"{why}; the criterion is defined on the {REGION[region]['basis']}",
                distance="", closest_type="", bound="")


def _result(region, flag, t, d, b, scope=""):
    info = REGION[region]
    rel = "within" if flag == KNOWN else "outside"
    unit = "patristic" if DISTANCE == "patristic" else "p"
    detail = (f"{info['label']} {info['alphabet']}: average {unit} distance = {d:.3f} to {t}, "
              f"{rel} its 2x SD bound {b:.3f}{scope}")
    if flag in FLAGGED:
        detail += "; single sequence, so tentative (NA) at most"
    return dict(flag=flag, detail=detail, distance=round(float(d), 4), closest_type=t,
                bound=round(float(b), 4))


def _averages(dists, labels):
    acc = {}
    for d, c in zip(dists, labels):
        if np.isfinite(d):
            acc.setdefault(c, []).append(d)
    return {c: float(np.mean(v)) for c, v in acc.items()}


def assess(region, group, qa, sub):
    """Flag one typed query region against the 2x SD criterion.

    qa  : the query aligned to the genogroup's reference alignment
    sub : the oriented nucleotide neighbourhood that was aligned (for VP1)

    The query is appended to the reference distance matrix and the whole thing
    is turned into patristic distances in one tree, because a patristic
    distance is a property of a tree and not of a pair: reading the query off
    one tree and the bounds off another would compare two different quantities.
    """
    if region == "orf1":
        ref, mem, Dref = rdrp_reference_matrix(group)
        if not mem:
            return _na(region, f"no {group} references")
        q = np.asarray(_encode(qa))
        cols = rdrp_region_columns(group)
        cover = float((q[cols] < 4).mean()) if len(cols) else 0.0
        if cover < RDRP_MIN_COVER:
            return _na(region, f"RdRp typing region {cover:.0%} covered (< {RDRP_MIN_COVER:.0%})")
        obs, _ = msa_distances(q[cols], ref.codes[mem][:, cols])
        recs = rdrp_alignment(group) + [("QUERY", "".join(np.asarray(list(qa))[cols]))]
        P = as_patristic(augment(Dref, obs), records=recs, alphabet="nt")
        d = P[-1, :-1]
        types = [genotype(ref.types[k]) for k in mem]
        group_avg = {group: float(np.nanmean(d))}
        flag, t, dist, b = classify(region, _averages(d, types), group_avg)
        # other genogroups' RdRp sits in different alignments, so the P-group
        # level is judged against the query's own genogroup only
        return _result(region, flag, t, dist, b, f" (P-group compared within {group} only)"
                       if flag == NOVEL_GROUP else "")

    refs, Dref = vp1_reference_matrix()
    if not refs:
        return _na(region, "no VP1 references")
    vp1 = complete_vp1(sub, [len(s) for *_, s in refs])
    if vp1 is None:
        return _na(region, "no complete VP1 open reading frame")
    obs = [p_distance(vp1, s, "aa")[0] for *_, s in refs]
    recs = None
    if DISTANCE == "ml":
        aligned_refs, aligned_q = add_to_vp1_alignment(vp1)
        recs = aligned_refs + [("QUERY", aligned_q)]
    P = as_patristic(augment(Dref, obs), records=recs, alphabet="aa")
    d = P[-1, :-1]
    # both levels come off the one tree: the shortlist the observed-distance
    # version needed to keep the cost down is unnecessary here, because the
    # matrix is cached and every reference is already in the tree
    flag, t, dist, b = classify(region, _averages(d, [genotype(r[2]) for r in refs]),
                                _averages(d, [r[1] for r in refs]))
    return _result(region, flag, t, dist, b)


SHORTLIST = 8


def _encode(s):
    from . import tree
    return tree.encode([s])[0]
