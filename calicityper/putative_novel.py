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

Sapovirus is not screened here. Its VP1 and RdRp reference sets are bundled and
used for typing, but the published RdRp cut-off (Zhao et al. 2024, BMC Genomics
25:1131) overlaps the body of the between-type distance distribution: only 1 of 20
testable P-types is separated from every other P-type in its genogroup by the
cut-off, and none of five singleton P-types was reported untypeable when its only
reference was removed. The cut-off therefore cannot establish that a query belongs
to none of the known P-types, which is what a novelty screen must do. VP1 passes
that test 6 of 6, but a VP1-only sapovirus screen would be a different rule from
the dual one the paper defines, so none is offered. See docs/SAPOVIRUS_STATUS.md.
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
# Opt-in constrained IQ-TREE (-g) placement against a cached reference backbone tree (S6).
FAST_NOVEL = os.environ.get("CALICITYPER_FAST_NOVEL", "").lower() in ("1", "true", "yes")

VP1_LEN_TOL = 0.10        # complete VP1 must be within 10% of the reference length range
RDRP_MIN_COVER = 0.95     # fraction of the RdRp typing region the query must cover
RDRP_OCCUPANCY = 0.90     # a column belongs to the typing region if >= 90% of references cover it

KNOWN, NOVEL_TYPE, NOVEL_GROUP, NA = ("known", "putative_novel_type",
                                  "putative_novel_group", "not_assessable")
FLAGGED = (NOVEL_TYPE, NOVEL_GROUP)
NOT_NOROVIRUS = ("not a norovirus: this screen is norovirus only. Sapovirus VP1 and "
                 "RdRp reference sets are bundled and used for typing, but the "
                 "published RdRp cut-off overlaps the between-type distance "
                 "distribution and cannot establish that a query belongs to none of "
                 "the known P-types, so no sapovirus novelty screen is offered. See "
                 "docs/SAPOVIRUS_STATUS.md")
SAPOVIRUS_RDRP_NA = (
    "sapovirus RdRp is not screened: the published 0.266 P-type cut-off overlaps the "
    "between-type distance distribution -- only 1 of 20 testable P-types is separated "
    "from every other P-type in its genogroup -- so a distance beyond the cut-off "
    "cannot be read as novelty. See docs/SAPOVIRUS_STATUS.md")


def sapovirus_vp1_flag(row):
    """(flag, detail) for a sapovirus VP1 call, from the published cut-offs.

    Deliberately not a 2x SD bound like the norovirus screen. One was built and
    measured on the bundled sapovirus VP1 references (2026-11-20): on
    nucleotide, the units the published scheme is defined in, the genotype-level
    bound came out at a median of 0.262, above the 0.161 cut-off for 37 of 50
    clusters, so it would have flagged strictly less than the cut-off the typing
    path already applies. (Cross-genogroup nucleotide p-distance also proved
    unmeasurable by pairwise alignment, which collapsed the tree the bound
    needs.) So the signal reported here IS the cut-off outcome: a VP1 further
    than the genotype cut-off from every reference is a possible novel genotype,
    which is what the typing path records as `untypeable`.

    Unlike the norovirus screen this adds no computation -- it reads the call
    the typing step already made -- and it is still screening, not a
    designation: the VP1 reference set is unevenly sampled, and a flag means
    "no represented genotype is within the cut-off", not "new genotype".
    """
    gt_cut, gg_cut = cg.SAPOVIRUS_CUTOFFS["orf2"]
    called = str(row.get("orf2_type") or "").strip()
    try:
        d = float(row.get("orf2_distance"))
    except (TypeError, ValueError):
        d = None
    if called == "unassigned" or not called:
        return NA, "no sapovirus VP1 reference overlaps the query enough to measure a distance"
    if called != "untypeable":
        return KNOWN, (f"VP1 p-distance {d:.3f} to the nearest reference is within the "
                       f"{gt_cut} genotype cut-off" if d is not None else
                       f"VP1 within the {gt_cut} genotype cut-off of {called}")
    if d is None:
        return NA, "VP1 is untypeable but no distance was recorded"
    if d < gg_cut:
        return NOVEL_TYPE, (f"VP1 p-distance {d:.3f} to the nearest reference is over the "
                            f"{gt_cut} genotype cut-off but within the {gg_cut} genogroup "
                            f"cut-off: possible novel genotype within a known genogroup")
    return NOVEL_GROUP, (f"VP1 p-distance {d:.3f} is over the {gg_cut} genogroup cut-off: "
                         f"the query sits outside every bundled sapovirus genogroup")


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


_WARNED = set()


def _warn_once(msg):
    if msg not in _WARNED:
        _WARNED.add(msg)
        print(f"calicityper: {msg}", file=__import__("sys").stderr)


def effective_measure(measure=None):
    """The measure that will actually be used, after any fallback."""
    m = "ml" if (measure or DISTANCE) == "patristic" else (measure or DISTANCE)
    if m != "ml":
        return m
    try:
        cg.exe("iqtree")
    except FileNotFoundError:
        return "nj"
    return m


def as_patristic(D, method=None, records=None, alphabet="nt", threads=1, constraint=None):
    """Measure distances between the taxa of `D` (and `records`) under `method`.

    `D` is the observed distance matrix; `records` is [(name, aligned seq)] in
    the same order, needed only for the ML tree, which is inferred from the
    alignment rather than from distances.
    """
    import numpy as _np
    m = (method or DISTANCE)
    m = "ml" if m == "patristic" else m
    if m == "p":
        return _np.asarray(D, dtype=float)
    if m == "nj":
        from .patristic import patristic_matrix
        return patristic_matrix(_np.asarray(D, dtype=float))
    if m == "ml":
        if records is None:
            raise ValueError("the ml measure needs the alignment, not only distances")
        from .patristic import ml_patristic
        try:
            return ml_patristic(records, alphabet, threads, constraint=constraint)
        except FileNotFoundError as e:
            # IQ-TREE is a newer dependency than the rest of the tool. Falling
            # back to the NJ tree is honest -- the measure is reported with the
            # result and its own bounds are used -- and better than refusing to
            # screen at all, but it is said out loud, because the ML tree is
            # what lifts the RdRp and genogroup screens.
            _warn_once(f"{e}\nFalling back to neighbour joining for the 2x SD screen "
                       f"(CALICITYPER_NOVEL_DISTANCE=nj). Sensitivity is lower: RdRp P-type "
                       f"50% rather than 63%, genogroup 21% rather than 52%.")
            if D is None:
                raise
            from .patristic import patristic_matrix
            return patristic_matrix(_np.asarray(D, dtype=float))
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
    m = effective_measure(measure)
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
    # The measure has to be named, and named after any IQ-TREE fallback: a
    # patristic distance off an ML tree and a raw p distance are different
    # quantities on different scales, and the bound is only read in its own
    # units. Saying "p" for all three (as this did) hid which one was used.
    unit = {"ml": "ML patristic", "nj": "NJ patristic", "p": "p"}[effective_measure()]
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


def _cached_backbone(tag, records, alphabet, threads=1):
    """Return path to a cached reference-only ML tree for IQ-TREE `-g`, or None."""
    if not FAST_NOVEL or effective_measure() != "ml" or len(records) < 4:
        return None
    key = hashlib.sha256("|".join(f"{n}:{len(s)}" for n, s in records).encode()).hexdigest()[:16]
    path = _cache_dir() / f"backbone_{tag}_{alphabet}_{key}.treefile"
    if path.exists() and path.stat().st_size > 0:
        return path
    from .patristic import backbone_newick
    try:
        nwk = backbone_newick(records, alphabet=alphabet, threads=threads)
    except FileNotFoundError:
        return None
    if not nwk:
        return None
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        tmp.write_text(nwk + "\n")
        os.replace(tmp, path)
    except OSError:
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass
        return nwk
    return path


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
        ref_recs = rdrp_alignment(group)
        recs = ref_recs + [("QUERY", "".join(np.asarray(list(qa))[cols]))]
        constraint = _cached_backbone(f"orf1_{group}", ref_recs, "nt")
        P = as_patristic(augment(Dref, obs), records=recs, alphabet="nt", constraint=constraint)
        d = P[-1, :-1]
        types = [genotype(ref.types[k]) for k in mem]
        group_avg = {group: float(np.nanmean(d))}
        flag, t, dist, b = classify(region, _averages(d, types), group_avg)
        # other genogroups' RdRp sits in different alignments, so the P-group
        # level is judged against the query's own genogroup only
        return _result(region, flag, t, dist, b, f" (P-group compared within {group} only)"
                       if flag == NOVEL_GROUP else "")

    refs = vp1_references()
    if not refs:
        return _na(region, "no VP1 references")
    vp1 = complete_vp1(sub, [len(s) for *_, s in refs])
    if vp1 is None:
        return _na(region, "no complete VP1 open reading frame")
    measure = effective_measure()
    if measure == "ml":
        aligned_refs, aligned_q = add_to_vp1_alignment(vp1)
        aa_by_acc = dict(aligned_refs)
        ordered_refs = [(r[0], aa_by_acc[r[0]]) for r in refs]
        recs = ordered_refs + [("QUERY", aligned_q)]
        constraint = _cached_backbone("orf2_vp1", ordered_refs, "aa")
        P = as_patristic(None, records=recs, alphabet="aa", constraint=constraint)
    else:
        refs, Dref = vp1_reference_matrix()
        obs = [p_distance(vp1, s, "aa")[0] for *_, s in refs]
        P = as_patristic(augment(Dref, obs), method=measure, records=None, alphabet="aa")
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
