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
from functools import lru_cache

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


def p_distance(x, y, alphabet):
    """Uncorrected p-distance after pairwise alignment, with pairwise deletion:
    mismatches / positions where neither sequence has a gap (as MEGA does).
    Returns (p, n_compared)."""
    aln = _aligner(alphabet).align(x, y)[0]
    n = diff = 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        sa, sb = x[a0:a1], y[b0:b1]
        n += len(sa)
        diff += sum(c != d for c, d in zip(sa, sb))
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
    """Per-cluster bound table for both regions and both levels.
    Used by scripts/build_putative_novel_bounds.py; the result ships as BOUNDS_FILE."""
    rows = []
    blocks = {"orf2": [vp1_distance_matrix()], "orf1": list(rdrp_distance_blocks())}
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
                                 alphabet=REGION[region]["alphabet"], **stats[c]))
            rows.append(dict(region=region, level=level, cluster="*pooled*", n_members="",
                             alphabet=REGION[region]["alphabet"], mean=pm, sd=psd,
                             sd_source="trimmed", bound=pm + 2 * psd,
                             n_pairs=sum(len(v) for v in pairs.values())))
    return rows


@lru_cache(maxsize=None)
def bounds():
    """{(region, level, cluster): bound} from the shipped table."""
    with open(BOUNDS_FILE) as fh:
        return {(r["region"], r["level"], r["cluster"]): float(r["bound"])
                for r in csv.DictReader(fh, delimiter="\t")}


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
                distance="", closest_type="")


def _result(region, flag, t, d, b, scope=""):
    info = REGION[region]
    rel = "within" if flag == KNOWN else "outside"
    detail = (f"{info['label']} {info['alphabet']}: average p = {d:.3f} to {t}, {rel} its "
              f"2x SD bound {b:.3f}{scope}")
    if flag in FLAGGED:
        detail += "; single sequence, so tentative (NA) at most"
    return dict(flag=flag, detail=detail, distance=round(float(d), 4), closest_type=t)


def _averages(dists, labels):
    acc = {}
    for d, c in zip(dists, labels):
        if np.isfinite(d):
            acc.setdefault(c, []).append(d)
    return {c: float(np.mean(v)) for c, v in acc.items()}


def assess(region, group, qa, sub):
    """Flag one typed query region.

    qa  : the query aligned to the genogroup's reference alignment
    sub : the oriented nucleotide neighbourhood that was aligned (for VP1)
    """
    ref, mem = _members(region, group)
    if not mem:
        return _na(region, f"no {group} references")
    q = np.asarray(_encode(qa))
    if region == "orf1":
        cols = rdrp_region_columns(group)
        cover = float((q[cols] < 4).mean()) if len(cols) else 0.0
        if cover < RDRP_MIN_COVER:
            return _na(region, f"RdRp typing region {cover:.0%} covered (< {RDRP_MIN_COVER:.0%})")
        p, _ = msa_distances(q[cols], ref.codes[mem][:, cols])
        types = [genotype(ref.types[k]) for k in mem]
        group_avg = {group: float(np.nanmean(p))}
        flag, t, d, b = classify(region, _averages(p, types), group_avg)
        # other genogroups' RdRp sits in different alignments, so the P-group
        # level is judged against the query's own genogroup only
        return _result(region, flag, t, d, b, f" (P-group compared within {group} only)"
                       if flag == NOVEL_GROUP else "")

    refs = vp1_references()
    vp1 = complete_vp1(sub, [len(s) for *_, s in refs])
    if vp1 is None:
        return _na(region, "no complete VP1 open reading frame")
    # candidate genotypes: those of the nearest references by nucleotide
    # distance in the alignment; then amino acid distances to all their members
    p, _ = msa_distances(q, ref.codes[mem])
    near = [mem[k] for k in np.argsort(np.nan_to_num(p, nan=9.0))[:SHORTLIST]]
    cand = {genotype(ref.types[k]) for k in near}
    pool = [r for r in refs if genotype(r[2]) in cand]
    d = [p_distance(vp1, s, "aa")[0] for *_, s in pool]
    type_avg = _averages(d, [genotype(r[2]) for r in pool])
    flag, t, dt, b = classify(region, type_avg, None)
    if flag != KNOWN:
        # outside every candidate genotype: now compare against all of them,
        # and at genogroup level against every genogroup
        d = [p_distance(vp1, s, "aa")[0] for *_, s in refs]
        flag, t, dt, b = classify(region, _averages(d, [genotype(r[2]) for r in refs]),
                                  _averages(d, [r[1] for r in refs]))
    return _result(region, flag, t, dt, b)


SHORTLIST = 8


def _encode(s):
    from . import tree
    return tree.encode([s])[0]
