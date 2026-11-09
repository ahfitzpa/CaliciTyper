"""Recombination detection and reporting.

Evidence-weighted annotation rather than a binary call. For each candidate
event this records which methods supported it, where the breakpoints are (in
both alignment and original query coordinates), which parental lineages are
implicated, whether segment phylogenies agree, and what was skipped and why.

The classical scanners (RDP, GENECONV, BootScan, MaxChi, Chimaera, SiScan,
3Seq) come from OpenRDP, which is run as an EXTERNAL program: it pins
numpy<2 while norotyper needs numpy>=2, so it cannot share this environment.
Install it in its own env and point $NOROTYPER_OPENRDP at the CLI. If it is
absent, the run still proceeds with whatever in-tree evidence is available and
says so in the provenance block.

Nothing here reaches "strong" on one method, and a sequence that is too short
or too ambiguous gets `insufficient_sequence` rather than a forced call.
"""
from __future__ import annotations

import csv
import io
import json
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import recomb_config as rcfg
from . import tree

VALID_CHARS = set("ATGCN-")
STATUSES = ("none", "weak", "possible", "likely", "strong", "insufficient_sequence")


class AlignmentError(ValueError):
    """The alignment cannot be used at all (as opposed to one bad query)."""


# --------------------------------------------------------------- QC and input

@dataclass
class Alignment:
    names: list
    seqs: list
    length: int

    def index(self, name):
        return self.names.index(name)


def read_alignment(path):
    """FASTA alignment -> Alignment, with ambiguity codes masked to N.

    OpenRDP accepts only A, T, G, C, N and gaps, so anything else (IUPAC R, Y,
    W ...) is masked. The count of masked characters is returned so the caller
    can record it rather than silently altering the input.
    """
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                names.append(line[1:].split()[0])
            elif line:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if not names:
        raise AlignmentError(f"no sequences found in {path}")
    if len(names) != len(seqs):
        raise AlignmentError(f"{path}: {len(names)} headers but {len(seqs)} sequences")

    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise AlignmentError(
            f"{path}: sequences have inconsistent lengths {sorted(lengths)}; "
            "the input must be an alignment, not unaligned reads")

    masked = 0
    clean = []
    for s in seqs:
        up = s.upper().replace("U", "T")
        new = re.sub(r"[^ATGCN-]", "N", up)
        masked += sum(1 for a, b in zip(up, new) if a != b)
        clean.append(new)

    dupes = [n for n in set(names) if names.count(n) > 1]
    if dupes:
        raise AlignmentError(f"{path}: duplicate sequence ids: {sorted(dupes)[:5]}")

    aln = Alignment(names, clean, lengths.pop())
    aln.masked_chars = masked
    return aln


def coordinate_map(seq):
    """Alignment column (0-based) -> 1-based ungapped position in that sequence.

    A gap column maps to the last real position seen, so a breakpoint landing
    in a gap is reported at the nearest preceding base rather than silently
    shifting. Needed because every scanner works in alignment coordinates but
    a user needs genome coordinates.
    """
    out = np.zeros(len(seq), dtype=np.int64)
    pos = 0
    for i, ch in enumerate(seq):
        if ch != "-":
            pos += 1
        out[i] = pos
    return out


def informative_sites(aln):
    """Count parsimony-informative columns: >=2 states each appearing >=2 times.

    This is the quantity that decides whether a scan can say anything at all.
    Gaps and Ns are not states. Vectorised over the four bases rather than
    looping per column, because this runs on whole-genome alignments.
    """
    arr = np.frombuffer("".join(aln.seqs).encode("ascii"), dtype=np.uint8)
    arr = arr.reshape(len(aln.seqs), aln.length)
    # counts[b, col] = how many sequences carry base b at that column
    counts = np.stack([(arr == ord(b)).sum(axis=0) for b in "ACGT"])
    return int(((counts >= 2).sum(axis=0) >= 2).sum())


def qc_query(aln, query, cfg):
    """Per-query quality gate. Returns (ok, reason, stats)."""
    i = aln.index(query)
    s = aln.seqs[i]
    ungapped = s.replace("-", "")
    n_amb = sum(1 for c in ungapped if c == "N")
    stats = dict(
        ungapped_length=len(ungapped),
        ambiguous_bases=n_amb,
        ambiguous_fraction=round(n_amb / max(1, len(ungapped)), 4),
        alignment_coverage=round(len(ungapped) / max(1, aln.length), 4),
    )
    if len(ungapped) < cfg.min_sequence_length:
        return False, (f"sequence is {len(ungapped)} nt, below min_sequence_length "
                       f"{cfg.min_sequence_length}"), stats
    if stats["ambiguous_fraction"] > cfg.max_ambiguous_fraction:
        return False, (f"{stats['ambiguous_fraction']:.1%} ambiguous bases, above "
                       f"max_ambiguous_fraction {cfg.max_ambiguous_fraction:.1%}"), stats
    return True, "", stats


# ------------------------------------------------------------------ detection

@dataclass
class Hit:
    """One significant result from one method, in alignment coordinates."""
    method: str
    start: int
    end: int
    recombinant: str
    parent1: str = ""
    parent2: str = ""
    pvalue: float = 1.0

    def midpoint(self):
        return (self.start + self.end) / 2.0

    def breakpoints(self, narrow_span):
        """Breakpoint estimates this hit supports, in alignment coordinates.

        Methods report two different things and conflating them misplaces the
        breakpoint. A window scanner (MaxChi, Chimaera, BootScan, SiScan)
        reports the WINDOW that contains the breakpoint, so the estimate is its
        centre. GENECONV and 3Seq report a converted FRAGMENT, which can span
        thousands of columns -- there the breakpoints are the fragment EDGES and
        the centre is meaningless.

        Rather than hard-code a per-method table, the span decides: a span
        narrower than `narrow_span` is treated as a window (centre), a wider one
        as a fragment (both edges).
        """
        if (self.end - self.start) <= narrow_span:
            return [self.midpoint()]
        return [float(self.start), float(self.end)]


def run_openrdp(alignment_path, cfg, workdir, timeout=7200):
    """Run OpenRDP and return (hits, used_methods, skipped, stderr).

    OpenRDP lives in its own environment (numpy<2), so it is invoked as a
    subprocess. A failure is reported, not raised: the caller still produces a
    result from whatever other evidence exists.
    """
    exe = cfg.openrdp()
    if not exe:
        return [], [], [("openrdp", "not found: set $NOROTYPER_OPENRDP or put it on PATH")], ""
    out_csv = Path(workdir) / "openrdp.csv"
    # The input path must come BEFORE -m: OpenRDP declares --methods with
    # nargs="*", so a trailing positional is swallowed as a method name and the
    # run dies with "invalid choice". Ordering it first is the whole fix.
    cmd = [exe, str(alignment_path), "-o", str(out_csv), "-s", str(cfg.seed),
           "-m", *list(cfg.methods)]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return [], [], [("openrdp", f"timed out after {timeout}s")], ""
    except OSError as e:
        return [], [], [("openrdp", f"could not execute: {e}")], ""
    if p.returncode != 0 or not out_csv.exists():
        msg = (p.stderr or p.stdout or "no output").strip().splitlines()
        return [], [], [("openrdp", f"exit {p.returncode}: {msg[-1] if msg else '?'}")], p.stderr
    hits = parse_openrdp(out_csv.read_text())
    used = sorted({h.method.lower() for h in hits})
    skipped = [(m, "ran but reported nothing significant")
               for m in cfg.methods if m not in used]
    return hits, used, skipped, p.stderr


def parse_openrdp(text):
    """Parse OpenRDP's CSV. Columns: Method,Start,End,Recombinant,Parent1,Parent2,Pvalue."""
    hits = []
    rdr = csv.DictReader(io.StringIO(text))
    for row in rdr:
        try:
            m = (row.get("Method") or "").strip()
            if not m:
                continue
            hits.append(Hit(
                method=m.lower().replace(" ", ""),
                start=int(float(row["Start"])),
                end=int(float(row["End"])),
                recombinant=(row.get("Recombinant") or "").strip(),
                parent1=(row.get("Parent1") or "").strip().lstrip("-") or "",
                parent2=(row.get("Parent2") or "").strip().lstrip("-") or "",
                pvalue=float(row["Pvalue"]),
            ))
        except (KeyError, ValueError, TypeError):
            continue                       # a malformed row must not sink the run
    return hits


def hits_for_query(hits, query):
    """Hits naming `query` as the recombinant, or as a parent.

    GENECONV reports pairwise fragments without a consistent orientation -- the
    query often appears in a parent column for an event that is really about
    it -- so both roles are collected, with the role recorded.
    """
    out = []
    for h in hits:
        if h.recombinant == query:
            out.append((h, "recombinant"))
        elif query in (h.parent1, h.parent2):
            out.append((h, "parent"))
    return out


# ----------------------------------------------------------------- clustering

def cluster_hits(hits_roles, cfg):
    """Group method hits into candidate events by breakpoint proximity.

    A significant hit from each of five methods at the same locus is one event,
    not five. Hits join a cluster when their midpoints fall within
    max_breakpoint_distance of it.
    """
    narrow = 2 * cfg.max_breakpoint_distance
    # (breakpoint estimate, hit, role) -- one hit can propose two breakpoints
    props = []
    for h, role in hits_roles:
        if h.pvalue > cfg.pvalue_threshold:
            continue
        for bp in h.breakpoints(narrow):
            props.append((bp, h, role))
    props.sort(key=lambda t: t[0])

    clusters = []
    for bp, h, role in props:
        for c in clusters:
            if abs(bp - c["midpoint"]) <= cfg.max_breakpoint_distance:
                c["hits"].append((h, role))
                c["bps"].append(bp)
                c["midpoint"] = sum(c["bps"]) / len(c["bps"])
                break
        else:
            clusters.append({"hits": [(h, role)], "bps": [bp], "midpoint": bp})

    # de-duplicate: the same hit must not count twice toward one event's
    # method support just because both its edges landed in the same cluster
    for c in clusters:
        seen, uniq = set(), []
        for h, role in c["hits"]:
            if id(h) not in seen:
                seen.add(id(h))
                uniq.append((h, role))
        c["hits"] = uniq
        c["methods"] = sorted({h.method for h, _ in c["hits"]})
    clusters.sort(key=lambda c: (-len(c["methods"]), -len(c["hits"])))
    return clusters


def correct_pvalues(hits, cfg, n_tests):
    """Bonferroni across the tests actually performed, or none."""
    if cfg.multiple_testing_correction != "bonferroni" or n_tests <= 1:
        return {id(h): h.pvalue for h in hits}
    return {id(h): min(1.0, h.pvalue * n_tests) for h in hits}


# ------------------------------------------------- segment phylogeny support

def _encode_subset(aln, cols):
    codes = tree.encode(aln.seqs)
    return codes[:, cols]


def segment_support(aln, query, breakpoint_col, cfg, rng=None):
    """Does the query's nearest relative change across the breakpoint?

    The question recombination asks phylogenetically is whether a sequence sits
    in one clade on the left of a breakpoint and a different one on the right.
    Rather than build and compare two full trees, this measures the quantity
    that matters and can be bootstrapped cheaply: the query's nearest neighbour
    in each segment, and how stable that choice is under resampling.

    Returns a dict with the two nearest neighbours, each one's bootstrap
    stability, and a 0-1 score that is zero unless the neighbour actually
    changes.
    """
    qi = aln.index(query)
    bp = int(max(0, min(aln.length - 1, breakpoint_col)))
    # Local windows, not whole segments. Comparing everything left of the
    # breakpoint against everything right of it confirms that SOME
    # recombination happened but cannot say where: for a genuine recombinant,
    # almost any split shows donor A against donor B, so a decoy breakpoint
    # scores as highly as the true one (measured: both 1.0 on the synthetic).
    # Flanking windows make the test specific to the breakpoint being assessed.
    w = cfg.segment_window if cfg.segment_window else aln.length
    left = np.arange(max(0, bp - w), bp)
    right = np.arange(bp, min(aln.length, bp + w))
    out = {"left_neighbour": "", "right_neighbour": "", "left_stability": 0.0,
           "right_stability": 0.0, "score": 0.0, "label": "none",
           "note": ""}
    if len(left) < cfg.min_segment_sites or len(right) < cfg.min_segment_sites:
        out["note"] = (f"segment shorter than min_segment_sites "
                       f"({len(left)} / {len(right)} columns)")
        return out

    rng = rng or np.random.default_rng(cfg.seed)
    others = [i for i in range(len(aln.names)) if i != qi]
    best = {}
    for side, cols in (("left", left), ("right", right)):
        codes = _encode_subset(aln, cols)
        pre = tree.indicators(codes)
        L = codes.shape[1]
        wins = defaultdict(int)
        reps = max(1, cfg.segment_bootstraps)
        for _ in range(reps):
            w = rng.multinomial(L, np.full(L, 1.0 / L))
            D = tree.tn93(codes, w, pre)
            d = D[qi, others]
            wins[others[int(np.argmin(d))]] += 1
        winner, n = max(wins.items(), key=lambda kv: kv[1])
        best[side] = (aln.names[winner], n / reps)

    out["left_neighbour"], out["left_stability"] = best["left"]
    out["right_neighbour"], out["right_stability"] = best["right"]
    changed = out["left_neighbour"] != out["right_neighbour"]
    stability = min(out["left_stability"], out["right_stability"])
    out["score"] = round(stability if changed else 0.0, 3)
    out["label"] = ("high" if out["score"] >= 0.7 else
                    "moderate" if out["score"] >= 0.5 else
                    "low" if changed else "none")
    if not changed:
        out["note"] = "same nearest neighbour on both sides: no topological change"
    return out


# ------------------------------------------------------- parental inference

def load_reference_metadata(path):
    """accession/id -> {type, ...} from TSV, CSV or JSON. Missing file is fine."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    if p.suffix.lower() == ".json":
        data = json.loads(p.read_text())
        return {str(k): (v if isinstance(v, dict) else {"type": v}) for k, v in data.items()}
    text = p.read_text()
    delim = "\t" if (p.suffix.lower() in (".tsv", ".tab") or "\t" in text.split("\n")[0]) else ","
    out = {}
    for row in csv.DictReader(io.StringIO(text), delimiter=delim):
        keys = [k for k in row if k and k.lower() in ("id", "accession", "name", "sequence_id")]
        if not keys:
            continue
        key = str(row[keys[0]]).strip()
        tkeys = [k for k in row if k and k.lower() in ("type", "genotype", "ptype", "p_type", "lineage")]
        out[key] = {"type": str(row[tkeys[0]]).strip() if tkeys else "", **row}
    return out


def _label(name, metadata):
    """Best available lineage label for a reference sequence name."""
    if name in metadata and metadata[name].get("type"):
        return metadata[name]["type"]
    base = name.split(".")[0]
    if base in metadata and metadata[base].get("type"):
        return metadata[base]["type"]
    # Without metadata, prefer a genotype-looking token ("GII.4", "GII.P16")
    # but never collapse to a bare genogroup -- "GII" as a parent label is
    # useless, so fall back to the sequence id itself.
    m = re.search(r"(G[IVX]+\.P?[A-Za-z0-9_]+)", name)
    return m.group(1) if m else name


def parental_distinctness(aln, p1, p2):
    """1 - similarity between two candidate parents, over the whole alignment.

    Near-identical "parents" are the classic false positive: any two copies of
    the same strain will look like a donor pair. This is the term that
    penalises that.
    """
    if not p1 or not p2 or p1 not in aln.names or p2 not in aln.names:
        return 0.0, "one or both candidate parents are not in the alignment"
    codes = tree.encode(aln.seqs)
    D = tree.tn93(codes)
    d = float(D[aln.index(p1), aln.index(p2)])
    # 0.10 substitutions/site is a generous "clearly different lineage" scale
    return round(min(1.0, d / 0.10), 3), ""


def infer_parents(cluster, aln, metadata, query=None):
    """Most frequently implicated parent pair across the cluster's hits.

    The query is excluded: GENECONV's orientation ambiguity means it often
    appears in a parent column for its own event, and reporting a sequence as
    its own donor is meaningless.
    """
    counts = defaultdict(int)
    for h, role in cluster["hits"]:
        for p in (h.parent1, h.parent2):
            if p and p in aln.names and p != query:
                counts[p] += 1
        if role == "parent" and h.recombinant in aln.names and h.recombinant != query:
            counts[h.recombinant] += 1          # GENECONV orientation ambiguity
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    major = ranked[0][0] if ranked else ""
    minor = ranked[1][0] if len(ranked) > 1 else ""
    amb = ""
    if len(ranked) > 2 and ranked[1][1] == ranked[2][1]:
        amb = (f"ambiguous minor parent: {ranked[1][0]} and {ranked[2][0]} are "
               f"implicated equally often")
    return major, minor, _label(major, metadata), _label(minor, metadata), amb


# ---------------------------------------------------------------- scoring

def method_support_score(methods, pvalues, cfg):
    """Independent methods agreeing, tempered by how strong each one is.

    Counting methods alone treats p=0.049 like p=1e-9. This takes the number of
    distinct methods as the dominant term and scales it by mean significance
    strength on a log scale.
    """
    n = len(set(methods))
    if n == 0:
        return 0.0
    breadth = min(1.0, n / 4.0)                   # 4+ independent methods saturates
    strengths = []
    for pv in pvalues:
        pv = max(float(pv), 1e-12)
        if pv > cfg.pvalue_threshold:
            continue
        # p at threshold -> 0, p <= 1e-6 -> 1
        strengths.append(min(1.0, np.log10(cfg.pvalue_threshold / pv) / np.log10(cfg.pvalue_threshold / 1e-6)))
    strength = float(np.mean(strengths)) if strengths else 0.0
    return round(0.7 * breadth + 0.3 * strength, 3)


def breakpoint_consistency_score(bps, cfg):
    """How tightly the methods agree on where the breakpoint is."""
    if len(bps) < 2:
        return 0.0 if len(bps) == 0 else 0.5     # a single estimate is unconfirmed
    spread = float(np.std(bps))
    return round(float(max(0.0, 1.0 - spread / max(1.0, cfg.max_breakpoint_distance))), 3)


def biological_plausibility_score(query_pos, aln_pos, aln_length, cfg):
    """Is the breakpoint somewhere recombination is actually known to happen?

    For norovirus the ORF1-ORF2 junction is the dominant hotspot. Breakpoints
    within 2% of an alignment terminus are penalised: fragment-reporting
    methods put edges there as an artefact, not as biology.
    """
    margin = max(50, int(0.02 * aln_length))
    if aln_pos <= margin or aln_pos >= aln_length - margin:
        return 0.05, "at alignment terminus: likely an edge artefact, not a breakpoint"
    best, where = 0.2, "outside the defined regions"
    for name, (lo, hi) in cfg.regions.items():
        if lo <= query_pos <= hi:
            val = 1.0 if "junction" in name.lower() else 0.6
            if val > best:
                best, where = val, name
    return round(best, 3), where


def score_event(components, cfg):
    w = cfg.scoring.weights()
    total = sum(w[k] * float(components.get(k, 0.0)) for k in w)
    return round(min(1.0, max(0.0, total)), 3)


def assign_status(n_methods, phylo_score, distinctness, cfg):
    """Map evidence onto a status label. Conservative by construction."""
    r = cfg.status_rules
    phylo_ok = phylo_score >= 0.5
    distinct_ok = distinctness >= 0.5
    if n_methods >= r.strong_min_methods and \
       (phylo_ok or not r.strong_require_phylo) and \
       (distinct_ok or not r.strong_require_distinct):
        return "strong"
    if n_methods >= r.likely_min_methods and (phylo_ok or not r.likely_require_phylo):
        return "likely"
    # One method plus a phylogeny hint is NOT "possible": on a non-recombinant
    # control that combination produced 15 spurious "possible" events. A single
    # method qualifies only with strong phylogenetic support AND distinct
    # parents; otherwise it stays "weak".
    if n_methods >= r.possible_min_methods or \
       (n_methods >= 1 and phylo_score >= 0.7 and distinct_ok):
        return "possible"
    if n_methods >= r.weak_min_methods:
        return "weak"
    return "none"


def confidence_label(score):
    return "high" if score >= 0.75 else "medium" if score >= 0.5 else "low"


# ------------------------------------------------------------- orchestration

TSV_COLUMNS = ["sequence_id", "status", "score", "confidence", "event_id",
               "breakpoint_alignment_start", "breakpoint_alignment_end",
               "breakpoint_query_start", "breakpoint_query_end",
               "supporting_methods", "candidate_major_parent",
               "candidate_minor_parent", "phylogenetic_support",
               "biological_region", "notes"]


def analyse_query(aln, query, hits, cfg, metadata, cmap=None):
    """Full recombination assessment for one query. Returns the result dict."""
    ok, reason, stats = qc_query(aln, query, cfg)
    n_inf = informative_sites(aln)
    result = {
        "sequence_id": query,
        "recombination": {"status": "none", "score": 0.0, "confidence": "low",
                          "alignment_length": aln.length,
                          "informative_sites": n_inf, "events": []},
        "qc": stats,
    }
    if not ok:
        result["recombination"]["status"] = "insufficient_sequence"
        result["recombination"]["notes"] = reason
        return result
    if n_inf < cfg.min_informative_sites:
        result["recombination"]["status"] = "insufficient_sequence"
        result["recombination"]["notes"] = (
            f"{n_inf} parsimony-informative sites, below min_informative_sites "
            f"{cfg.min_informative_sites}")
        return result

    cmap = cmap if cmap is not None else coordinate_map(aln.seqs[aln.index(query)])
    hr = hits_for_query(hits, query)
    clusters = cluster_hits(hr, cfg)

    events, best_status, best_score = [], "none", 0.0
    for n, c in enumerate(clusters, start=1):
        methods = c["methods"]
        pvals = [h.pvalue for h, _ in c["hits"]]
        bps = c["bps"]
        aln_start, aln_end = int(min(bps)), int(max(bps))
        q_start = int(cmap[min(aln_start, aln.length - 1)])
        q_end = int(cmap[min(aln_end, aln.length - 1)])

        phylo = ({"score": 0.0, "label": "not computed", "note": "build_segment_trees disabled"}
                 if not cfg.build_segment_trees
                 else segment_support(aln, query, int(c["midpoint"]), cfg))
        major, minor, major_lab, minor_lab, amb = infer_parents(c, aln, metadata, query=query)
        distinct, dnote = parental_distinctness(aln, major, minor)
        plaus, region = biological_plausibility_score(q_start, int(c["midpoint"]), aln.length, cfg)

        comps = {
            "method_support": method_support_score(methods, pvals, cfg),
            "phylogenetic_support": phylo["score"],
            "breakpoint_consistency": breakpoint_consistency_score(bps, cfg),
            "parental_distinctness": distinct,
            "biological_plausibility": plaus,
        }
        score = score_event(comps, cfg)
        status = assign_status(len(methods), phylo["score"], distinct, cfg)
        notes = "; ".join(x for x in (phylo.get("note", ""), amb, dnote) if x)

        events.append({
            "event_id": f"event_{n}",
            "breakpoints": [{"alignment_start": aln_start, "alignment_end": aln_end,
                             "query_start": q_start, "query_end": q_end}],
            "supporting_methods": methods,
            "method_pvalues": {h.method: float(h.pvalue) for h, _ in c["hits"]},
            "candidate_major_parent": major_lab,
            "candidate_minor_parent": minor_lab,
            "candidate_major_parent_id": major,
            "candidate_minor_parent_id": minor,
            "parental_distinctness": distinct,
            "phylogenetic_support": phylo["label"],
            "phylogenetic_detail": phylo,
            "breakpoint_consistency": comps["breakpoint_consistency"],
            "biological_region": region,
            "biological_plausibility": plaus,
            "score": score,
            "status": status,
            "score_components": comps,
            "notes": notes,
        })
        rank = STATUSES.index(status) if status in STATUSES else 0
        if rank > (STATUSES.index(best_status) if best_status in STATUSES else 0) or \
           (status == best_status and score > best_score):
            best_status, best_score = status, score
        best_score = max(best_score, score) if status == best_status else best_score

    events.sort(key=lambda e: -e["score"])
    rec = result["recombination"]
    rec["events"] = events
    if events:
        top = events[0]
        rec["status"], rec["score"] = top["status"], top["score"]
        rec["confidence"] = confidence_label(top["score"])
    return result


def analyse(alignment_path, queries, cfg, metadata_path=None, workdir=None,
            openrdp_timeout=7200):
    """Analyse one or more queries in an alignment. Returns (results, provenance)."""
    aln = read_alignment(alignment_path)
    metadata = load_reference_metadata(metadata_path)
    missing = [q for q in queries if q not in aln.names]
    if missing:
        raise AlignmentError(f"query id(s) not in alignment: {missing[:5]}")

    tmp = workdir or tempfile.mkdtemp(prefix="norotyper_recomb_")
    Path(tmp).mkdir(parents=True, exist_ok=True)
    clean = Path(tmp) / "alignment_clean.fasta"
    with open(clean, "w") as fh:
        for n, s in zip(aln.names, aln.seqs):
            fh.write(f">{n}\n{s}\n")

    hits, used, skipped, _ = run_openrdp(clean, cfg, tmp, timeout=openrdp_timeout)
    if cfg.run_gard:
        skipped.append(("gard", "not implemented yet (stage 3)"))

    n_tests = max(1, len(used) * max(1, len(aln.names)))
    corrected = correct_pvalues(hits, cfg, n_tests)
    for h in hits:
        h.pvalue = corrected.get(id(h), h.pvalue)

    results = []
    for q in queries:
        try:
            results.append(analyse_query(aln, q, hits, cfg, metadata))
        except Exception as e:                       # one bad query must not sink the rest
            results.append({"sequence_id": q,
                            "recombination": {"status": "none", "score": 0.0,
                                              "confidence": "low", "events": [],
                                              "notes": f"failed: {type(e).__name__}: {e}"},
                            "qc": {}})

    provenance = {
        "alignment_file": str(alignment_path),
        "alignment_sequences": len(aln.names),
        "alignment_length": aln.length,
        "ambiguity_chars_masked": getattr(aln, "masked_chars", 0),
        "reference_metadata": str(metadata_path or ""),
        "tools_used": ([f"OpenRDP ({','.join(used)})"] if used else []),
        "tools_skipped": [{"tool": t, "reason": r} for t, r in skipped],
        "parameters": {
            "pvalue_threshold": cfg.pvalue_threshold,
            "multiple_testing_correction": cfg.multiple_testing_correction,
            "n_tests_for_correction": n_tests,
            "max_breakpoint_distance": cfg.max_breakpoint_distance,
            "minimum_sequence_length": cfg.min_sequence_length,
            "min_informative_sites": cfg.min_informative_sites,
            "segment_window": cfg.segment_window,
            "segment_bootstraps": cfg.segment_bootstraps,
            "seed": cfg.seed,
            "methods_requested": list(cfg.methods),
            "scoring_weights": cfg.scoring.weights(),
        },
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return results, provenance


# -------------------------------------------------------------------- writers

def write_json(results, provenance, path):
    doc = {"tool": "norotyper recomb", "provenance": provenance, "results": results}
    Path(path).write_text(json.dumps(doc, indent=2))
    return path


def write_tsv(results, path):
    rows = []
    for r in results:
        rec = r["recombination"]
        if not rec["events"]:
            rows.append({"sequence_id": r["sequence_id"], "status": rec["status"],
                         "score": rec["score"], "confidence": rec["confidence"],
                         "notes": rec.get("notes", "")})
            continue
        for ev in rec["events"]:
            bp = ev["breakpoints"][0]
            rows.append({
                "sequence_id": r["sequence_id"], "status": ev["status"],
                "score": ev["score"], "confidence": confidence_label(ev["score"]),
                "event_id": ev["event_id"],
                "breakpoint_alignment_start": bp["alignment_start"],
                "breakpoint_alignment_end": bp["alignment_end"],
                "breakpoint_query_start": bp["query_start"],
                "breakpoint_query_end": bp["query_end"],
                "supporting_methods": ",".join(ev["supporting_methods"]),
                "candidate_major_parent": ev["candidate_major_parent"],
                "candidate_minor_parent": ev["candidate_minor_parent"],
                "phylogenetic_support": ev["phylogenetic_support"],
                "biological_region": ev["biological_region"],
                "notes": ev["notes"],
            })
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path
