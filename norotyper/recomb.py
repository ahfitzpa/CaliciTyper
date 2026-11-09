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
    cmd = [exe, "-o", str(out_csv), "-s", str(cfg.seed),
           "-m", *list(cfg.methods), str(alignment_path)]
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
