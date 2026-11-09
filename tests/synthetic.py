"""Synthetic alignments with known recombination history, for testing.

Real norovirus reference genomes cannot serve as negative controls: norovirus
recombines constantly around the ORF1-ORF2 junction, so a published reference
strain IS frequently a recombinant and a "possible" call on one may be correct.
Specificity therefore has to be measured against sequences that are
non-recombinant by construction, which is what these generators provide.

Every generator is seeded, so a test failure is reproducible.
"""
from __future__ import annotations

import numpy as np

BASES = "ACGT"
ORF1_ORF2_JUNCTION = 5100          # approximate junction on a GII-like genome


def _random_genome(rng, length):
    return "".join(rng.choice(list(BASES), size=length))


def _mutate(seq, rng, rate):
    arr = np.array(list(seq))
    n = len(arr)
    idx = rng.choice(n, size=max(1, int(n * rate)), replace=False)
    for i in idx:
        arr[i] = rng.choice([b for b in BASES if b != arr[i]])
    return "".join(arr)


def _panel(rng, n_lineages, length, divergence):
    """n_lineages distinct lineages, each with two near-identical members."""
    root = _random_genome(rng, length)
    lineages = {}
    for k in range(n_lineages):
        anc = _mutate(root, rng, divergence)
        lineages[f"lineage{k}"] = [_mutate(anc, rng, divergence / 10) for _ in range(2)]
    return lineages


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    return str(path)


def _records(lineages):
    out = []
    for lname, seqs in lineages.items():
        for i, s in enumerate(seqs):
            out.append((f"{lname}_ref{i}", s))
    return out


def non_recombinant(path, seed=1, length=7500):
    """A query descended from ONE lineage. No recombination exists by construction."""
    rng = np.random.default_rng(seed)
    lin = _panel(rng, 4, length, 0.08)
    q = _mutate(lin["lineage0"][0], rng, 0.01)
    recs = _records(lin) + [("QUERY", q)]
    return write_fasta(path, recs), {"truth": "none", "breakpoint": None}


def clear_recombinant(path, seed=2, length=7500, breakpoint=ORF1_ORF2_JUNCTION):
    """Query is lineage0 left of the breakpoint and lineage2 right of it."""
    rng = np.random.default_rng(seed)
    lin = _panel(rng, 4, length, 0.08)
    a, b = lin["lineage0"][0], lin["lineage2"][0]
    q = _mutate(a[:breakpoint] + b[breakpoint:], rng, 0.002)
    recs = _records(lin) + [("QUERY", q)]
    return write_fasta(path, recs), {"truth": "recombinant", "breakpoint": breakpoint}


def weak_signal(path, seed=3, length=7500, tract=150):
    """A very short imported tract: at most marginal, must not reach 'likely'."""
    rng = np.random.default_rng(seed)
    lin = _panel(rng, 4, length, 0.08)
    a, b = lin["lineage0"][0], lin["lineage3"][0]
    s = ORF1_ORF2_JUNCTION
    q = _mutate(a[:s] + b[s:s + tract] + a[s + tract:], rng, 0.002)
    recs = _records(lin) + [("QUERY", q)]
    return write_fasta(path, recs), {"truth": "weak", "breakpoint": s}


def conflicting(path, seed=4, length=7500):
    """Two imported tracts from different donors, far apart: conflicting evidence."""
    rng = np.random.default_rng(seed)
    lin = _panel(rng, 4, length, 0.08)
    a, b, c = lin["lineage0"][0], lin["lineage1"][0], lin["lineage2"][0]
    q = a[:1500] + b[1500:2500] + a[2500:5200] + c[5200:6200] + a[6200:]
    q = _mutate(q, rng, 0.002)
    recs = _records(lin) + [("QUERY", q)]
    return write_fasta(path, recs), {"truth": "multiple", "breakpoints": [1500, 2500, 5200, 6200]}


def orf_mismatch_only(path, seed=5, length=7500):
    """Distinct ORF1 and ORF2 ancestry but with donors too similar to resolve.

    Models the case the spec warns about: an ORF1/ORF2 type mismatch is a clue,
    not proof. Parental distinctness should be low and the status must not be
    'strong'.
    """
    rng = np.random.default_rng(seed)
    root = _random_genome(rng, length)
    near_a = _mutate(root, rng, 0.004)          # two nearly identical donors
    near_b = _mutate(root, rng, 0.004)
    others = [_mutate(root, rng, 0.08) for _ in range(4)]
    s = ORF1_ORF2_JUNCTION
    q = _mutate(near_a[:s] + near_b[s:], rng, 0.001)
    recs = [("nearA_ref0", near_a), ("nearA_ref1", _mutate(near_a, rng, 0.001)),
            ("nearB_ref0", near_b), ("nearB_ref1", _mutate(near_b, rng, 0.001))]
    recs += [(f"other_ref{i}", s_) for i, s_ in enumerate(others)]
    recs += [("QUERY", q)]
    return write_fasta(path, recs), {"truth": "ambiguous_parents", "breakpoint": s}


def too_short(path, seed=6, length=900):
    """Below min_sequence_length: must return insufficient_sequence."""
    rng = np.random.default_rng(seed)
    lin = _panel(rng, 3, length, 0.08)
    a, b = lin["lineage0"][0], lin["lineage2"][0]
    q = a[:400] + b[400:]
    recs = _records(lin) + [("QUERY", q)]
    return write_fasta(path, recs), {"truth": "insufficient_sequence", "breakpoint": 400}


ALL = {
    "non_recombinant": non_recombinant,
    "clear_recombinant": clear_recombinant,
    "weak_signal": weak_signal,
    "conflicting": conflicting,
    "orf_mismatch_only": orf_mismatch_only,
    "too_short": too_short,
}
