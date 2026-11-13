"""Patristic distances from a neighbour-joining tree.

Chhabra et al. 2019 apply the 2x SD criterion to *patristic* distances -- the
sum of branch lengths along the path between two tips of a phylogeny -- not to
the observed pairwise distances this package measured before. The two differ in
a way that matters for a clustering criterion: an observed distance saturates
as substitutions accumulate at the same sites, while branch lengths fitted to a
tree add up along a path, so distant clusters separate further and the spread
*within* a tight cluster is re-expressed in the same units.

The tree is inferred by maximum likelihood with IQ-TREE, as the paper infers
its trees by maximum likelihood (with PhyML). `-fast` is used, which is a
single tree search rather than the full multi-start default: on the bundled
sets that is about a second for a genogroup's RdRp and six for VP1, so the
screen can afford one tree per query, which is what putting the query and the
references in the same tree requires. Models are fixed -- GTR+G4 for
nucleotides, LG+G4 for amino acids -- rather than selected per run, so two runs
of the same data give the same bound.

Neighbour joining is kept as a fallback for when IQ-TREE is not installed, and
as a comparison: it is a least-squares fit to additivity rather than a
likelihood optimum, so its branch lengths are close but not identical.

Patristic distance between two tips depends on which taxa are in the tree, so
the bounds and the query must be read off ONE tree. `patristic_matrix` is
therefore given the whole matrix -- references plus the query -- at once, and
the bound table is rebuilt from the same function.
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

MODEL = {"nt": "GTR+G4", "aa": "LG+G4"}


def nj_edges(D):
    """Neighbour-joining on a square distance matrix.

    Returns (edges, n_nodes) where edges is a list of (a, b, length) over node
    indices: tips are 0..n-1, internal nodes n.. . Negative branch lengths,
    which NJ can produce when the data are not additive, are clamped to zero --
    a negative length would make a patristic distance shorter than the path it
    contains.
    """
    D = np.array(D, dtype=float, copy=True)
    n = D.shape[0]
    if n < 2:
        return [], n
    np.fill_diagonal(D, 0.0)
    active = list(range(n))
    nxt = n
    edges = []
    # work on a growing matrix indexed by node id
    size = 2 * n
    M = np.full((size, size), np.nan)
    M[:n, :n] = D
    while len(active) > 2:
        idx = np.array(active)
        sub = M[np.ix_(idx, idx)]
        r = len(active)
        totals = sub.sum(axis=1)
        # Q = (r - 2) d_ij - sum_i - sum_j, minimised off the diagonal
        Q = (r - 2) * sub - totals[:, None] - totals[None, :]
        np.fill_diagonal(Q, np.inf)
        a, b = np.unravel_index(np.argmin(Q), Q.shape)
        i, j = idx[a], idx[b]
        dij = M[i, j]
        li = 0.5 * dij + (totals[a] - totals[b]) / (2 * (r - 2))
        lj = dij - li
        edges.append((i, nxt, max(0.0, li)))
        edges.append((j, nxt, max(0.0, lj)))
        for k in active:
            if k in (i, j):
                continue
            M[nxt, k] = M[k, nxt] = 0.5 * (M[i, k] + M[j, k] - dij)
        M[nxt, nxt] = 0.0
        active = [k for k in active if k not in (i, j)] + [nxt]
        nxt += 1
    if len(active) == 2:
        i, j = active
        edges.append((i, j, max(0.0, M[i, j])))
    return edges, nxt


def patristic_matrix(D):
    """All-pairs patristic distances between the tips of the NJ tree of `D`."""
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    if n < 2:
        return np.zeros((n, n))
    edges, n_nodes = nj_edges(D)
    adj = [[] for _ in range(n_nodes)]
    for a, b, w in edges:
        adj[a].append((b, w))
        adj[b].append((a, w))
    out = np.zeros((n, n))
    for t in range(n):
        # the graph is a tree, so one traversal gives every distance from t
        dist = np.full(n_nodes, np.nan)
        dist[t] = 0.0
        stack = [t]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if np.isnan(dist[v]):
                    dist[v] = dist[u] + w
                    stack.append(v)
        out[t] = dist[:n]
    np.fill_diagonal(out, 0.0)
    return out


# --------------------------------------------------------------- ML trees

def _edges_from_newick(newick, n_tips):
    """Parse a newick string into (edges, n_nodes) over integer node ids.

    Tips are named t0..t{n-1} by the caller, so a tip's id is its number and
    internal nodes are numbered after them. Written directly rather than via
    Bio.Phylo because all that is needed is the edge list the traversal uses.
    """
    edges = []
    nxt = [n_tips]
    pos = [0]
    s = newick.strip().rstrip(";")

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            me = nxt[0]
            nxt[0] += 1
            while True:
                child, length = node()
                edges.append((child, me, length))
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                break
            assert s[pos[0]] == ")", s[pos[0]]
            pos[0] += 1
        else:
            m = re.match(r"[^,():]+", s[pos[0]:])
            label = m.group(0)
            pos[0] += len(label)
            me = int(label[1:]) if label.startswith("t") else int(label)
        # an optional support label then an optional branch length
        m = re.match(r"[^,():]*", s[pos[0]:])
        if m.group(0):
            pos[0] += len(m.group(0))
        length = 0.0
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            m = re.match(r"-?[0-9.eE+-]+", s[pos[0]:])
            length = float(m.group(0))
            pos[0] += len(m.group(0))
        return me, length

    root, _ = node()
    return edges, nxt[0]


def _patristic_from_edges(edges, n_nodes, n_tips):
    adj = [[] for _ in range(n_nodes)]
    for a, b, w in edges:
        w = max(0.0, w)
        adj[a].append((b, w))
        adj[b].append((a, w))
    out = np.zeros((n_tips, n_tips))
    for t in range(n_tips):
        dist = np.full(n_nodes, np.nan)
        dist[t] = 0.0
        stack = [t]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if np.isnan(dist[v]):
                    dist[v] = dist[u] + w
                    stack.append(v)
        out[t] = np.nan_to_num(dist[:n_tips], nan=0.0)
    np.fill_diagonal(out, 0.0)
    return out


def ml_patristic(records, alphabet="nt", threads=1, iqtree=None):
    """Patristic distances from an IQ-TREE maximum likelihood tree.

    `records` is [(name, aligned sequence)], already aligned and all of one
    length. Returns a matrix in the order given.
    """
    from . import config as cg
    exe = iqtree or os.environ.get("CALICITYPER_IQTREE") or cg.exe("iqtree")
    n = len(records)
    if n < 3:
        return np.zeros((n, n))
    with tempfile.TemporaryDirectory(prefix="calicityper_ml_") as td:
        fa = Path(td) / "aln.fasta"
        with open(fa, "w") as fh:
            for i, (_, seq) in enumerate(records):
                fh.write(f">t{i}\n{seq}\n")
        cmd = [str(exe), "-s", str(fa), "-m", MODEL[alphabet], "-fast", "-nt", str(max(1, threads)),
               "--prefix", str(Path(td) / "out"), "--quiet", "-redo"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        tree = Path(td) / "out.treefile"
        if p.returncode != 0 or not tree.exists():
            raise RuntimeError(f"IQ-TREE failed (exit {p.returncode}) on {n} sequences:\n"
                               f"{(p.stderr or p.stdout or '')[-1500:]}")
        newick = tree.read_text()
    edges, n_nodes = _edges_from_newick(newick, n)
    return _patristic_from_edges(edges, n_nodes, n)
