"""Patristic distances from a neighbour-joining tree.

Chhabra et al. 2019 apply the 2x SD criterion to *patristic* distances -- the
sum of branch lengths along the path between two tips of a phylogeny -- not to
the observed pairwise distances this package measured before. The two differ in
a way that matters for a clustering criterion: an observed distance saturates
as substitutions accumulate at the same sites, while branch lengths fitted to a
tree add up along a path, so distant clusters separate further and the spread
*within* a tight cluster is re-expressed in the same units.

One honest difference from the paper remains, and it is in the tree, not the
distance: the paper infers maximum likelihood trees with PhyML, while this
builds a neighbour-joining tree from the same distance matrix. NJ is a
least-squares fit to additivity rather than a likelihood optimum, so the branch
lengths are close but not identical, and a bound rebuilt from an ML tree would
differ slightly. NJ is used because it needs no extra program -- the typing
path already builds NJ trees -- and because it is deterministic. If PhyML or
FastTree is available, a tree from either can be substituted by replacing
`patristic_matrix`.

Patristic distance between two tips depends on which taxa are in the tree, so
the bounds and the query must be read off ONE tree. `patristic_matrix` is
therefore given the whole matrix -- references plus the query -- at once, and
the bound table is rebuilt from the same function.
"""
from __future__ import annotations

import numpy as np


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
