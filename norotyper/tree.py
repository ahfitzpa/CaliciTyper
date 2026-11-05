"""Distance, neighbour-joining and bootstrap support, vectorised with numpy.

Follows the RIVM typing tool: profile alignment -> neighbour-joining tree from HKY85-type distances ->
100 bootstrap replicates; a query is assigned to a type only if it clusters with that type's reference
sequences in at least `min_support` % of the replicates.

NOTE: RIVM uses PAUP*; here the distance is the Tamura-Nei (1993) analytic distance, which contains HKY85
as a special case. Agreement with the RIVM output is measured by scripts/validate_vs_rivm.py.
"""
import numpy as np

MAX_DIST = 5.0


def encode(seqs):
    out = np.full((len(seqs), len(seqs[0])), 4, dtype=np.int8)
    for i, s in enumerate(seqs):
        b = np.frombuffer(s.upper().replace("U", "T").encode(), dtype=np.uint8)
        for k, ch in enumerate(b"ACGT"):
            out[i, b == ch] = k
    return out


def indicators(codes):
    """One-off precomputation for tn93: base indicator matrices (float32; counts stay exact integers)."""
    X = [(codes == k).astype(np.float32) for k in range(4)]
    return X, X[0] + X[1] + X[2] + X[3]


def tn93(codes, w=None, pre=None):
    """Pairwise Tamura-Nei distance, pairwise deletion of gaps/ambiguities, site weights w (bootstrap)."""
    n, L = codes.shape
    w = np.ones(L) if w is None else w
    X, valid = pre if pre is not None else indicators(codes)
    wf = np.asarray(w, dtype=np.float32)
    Xw = [x * wf for x in X]
    V = ((valid * wf) @ valid.T).astype(float)
    match = (Xw[0] @ X[0].T + Xw[1] @ X[1].T + Xw[2] @ X[2].T + Xw[3] @ X[3].T).astype(float)
    ag = (Xw[0] @ X[2].T + Xw[2] @ X[0].T).astype(float)
    ct = (Xw[1] @ X[3].T + Xw[3] @ X[1].T).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        P1, P2, Q = ag / V, ct / V, (V - match - ag - ct) / V
    cnt = np.array([float((x @ wf).sum()) for x in X]) + 1e-9
    gA, gC, gG, gT = cnt / cnt.sum()
    gR, gY = gA + gG, gC + gT
    a1 = np.maximum(1 - gR * P1 / (2 * gA * gG) - Q / (2 * gR), 1e-12)
    a2 = np.maximum(1 - gY * P2 / (2 * gC * gT) - Q / (2 * gY), 1e-12)
    a3 = np.maximum(1 - Q / (2 * gR * gY), 1e-12)
    d = (-2 * gA * gG / gR * np.log(a1) - 2 * gC * gT / gY * np.log(a2)
         - 2 * (gR * gY - gA * gG * gY / gR - gC * gT * gR / gY) * np.log(a3))
    d = np.where(np.isfinite(d) & (V > 0), np.minimum(d, MAX_DIST), MAX_DIST)
    np.fill_diagonal(d, 0.0)
    return (d + d.T) / 2


def nj_splits(D):
    """Neighbour joining; returns the set of leaf bitmasks of every internal node (clade) of the tree."""
    D = np.array(D, dtype=float)
    m = D.shape[0]
    masks = [1 << i for i in range(m)]
    full = (1 << m) - 1
    splits = set()
    while m > 3:
        r = D.sum(1)
        Q = (m - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(np.argmin(Q)), m)
        if i > j:
            i, j = j, i
        new = 0.5 * (D[i, :] + D[j, :] - D[i, j])
        mask = masks[i] | masks[j]
        splits.add(mask)
        D[i, :] = new; D[:, i] = new; D[i, i] = 0.0
        keep = [k for k in range(m) if k != j]
        D = D[np.ix_(keep, keep)]
        masks[i] = mask
        del masks[j]
        m -= 1
    splits.update(masks)
    return splits, full


def is_clade(splits, full, mask):
    if mask == 0 or mask == full or bin(mask).count("1") == 1 or bin(full ^ mask).count("1") == 1:
        return True
    return mask in splits or (full ^ mask) in splits


try:                                   # optional accelerator: pip/conda install numba
    from numba import njit
    HAVE_NUMBA = True
except Exception:                      # pragma: no cover
    HAVE_NUMBA = False

if HAVE_NUMBA:
    @njit(cache=True)
    def _nj_clades(D, G):
        """NJ on distance matrix D (n x n). G: (g x n) 0/1 membership of candidate groups.
        Returns bool[g]: is group g (or its complement) a clade of the unrooted NJ tree."""
        n = D.shape[0]
        act = np.arange(n)
        m = n
        M = np.zeros((n, n), dtype=np.uint8)           # membership of each active node (row = node id)
        for i in range(n):
            M[i, i] = 1
        rec = np.zeros((2 * n, n), dtype=np.uint8)
        nrec = 0
        Dm = D.copy()
        r = np.zeros(n)
        while m > 3:
            for a in range(m):
                s_ = 0.0
                for b in range(m):
                    s_ += Dm[act[a], act[b]]
                r[a] = s_
            best = 1e300; bi = 0; bj = 1
            for a in range(m):
                for b in range(a + 1, m):
                    q = (m - 2) * Dm[act[a], act[b]] - r[a] - r[b]
                    if q < best:
                        best = q; bi = a; bj = b
            i = act[bi]; j = act[bj]
            dij = Dm[i, j]
            for c in range(m):
                k = act[c]
                if k != i and k != j:
                    v = 0.5 * (Dm[i, k] + Dm[j, k] - dij)
                    Dm[i, k] = v; Dm[k, i] = v
            for t in range(n):
                M[i, t] = M[i, t] | M[j, t]
            rec[nrec, :] = M[i, :]; nrec += 1
            act[bj] = act[m - 1]
            m -= 1
        for a in range(m):
            rec[nrec, :] = M[act[a], :]; nrec += 1
        out = np.zeros(G.shape[0], dtype=np.bool_)
        for g in range(G.shape[0]):
            sz = 0
            for t in range(n):
                sz += G[g, t]
            if sz <= 1 or n - sz <= 1:
                out[g] = True
                continue
            for u in range(nrec):
                same = True; comp = True
                for t in range(n):
                    if rec[u, t] != G[g, t]:
                        same = False
                    if rec[u, t] == G[g, t]:
                        comp = False
                    if not same and not comp:
                        break
                if same or comp:
                    out[g] = True
                    break
        return out


def _clades_python(D, G):
    splits, full = nj_splits(D)
    res = []
    for row in G:
        mask = sum(1 << i for i in np.flatnonzero(row))
        res.append(is_clade(splits, full, mask))
    return np.array(res)


def assign(names, types, query, codes, nrep=100, seed=1, min_support=70.0, max_refs_per_type=None):
    """names: aligned sequence names (references + query); types: {ref_name: type}; query: its name.
    Returns dict(type, support, supports, n_sites, n_refs)."""
    qi = names.index(query)
    refs = [n for n in names if n != query and n in types]
    if max_refs_per_type:
        seen, keep = {}, []
        for n in refs:
            seen[types[n]] = seen.get(types[n], 0) + 1
            if seen[types[n]] <= max_refs_per_type:
                keep.append(n)
        refs = keep
    order = [query] + refs
    idx = [names.index(n) for n in order]
    codes = codes[idx]
    n = len(order)
    tnames = sorted({types[x] for x in order[1:]})
    if not tnames:
        return dict(type=None, support=0.0, supports={}, n_sites=codes.shape[1], n_refs=0)
    G = np.zeros((len(tnames), n), dtype=np.uint8)
    for pos, nm in enumerate(order[1:], start=1):
        G[tnames.index(types[nm]), pos] = 1
    G[:, 0] = 1                                      # query + all references of the type
    L = codes.shape[1]
    # collapse identical alignment columns: the bootstrap only needs a weight per distinct column pattern
    pats, inv = np.unique(codes.T, axis=0, return_inverse=True)
    pcodes = pats.T.copy()
    inv = inv.ravel()
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(tnames))
    pre = indicators(pcodes)
    clades = _nj_clades if HAVE_NUMBA else _clades_python
    for rep in range(nrep):
        w = rng.multinomial(L, np.full(L, 1.0 / L))
        wp = np.bincount(inv, weights=w, minlength=len(pats))
        counts += clades(tn93(pcodes, wp, pre), G)
    sup = {t: 100.0 * c / nrep for t, c in zip(tnames, counts)}
    best = max(sup, key=sup.get)
    ranked = sorted(sup.items(), key=lambda kv: -kv[1])[:3]
    return dict(type=best if sup[best] >= min_support else None, support=sup[best],
                best_type=best, supports=dict(ranked), n_sites=L, n_refs=len(refs))


def warmup():
    """Compile/load the numba kernel once (in the parent process, before workers are forked)."""
    if HAVE_NUMBA:
        D = np.random.default_rng(0).random((6, 6)); D = D + D.T; np.fill_diagonal(D, 0)
        _nj_clades(D, np.ones((1, 6), dtype=np.uint8))
