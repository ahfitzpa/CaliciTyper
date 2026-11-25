"""Distance, neighbour-joining and bootstrap support, vectorised with numpy.

Follows the RIVM typing tool: profile alignment -> neighbour-joining tree from HKY85-type distances ->
100 bootstrap replicates; a query is assigned to a type only if it clusters with that type's reference
sequences in at least `min_support` % of the replicates.

NOTE: RIVM uses PAUP*; here the distance is the Tamura-Nei (1993) analytic distance, which contains HKY85
as a special case. Agreement with the RIVM output is measured by validate_vs_rivm.py
in the benchmarking repository.
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


# Upper-triangle pair-indicator ceiling (~80 MB float32 total across 4 matrices).
# All bundled CaliciTyper reference sets are well below this (n <= 120, P <= 1500).
_MAX_BATCH_PAIR_ELEMS = 20_000_000


def pair_indicators(codes):
    """Precompute symmetric upper-triangle base-pair indicators for batched TN93.

    Returns (n, iu, ju, v_pair, m_pair, ag_pair, ct_pair, col_counts) or None if
    the alignment exceeds _MAX_BATCH_PAIR_ELEMS (caller falls back to tn93()).
    """
    n, P = codes.shape
    if n < 2:
        return None
    M = n * (n - 1) // 2
    if M * P > _MAX_BATCH_PAIR_ELEMS:
        return None
    iu, ju = np.triu_indices(n, k=1)
    ci = codes[iu]
    cj = codes[ju]
    valid_i = ci < 4
    v_pair = (valid_i & (cj < 4)).astype(np.float32)
    m_pair = ((ci == cj) & valid_i).astype(np.float32)
    ag_pair = (((ci == 0) & (cj == 2)) | ((ci == 2) & (cj == 0))).astype(np.float32)
    ct_pair = (((ci == 1) & (cj == 3)) | ((ci == 3) & (cj == 1))).astype(np.float32)
    col_counts = np.stack([(codes == k).sum(axis=0) for k in range(4)], axis=0).astype(np.float32)
    return n, iu, ju, v_pair, m_pair, ag_pair, ct_pair, col_counts


def tn93_batch(codes, W, pair_pre=None, pre=None):
    """Batched pairwise Tamura-Nei distances for a (B, P) site-weight matrix W.

    Computes all B (n x n) distance matrices in four 2D BLAS SGEMM calls over
    upper-triangle pairs, matching tn93(codes, W[b]) to floating-point equality.
    """
    Wf = np.asarray(W, dtype=np.float32)
    B = Wf.shape[0]
    if pair_pre is None:
        pair_pre = pair_indicators(codes)
    if pair_pre is None:
        pre = pre if pre is not None else indicators(codes)
        return np.stack([tn93(codes, Wf[b], pre) for b in range(B)], axis=0)
    n, iu, ju, v_pair, m_pair, ag_pair, ct_pair, col_counts = pair_pre
    WT = np.ascontiguousarray(Wf.T, dtype=np.float32)
    V_u = (v_pair @ WT).T.astype(np.float64)
    M_u = (m_pair @ WT).T.astype(np.float64)
    AG_u = (ag_pair @ WT).T.astype(np.float64)
    CT_u = (ct_pair @ WT).T.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        P1 = AG_u / V_u
        P2 = CT_u / V_u
        Q = (V_u - M_u - AG_u - CT_u) / V_u
    cnt = (Wf @ col_counts.T).astype(np.float64) + 1e-9
    g = cnt / cnt.sum(axis=1, keepdims=True)
    gA, gC, gG, gT = g[:, 0:1], g[:, 1:2], g[:, 2:3], g[:, 3:4]
    gR, gY = gA + gG, gC + gT
    a1 = np.maximum(1.0 - gR * P1 / (2.0 * gA * gG) - Q / (2.0 * gR), 1e-12)
    a2 = np.maximum(1.0 - gY * P2 / (2.0 * gC * gT) - Q / (2.0 * gY), 1e-12)
    a3 = np.maximum(1.0 - Q / (2.0 * gR * gY), 1e-12)
    du = (-2.0 * gA * gG / gR * np.log(a1)
          - 2.0 * gC * gT / gY * np.log(a2)
          - 2.0 * (gR * gY - gA * gG * gY / gR - gC * gT * gR / gY) * np.log(a3))
    du = np.where(np.isfinite(du) & (V_u > 0), np.minimum(du, MAX_DIST), MAX_DIST)
    out = np.zeros((B, n, n), dtype=np.float64)
    out[:, iu, ju] = du
    out[:, ju, iu] = du
    return out


def nj_splits(D):
    """Neighbour joining; returns the set of leaf bitmasks of every internal node (clade) of the tree."""
    D = np.array(D, dtype=float, copy=True)
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
        masks[i] = mask
        last = m - 1
        if j != last:
            D[j, :] = D[last, :]
            D[:, j] = D[:, last]
            masks[j] = masks[last]
        D = D[:last, :last]
        masks.pop()
        m = last
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
        # int(i): np.flatnonzero yields np.int64, and `1 << np.int64(63)` wraps
        # instead of promoting, so masks silently went wrong above 63 taxa.
        mask = sum(1 << int(i) for i in np.flatnonzero(row))
        res.append(is_clade(splits, full, mask))
    return np.array(res)


def nearest_reference(names, query, codes):
    """(name, identity) of the reference most identical to the query.

    Identity is counted over the columns where BOTH have a base, so a short
    query is not penalised for the reference's extra length -- it answers
    "how close is the query to this reference where they overlap", which is
    what a reader wants beside a genotype call.
    """
    qi = names.index(query)
    q = codes[qi]
    best, best_id, best_n = "", 0.0, 0
    for i, name in enumerate(names):
        if i == qi:
            continue
        r = codes[i]
        both = (q < 4) & (r < 4)
        n = int(both.sum())
        if n < MIN_IDENTITY_SITES:
            continue
        ident = float((q[both] == r[both]).sum()) / n
        if ident > best_id:
            best, best_id, best_n = name, ident, n
    return best, round(100 * best_id, 2), best_n


def second_nearest_type(names, types, query, codes, exclude_type):
    """(type, p_distance, n_sites) of the closest reference from any OTHER type.

    Quantifies assignment margin for distance-based classifiers (sapovirus):
    a small margin or a second type also below the genotype cut-off indicates
    ambiguity between closely related types (common in sapovirus RdRp).
    """
    qi = names.index(query)
    q = codes[qi]
    best_type, best_d, best_n = "", None, 0
    for i, name in enumerate(names):
        if i == qi:
            continue
        t = types.get(name, "")
        if not t or t == exclude_type:
            continue
        r = codes[i]
        both = (q < 4) & (r < 4)
        n = int(both.sum())
        if n < MIN_IDENTITY_SITES:
            continue
        d = float((q[both] != r[both]).sum()) / n
        if best_d is None or d < best_d:
            best_type, best_d, best_n = t, d, n
    if best_d is None:
        return "", None, 0
    return best_type, round(best_d, 4), best_n


MIN_IDENTITY_SITES = 50


# The acceptance threshold stays at 70%. Calibration across three datasets put
# every error at support <= 71, so 75 would be the lowest threshold with no
# errors anywhere -- but the whole margin is a single observation at 71%, and
# raising it would withhold correct calls to catch that one. Decision taken
# 2026-11-17: keep 70 and report support beside every call so a borderline
# result is visible. Revisit only with more than one observation in the gap.
FAST_BOOT_CHECKPOINT = 25
FAST_BOOT_LOW_CUTOFF = 8


def _sample_pattern_weights(rng, L, probs, inv, n_pats, count):
    """Draw `count` multinomial bootstrap weight vectors collapsed to column patterns."""
    W_raw = rng.multinomial(L, probs, size=count)
    return np.asarray(
        [np.bincount(inv, weights=w, minlength=n_pats) for w in W_raw],
        dtype=np.float32,
    )


def assign(names, types, query, codes, nrep=100, seed=1, min_support=70.0,
           max_refs_per_type=None, fast_boot=False):
    """names: aligned sequence names (references + query); types: {ref_name: type}; query: its name.
    Returns dict(type, support, supports, n_sites, n_refs)."""
    if query not in names:
        raise ValueError(f"{query!r} is not in names")
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
    n_pats = len(pats)
    rng = np.random.default_rng(seed)
    probs = np.full(L, 1.0 / L)
    counts = np.zeros(len(tnames))
    pair_pre = pair_indicators(pcodes)
    pre = None if pair_pre is not None else indicators(pcodes)
    clades = _nj_clades if HAVE_NUMBA else _clades_python

    def _run_batch(batch_n):
        Wp = _sample_pattern_weights(rng, L, probs, inv, n_pats, batch_n)
        Ds = tn93_batch(pcodes, Wp, pair_pre=pair_pre, pre=pre)
        for rep in range(batch_n):
            counts[:] += clades(Ds[rep], G)

    if fast_boot and nrep > FAST_BOOT_CHECKPOINT:
        b0 = FAST_BOOT_CHECKPOINT
        _run_batch(b0)
        best_c = int(counts.max()) if len(counts) else 0
        stop_high = (best_c == b0 and min_support <= 95.0)
        stop_low = (best_c <= FAST_BOOT_LOW_CUTOFF and min_support >= 70.0)
        if stop_high or stop_low:
            used_rep = b0
        else:
            _run_batch(nrep - b0)
            used_rep = nrep
    else:
        if nrep > 0:
            _run_batch(nrep)
        used_rep = max(1, nrep)

    sup = {t: 100.0 * c / used_rep for t, c in zip(tnames, counts)}
    best = max(sup, key=sup.get)
    ranked = sorted(sup.items(), key=lambda kv: -kv[1])[:3]
    return dict(type=best if sup[best] >= min_support else None, support=sup[best],
                best_type=best, supports=dict(ranked), all_supports=sup,
                n_sites=L, n_refs=len(refs), n_replicates=used_rep)


def warmup():
    """Compile/load the numba kernel once (in the parent process, before workers are forked)."""
    if HAVE_NUMBA:
        D = np.random.default_rng(0).random((6, 6)); D = D + D.T; np.fill_diagonal(D, 0)
        _nj_clades(D, np.ones((1, 6), dtype=np.uint8))
