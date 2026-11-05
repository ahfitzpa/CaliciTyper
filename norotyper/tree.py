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


def tn93(codes, w=None):
    """Pairwise Tamura-Nei distance, pairwise deletion of gaps/ambiguities, site weights w (bootstrap)."""
    n, L = codes.shape
    w = np.ones(L) if w is None else w.astype(float)
    X = [(codes == k).astype(float) for k in range(4)]
    Xw = [x * w for x in X]
    valid, validw = X[0] + X[1] + X[2] + X[3], Xw[0] + Xw[1] + Xw[2] + Xw[3]
    V = validw @ valid.T
    match = sum(Xw[k] @ X[k].T for k in range(4))
    ag = Xw[0] @ X[2].T + Xw[2] @ X[0].T
    ct = Xw[1] @ X[3].T + Xw[3] @ X[1].T
    with np.errstate(divide="ignore", invalid="ignore"):
        P1, P2, Q = ag / V, ct / V, (V - match - ag - ct) / V
    cnt = np.array([(Xw[k] * 1.0).sum() for k in range(4)]) + 1e-9
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
    groups = {}
    for pos, n in enumerate(order[1:], start=1):
        groups.setdefault(types[n], 0)
        groups[types[n]] |= 1 << pos
    if not groups:
        return dict(type=None, support=0.0, supports={}, n_sites=codes.shape[1], n_refs=0)
    L = codes.shape[1]
    rng = np.random.default_rng(seed)
    counts = dict.fromkeys(groups, 0)
    for rep in range(nrep):
        w = rng.multinomial(L, np.full(L, 1.0 / L))
        splits, full = nj_splits(tn93(codes, w))
        for t, gm in groups.items():
            if is_clade(splits, full, gm | 1):
                counts[t] += 1
    sup = {t: 100.0 * c / nrep for t, c in counts.items()}
    best = max(sup, key=sup.get)
    ranked = sorted(sup.items(), key=lambda kv: -kv[1])[:3]
    return dict(type=best if sup[best] >= min_support else None, support=sup[best],
                best_type=best, supports=dict(ranked), n_sites=L, n_refs=len(refs))
