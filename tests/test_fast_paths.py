"""The accelerated code paths must give the same answers as the plain ones."""
import numpy as np
import pytest
from calicityper import tree


@pytest.mark.skipif(not tree.HAVE_NUMBA, reason="numba not installed")
def test_numba_nj_matches_python():
    rng = np.random.default_rng(1)
    for _ in range(25):
        n = int(rng.integers(6, 50)); X = rng.random((n, 4))
        D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
        sp, _ = tree.nj_splits(D); sp = list(sp)
        rows = []
        for m in rng.choice(len(sp), min(5, len(sp)), replace=False):
            r = [(sp[m] >> i) & 1 for i in range(n)]; rows.append(r)          # a true clade
            r2 = list(r); r2[int(rng.integers(n))] ^= 1; rows.append(r2)      # usually not a clade
        G = np.array(rows, dtype=np.uint8)
        assert (tree._nj_clades(D, G) == tree._clades_python(D, G)).all()


def test_tn93_precomputed_equals_direct():
    rng = np.random.default_rng(2)
    codes = rng.integers(0, 5, size=(8, 300)).astype(np.int8)
    w = rng.multinomial(300, np.full(300, 1 / 300))
    assert np.allclose(tree.tn93(codes, w), tree.tn93(codes, w, tree.indicators(codes)))


def test_clades_python_matches_numba_above_63_taxa():
    """Regression: masks are built with `1 << i` over np.flatnonzero, which yields
    np.int64. `1 << np.int64(63)` wraps instead of promoting to a Python int, so
    the pure-Python fallback silently produced wrong clade masks for trees with
    64+ taxa -- and GII reference trees have 72-98. Numba is unaffected (it uses
    its own bitset), so the two paths disagreeing is the signal.
    """
    import numpy as np
    import pytest

    from calicityper import tree

    if not tree.HAVE_NUMBA:
        pytest.skip("needs numba to compare the two implementations")

    rng = np.random.default_rng(7)
    for n in (64, 72, 98):
        # two tight clusters so that real clades exist and some tests return True
        base = np.vstack([rng.normal(0, 0.01, (n // 2, 8)) + 0,
                          rng.normal(0, 0.01, (n - n // 2, 8)) + 5])
        D = np.linalg.norm(base[:, None, :] - base[None, :, :], axis=-1)
        np.fill_diagonal(D, 0.0)

        half = n // 2
        G = np.zeros((4, n), dtype=np.uint8)
        G[0, :half] = 1                 # the real first cluster: should be a clade
        G[1, half:] = 1                 # the real second cluster: should be a clade
        G[2, [0, n - 1]] = 1            # one from each cluster: should not be
        G[3, [0, 1, n - 1]] = 1         # straddles the split: should not be

        py = np.asarray(tree._clades_python(D, G)).astype(bool)
        nb = np.asarray(tree._nj_clades(D, G)).astype(bool)
        assert np.array_equal(py, nb), f"n={n}: python={py} numba={nb}"
        assert py[0] and py[1], f"n={n}: the two true clusters should be clades, got {py}"
        assert not py[2], f"n={n}: cross-cluster pair should not be a clade, got {py}"
