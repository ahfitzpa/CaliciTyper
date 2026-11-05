"""The accelerated code paths must give the same answers as the plain ones."""
import numpy as np
import pytest
from norotyper import tree


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
