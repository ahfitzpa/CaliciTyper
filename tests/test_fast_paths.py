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


def test_tn93_batch_matches_single_bit_for_bit():
    rng = np.random.default_rng(11)
    codes = rng.integers(0, 5, size=(16, 250)).astype(np.int8)
    W = rng.multinomial(250, np.full(250, 1.0 / 250), size=20).astype(np.float32)
    pre = tree.indicators(codes)
    expected = np.stack([tree.tn93(codes, W[b], pre) for b in range(20)], axis=0)
    got = tree.tn93_batch(codes, W)
    assert np.allclose(expected, got, atol=1e-12, rtol=1e-12)


def test_assign_fast_boot_stops_early_and_preserves_clear_calls():
    # Two well-separated clusters: query matches cluster A cleanly -> 25/25 early stop.
    rng = np.random.default_rng(3)
    L = 300
    base_a = rng.integers(0, 4, size=L, dtype=np.int8)
    base_b = (base_a + 1) % 4
    base_c = (base_a + 2) % 4
    seqs = [base_a.copy(), base_a.copy(), base_a.copy(),
            base_b.copy(), base_b.copy(), base_c.copy(), base_c.copy()]
    # small divergence within clusters
    seqs[1][:5] = (seqs[1][:5] + 1) % 4
    seqs[2][5:10] = (seqs[2][5:10] + 1) % 4
    seqs[4][:5] = (seqs[4][:5] + 1) % 4
    seqs[6][:5] = (seqs[6][:5] + 1) % 4
    codes = np.vstack(seqs)
    names = ["QUERY", "a1", "a2", "b1", "b2", "c1", "c2"]
    types = {"a1": "GII.4", "a2": "GII.4", "b1": "GII.2", "b2": "GII.2", "c1": "GII.3", "c2": "GII.3"}

    full = tree.assign(names, types, "QUERY", codes, nrep=100, seed=1, fast_boot=False)
    fast = tree.assign(names, types, "QUERY", codes, nrep=100, seed=1, fast_boot=True)
    assert full["type"] == "GII.4" and full["n_replicates"] == 100
    assert fast["type"] == "GII.4" and fast["n_replicates"] == tree.FAST_BOOT_CHECKPOINT


def test_cli_and_typing_do_not_eagerly_import_pandas():
    import subprocess, sys
    code = (
        "import sys; import calicityper.cli, calicityper.typing; "
        "assert 'pandas' not in sys.modules, 'pandas eagerly imported'"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


def test_polcap_preset_and_region_specific_min_sites():
    import argparse
    from calicityper import cli, config as cg, typing as ty

    orig = (cg.MIN_INFORMATIVE_SITES, cg.MIN_INFORMATIVE_SITES_ORF1, cg.MIN_INFORMATIVE_SITES_ORF2)
    try:
        ns = argparse.Namespace(preset="polcap", min_sites=None, min_sites_orf1=None, min_sites_orf2=None)
        cli._apply_runtime_overrides(ns)
        assert cg.min_sites_for("orf1") == 170
        assert cg.min_sites_for("orf2") == 250

        row = dict(query="q1", genogroup="GII", note="")
        ty.apply_result(row, "orf1", dict(type="GII.P16", best_type="GII.P16", support=99.0, n_sites=210), 70.0)
        ty.apply_result(row, "orf2", dict(type="GII.4", best_type="GII.4", support=99.0, n_sites=210), 70.0)
        assert row["orf1_accepted"] == "accepted"
        assert row["orf2_accepted"].startswith("rejected:") and "< 250" in row["orf2_accepted"]

        # Explicit --min-sites-orf1 overrides preset
        ns2 = argparse.Namespace(preset="polcap", min_sites=None, min_sites_orf1=190, min_sites_orf2=None)
        cli._apply_runtime_overrides(ns2)
        assert cg.min_sites_for("orf1") == 190
        assert cg.min_sites_for("orf2") == 250
    finally:
        cg.MIN_INFORMATIVE_SITES, cg.MIN_INFORMATIVE_SITES_ORF1, cg.MIN_INFORMATIVE_SITES_ORF2 = orig


def test_write_tsv_and_read_fasta_roundtrip(tmp_path):
    from calicityper import cli, typing as ty

    fa = tmp_path / "in.fa"
    fa.write_text(">seq1 extra desc\nACGT\nTGCA\n>seq2\nGGGG\n")
    assert ty.read_fasta(fa) == [("seq1", "ACGTTGCA"), ("seq2", "GGGG")]

    out = tmp_path / "out.tsv"
    cli._write_tsv([{"a": 1, "b": None}, {"a": 2, "c": "x"}], str(out))
    lines = out.read_text().splitlines()
    assert lines == ["a\tb\tc", "1\t\t", "2\t\tx"]


def test_ml_patristic_passes_constraint_to_iqtree(tmp_path, monkeypatch):
    from calicityper import patristic
    seen_cmd = {}

    def fake_run(cmd, capture_output, text):
        seen_cmd["cmd"] = list(cmd)
        prefix = cmd[cmd.index("--prefix") + 1]
        (tmp_path / "dummy").write_text("")
        from pathlib import Path
        Path(prefix + ".treefile").write_text("(t0:0.1,t1:0.2,(t2:0.1,t3:0.1):0.05);")
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        return R()

    monkeypatch.setattr(patristic.subprocess, "run", fake_run)
    recs = [("r0", "ACGT"), ("r1", "ACGG"), ("r2", "TCGT"), ("q", "TCGG")]
    P = patristic.ml_patristic(recs, alphabet="nt", iqtree="fake_iqtree",
                               constraint="(t0:0.1,t1:0.2,t2:0.1);")
    assert "-g" in seen_cmd["cmd"]
    assert P.shape == (4, 4)
    assert np.isclose(P[0, 1], 0.3)
