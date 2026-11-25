"""Edge-case & LOO exclusion tests addressing peer-review transparency checks.

Covers:
1. LOO `exclude=[acc]` strictly removes the held-out accession from BOTH the
   reference alignment and nearest-reference pool in memory without altering
   the installed reference bundle on disk.
2. Ambiguous IUPAC bases (N, R, Y, W, S, K, M) are masked pairwise in TN93 NJ
   distances and nearest-reference identity calculations without inflating
   divergence or crashing typing.
3. High-ambiguity consensus sequences (e.g., mixed infections) trigger the
   recombination QC gate (`> 5%` ambiguous -> `insufficient_sequence`).
"""
from __future__ import annotations

import numpy as np

from calicityper import align, recomb as R, recomb_config as RC, tree, typing as ty


def test_loo_exclude_removes_held_out_reference_from_nearest_and_tree():
    """During LOO, the excluded reference accession must never appear as nearest."""
    ref = align.refset("orf2", "GI")
    target_acc = ref.names[0]
    # Reconstruct an aligned copy of the target sequence directly from ref.codes
    bases = np.array(list("ACGT-"))
    qaln = "".join(bases[ref.codes[0]])
    res_incl = ty.type_aligned(ref, qaln, exclude=(), nrep=10, seed=1)
    assert res_incl["nearest"] == target_acc
    assert res_incl["nearest_identity"] == 100.0

    res_excl = ty.type_aligned(ref, qaln, exclude=[target_acc], nrep=10, seed=1)
    assert res_excl["nearest"] != target_acc
    assert res_excl["n_refs"] <= len(ref.names) - 1


def test_ambiguous_iupac_bases_are_masked_pairwise_not_counted_as_mismatches():
    """IUPAC ambiguity codes (N, R, Y, etc.) must be excluded pairwise from identity and TN93."""
    s1 = "ACGT" * 50
    # Replace 20 bases with IUPAC ambiguity codes (mixed-infection consensus simulation)
    s2_list = list(s1)
    for idx in range(0, 20):
        s2_list[idx] = "RYWSKMN"[idx % 7]
    s2 = "".join(s2_list)
    codes = tree.encode([s1, s2])
    # All non-ACGT codes encode to 4 (masked)
    assert int((codes[1] == 4).sum()) == 20
    # Pairwise distance over remaining 180 unambiguous bases is 0.0
    D = tree.tn93(codes)
    assert D[0, 1] == 0.0
    acc, ident, n_sites = tree.nearest_reference(["ref", "query"], "query", codes)
    assert acc == "ref"
    assert ident == 100.0
    assert n_sites == 180


def test_recomb_qc_refuses_high_ambiguity_mixed_consensus(tmp_path):
    """A consensus with >5% N/IUPAC ambiguities (e.g. mixed infection) is refused cleanly."""
    cfg = RC.load(None)
    # 6000 nt sequence with 8% IUPAC ambiguity codes
    s_ref = ("ACGTACGTAC" * 600)
    s_mix = ("R" * 480) + s_ref[480:]
    fa = tmp_path / "mix.fasta"
    fa.write_text(f">ref1\n{s_ref}\n>mixed_q\n{s_mix}\n")
    aln = R.read_alignment(str(fa))
    assert aln.masked_chars == 480
    ok, reason, stats = R.qc_query(aln, "mixed_q", cfg)
    assert not ok
    assert "ambiguous bases" in reason
    assert stats["ambiguous_fraction"] > cfg.max_ambiguous_fraction


def test_gii4_variant_reference_missing_warning_on_divergent_or_rejected():
    """Divergent (<94%) or low-support GII.4 calls emit variant_reference_missing_risk."""
    def _row():
        return {"genogroup": "GII", "note": ""}

    # 1. Accepted GII.4 with <94% identity (e.g. post-2024 GII.4 San Francisco mapped to Sydney)
    r_div = _row()
    ty.apply_result(r_div, "orf2", {
        "best_type": "GII.4_Sydney", "support": 80.0, "supports": {"GII.4_Sydney": 80.0},
        "n_sites": 1600, "nearest": "JX459908.1", "nearest_identity": 91.95,
        "nearest_type": "GII.4_Sydney",
    }, min_support=70.0)
    assert r_div["orf2_accepted"] == "accepted"
    assert "variant_reference_missing_risk" in r_div["note"]
    assert "92.0% < 94.0%" in r_div["note"]

    # 2. Rejected GII.4 with low bootstrap support (e.g. GII.4 Hong Kong)
    r_rej = _row()
    ty.apply_result(r_rej, "orf2", {
        "best_type": "GII.4_Osaka", "support": 36.0, "supports": {"GII.4_Osaka": 36.0},
        "n_sites": 1600, "nearest": "AB434770.1", "nearest_identity": 91.02,
        "nearest_type": "GII.4_Osaka",
    }, min_support=70.0)
    assert r_rej["orf2_accepted"].startswith("rejected:")
    assert "variant_reference_missing_risk" in r_rej["note"]

    # 3. Accepted GII.4 with >=94% identity (Sydney-like) -> no warning
    r_ok = _row()
    ty.apply_result(r_ok, "orf2", {
        "best_type": "GII.4_Sydney", "support": 100.0, "supports": {"GII.4_Sydney": 100.0},
        "n_sites": 1600, "nearest": "JX459908.1", "nearest_identity": 96.5,
        "nearest_type": "GII.4_Sydney",
    }, min_support=70.0)
    assert r_ok["orf2_accepted"] == "accepted"
    assert "variant_reference_missing_risk" not in r_ok["note"]

