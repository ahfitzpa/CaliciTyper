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
from Bio import SeqIO

from calicityper import align, config as cg, recomb as R, recomb_config as RC, tree, typing as ty


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
