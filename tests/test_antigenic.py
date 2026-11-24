"""GII.4 antigenic-site annotation.

The annotation is only as good as the residue-to-column mapping, so most of
these test that mapping, its refusals, coevolving pairs, loop indels, and
coordinate anchoring.
"""
import subprocess
import sys

import pytest

from calicityper import antigenic as ag

REF_LEN = 540          # canonical post-2002 full-length GII.4 VP1


def seq(length=REF_LEN, fill="A"):
    return fill * length


def test_site_residues_match_the_published_definitions():
    """Parra et al. 2026 J Virol, quoted in inputs/gii4_antigenic_sites.tsv."""
    assert ag.SITES["A"] == (294, 295, 296, 297, 298, 368, 372, 373)
    assert ag.SITES["C"] == (339, 340, 341, 375, 376, 377, 378)
    assert ag.SITES["D"] == (393, 394, 395, 396, 397, 398)
    assert ag.SITES["E"] == (407, 408, 409, 410, 411, 412, 413, 414)
    assert ag.SITES["G"] == (352, 355, 356, 357, 359, 364)
    assert ag.VARIABLE["E"] == (407, 412, 413, 414)
    assert ag.CONSERVED["I"] == (402, 403, 504, 505, 506)


def test_residues_map_to_columns_through_the_reference_not_raw_indices():
    """A gap before a residue must shift its column, or every site is wrong."""
    row = "-" * 10 + seq()
    cols = ag.residue_columns(row, [1, 294, 506])
    assert cols[1] == 10
    assert cols[294] == 10 + 293
    assert cols[506] == 10 + 505


def test_a_trimmed_alignment_is_refused():
    """The P2-only alignment case: its residue 1 is not VP1 residue 1."""
    trimmed = seq(260)
    with pytest.raises(ag.AntigenicError, match="full-length VP1 numbering"):
        ag.residue_columns(trimmed, ag.SITES["A"])


def test_site_c_alone_does_not_call_a_shift():
    """Kendra et al. 2021: site C is excluded from the site criterion."""
    ref = seq()
    q = list(ref)
    for r in ag.SITES["C"]:
        q[r - 1] = "W"
    for i in range(100, 160):
        q[i] = "W"
    ann = ag.annotate("q", "".join(q), "ref", ref)
    assert ann.aa_differences >= ag.MIN_AA_CHANGES
    assert ann.site_changes["C"]["n_changed"] == len(ag.SITES["C"])
    assert not ann.candidate_antigenic_shift


def test_shift_needs_both_enough_changes_and_two_of_ADEG():
    ref = seq()
    q = list(ref)
    for r in (ag.VARIABLE["A"][0], ag.VARIABLE["G"][0]):
        q[r - 1] = "W"
    few = ag.annotate("q", "".join(q), "ref", ref)
    assert not few.candidate_antigenic_shift, "two sites but far fewer than 18 changes"

    for i in range(100, 160):
        q[i] = "W"
    both = ag.annotate("q", "".join(q), "ref", ref)
    assert both.aa_differences >= ag.MIN_AA_CHANGES
    assert both.candidate_antigenic_shift
    assert "site C excluded" in both.reason


def test_one_site_is_not_enough_however_divergent():
    ref = seq()
    q = list(ref)
    for r in ag.SITES["A"]:
        q[r - 1] = "W"
    for i in range(100, 200):
        q[i] = "W"
    ann = ag.annotate("q", "".join(q), "ref", ref)
    assert ann.aa_differences >= ag.MIN_AA_CHANGES
    assert not ann.candidate_antigenic_shift


def test_the_caveat_travels_with_every_row():
    """The ~22% false-positive rate must not be separable from the result."""
    ref = seq()
    row = ag.to_row(ag.annotate("q", ref, "ref", ref))
    assert "22%" in row["caveat"]
    assert "not inferable from sequence" in row["caveat"]


def test_mismatched_alignment_lengths_are_refused():
    with pytest.raises(ag.AntigenicError, match="not from one alignment"):
        ag.annotate("q", seq(500), "ref", seq(540))


def test_immunotype_lookup_covers_all_28_parra_human_genotypes():
    assert len(ag.IMMUNOTYPE) == 28
    assert set(ag.IMMUNOTYPE.values()) == set("ABCDEFGHIJKL")
    assert ag.immunotype("GII.4") == ("G", "evolving")
    assert ag.immunotype("GII.2") == ("E", "static")
    assert ag.immunotype("GII.16") == ("D", "not_assessed")
    assert ag.immunotype("GI.3")[0] == ag.immunotype("GI.8")[0] == ag.immunotype("GI.9")[0] == "C"
    assert ag.immunotype("GII.17")[0] == ag.immunotype("GII.13")[0] == "J"
    assert ag.immunotype("GII.99") == (None, None)


def test_regime_is_the_tables_own_wording_not_a_binary_derivation():
    """`not_assessed` must survive: it covers 17 of the table's 31 genotypes.

    Deriving static/evolving from the one-entry EVOLVING tuple reported
    `static` for every unassessed genotype, asserting a regime Parra et al.
    never assigned.
    """
    for gt in ("GII.20", "GII.8", "GII.15", "GII.22"):
        assert ag.immunotype(gt)[1] == "not_assessed", gt
    assert ag.immunotype("GII.17")[1] == "static_emerging_candidate"
    assert set(ag.REGIME.values()) <= {
        "evolving", "static", "static_emerging_candidate", "not_assessed"}
    assert ag.immunotype("GII.4_Sydney") == ("G", "evolving")


def test_load_gii4_variant_references_orders_sydney_first(tmp_path):
    vp1_dir = tmp_path / "reference_sequences_VP1"
    vp1_dir.mkdir(parents=True)
    (vp1_dir / "GII.orf").write_text(
        f">GII|GII.4_US95_96|AJ004864.1\n{'GCT' * 539}TAA\n"
        f">GII|GII.4_Sydney|JX459908.1\n{'GCT' * 540}TAA\n"
    )
    refs = ag.load_gii4_variant_references(data_dir=tmp_path)
    assert refs[0][0] == "JX459908.1|Sydney"
    assert len(refs[0][1]) == 540
    assert refs[1][0] == "AJ004864.1|US95_96"
    assert len(refs[1][1]) == 539


def test_all_eight_coevolving_pairs_are_evaluated():
    """Every Kendra et al. 2021 Table 1 pair (including non-SITES 119-145, 309-459,
    395-504, 397-448, 465-483) must be mapped and detected."""
    ref = seq()
    q = list(ref)
    for a, b in ag.COEVOLVING:
        q[a - 1] = "Y"
        q[b - 1] = "Y"
    ann = ag.annotate("q", "".join(q), "ref", ref)
    assert ann.coevolving_pairs_changed == [f"{a}-{b}" for a, b in ag.COEVOLVING]


def test_internal_loop_insertion_caught_while_terminal_gaps_ignored():
    """GII.4 San Francisco Epitope A 294-295 insertion counts; terminal gaps do not."""
    ref_540 = "A" * 294 + "-" + "A" * (540 - 294)
    q_sf = "A" * 294 + "V" + "A" * (540 - 294 - 15) + "-" * 15
    ann = ag.annotate("sf", q_sf, "sydney", ref_540)
    assert ann.aa_differences == 1
    assert "294_295insV" in ann.site_changes["A"]["changes"]
    assert "[+V]" in ann.haplotype()


def test_540aa_anchor_prevents_us95_96_coordinate_shift(tmp_path):
    """Comparing against 539-aa US95_96 (gap at 394) using --anchor JX459908.1."""
    sydney = "A" * 540
    us95 = "A" * 393 + "-" + "A" * (540 - 394)
    query = list(us95)
    query[503] = "K"  # Site I residue 504
    fa = tmp_path / "gii4.faa"
    fa.write_text(f">JX459908.1|Sydney\n{sydney}\n>AJ004864.1|US95_96\n{us95}\n>q\n{''.join(query)}\n")
    out = tmp_path / "out.tsv"
    r = subprocess.run(
        [sys.executable, "-m", "calicityper.cli", "antigenic", str(fa),
         "--reference", "AJ004864.1|US95_96", "--anchor", "JX459908.1|Sydney",
         "--skip-input-check", "-o", str(out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "504A>K" in out.read_text()


def test_hbga_pocket_qc_alert_and_informational_sites():
    ref = seq()
    q = list(ref)
    q[374 - 1] = "G"
    q[309 - 1] = "N"
    ann = ag.annotate("ont_q", "".join(q), "ref", ref)
    assert ann.hbga_pocket_qc_alert is True
    assert ann.informational_changes["HBGA_CORE"] == "374A>G"
    assert ann.informational_changes["H"] == "309A>N"
    assert not ann.candidate_antigenic_shift


def test_cli_refuses_a_non_gii4_genotype(tmp_path):
    fa = tmp_path / "a.faa"
    fa.write_text(f">r\n{seq()}\n>q\n{seq()}\n")
    r = subprocess.run([sys.executable, "-m", "calicityper.cli", "antigenic",
                        str(fa), "--genotype", "GII.17", "-o", str(tmp_path / "o.tsv")],
                       capture_output=True, text=True)
    assert r.returncode == 2
    assert "530 and 550" in r.stderr


def test_cli_writes_a_tsv(tmp_path):
    ref = seq()
    q = list(ref)
    for r in (ag.VARIABLE["A"][0], ag.VARIABLE["D"][0]):
        q[r - 1] = "W"
    for i in range(100, 160):
        q[i] = "W"
    fa = tmp_path / "a.faa"
    fa.write_text(f">ref\n{ref}\n>query\n{''.join(q)}\n")
    out = tmp_path / "antigenic.tsv"
    r = subprocess.run([sys.executable, "-m", "calicityper.cli", "antigenic",
                        str(fa), "--genotype", "GII.4", "--skip-input-check",
                        "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = out.read_text()
    assert text.split("\n")[0].split("\t") == ag.TSV_COLUMNS
    assert "query" in text and "true" in text


def test_entropy_is_zero_for_an_invariant_column_and_two_bits_for_four():
    assert ag.shannon_entropy("AAAA") == 0.0
    assert repr(ag.shannon_entropy("AAAA")) == "0.0", "negative zero reads wrong"
    assert ag.shannon_entropy("ACDE") == 2.0


def test_entropy_ignores_gaps_rather_than_counting_them_as_a_state():
    assert ag.shannon_entropy("AAAA----") == 0.0
    assert ag.shannon_entropy("AAAAXXXX") == 0.0
    assert ag.shannon_entropy("----") is None


def test_entropy_is_a_property_of_the_alignment_not_of_a_pair():
    ref = seq()
    a = list(ref); a[293] = "W"
    b = list(ref); b[293] = "Y"
    rows = [ref, "".join(a), "".join(b)]
    cols = ag.residue_columns(ref, [294, 295])
    ent = ag.site_entropy(rows, cols)
    assert ent[294] > 0
    assert ent[295] == 0.0
    one = ag.annotate("a", rows[1], "ref", ref, entropy=ent)
    two = ag.annotate("b", rows[2], "ref", ref, entropy=ent)
    assert one.site_entropy == two.site_entropy


def test_conserved_cutoff_is_the_published_one():
    assert ag.CONSERVED_ENTROPY == 0.3


def test_cli_writes_a_provenance_sidecar(tmp_path):
    ref = seq()
    q = list(ref)
    for r in (ag.VARIABLE["A"][0], ag.VARIABLE["G"][0]):
        q[r - 1] = "W"
    fa = tmp_path / "a.faa"
    fa.write_text(f">ref\n{ref}\n>q\n{''.join(q)}\n")
    out = tmp_path / "ag.tsv"
    r = subprocess.run([sys.executable, "-m", "calicityper.cli", "antigenic",
                        str(fa), "--skip-input-check", "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    import json
    rec = json.loads((tmp_path / "ag.tsv.provenance.json").read_text())
    assert rec["subcommand"] == "antigenic"
    assert rec["shift_rule"]["min_aa_differences"] == 18
    assert rec["shift_rule"]["site_C_excluded"] is True
    assert rec["entropy"]["conserved_at_or_below"] == 0.3
    assert "not that corpus" in rec["entropy"]["measured_over"]
    assert "Parra" in rec["site_definitions"] and "Kendra" in rec["shift_rule"]["source"]


# --- input guard: norovirus GII VP1 only -----------------------------------

def _first_vp1(virus):
    """One translated VP1 per genogroup from the bundled reference sets."""
    out = {}
    for v, gg, acc, _kmers in ag.vp1_genogroup_panel():
        if v == virus:
            out.setdefault(gg, acc)
    return out


def _translated(virus, genogroup):
    from Bio import SeqIO

    from calicityper import config as cg

    ref_dir = cg.REF_DIR["orf2"] if virus == "norovirus" else cg.SAPOVIRUS_REF_DIR["orf2"]
    for rec in SeqIO.parse(str(ref_dir / "all.orf"), "fasta"):
        if rec.id.split("|")[0] == genogroup:
            return rec.id.split("|")[2], ag._to_vp1_aa(str(rec.seq))
    raise AssertionError(f"no {virus} {genogroup} reference")


def test_guard_accepts_a_norovirus_gii_vp1():
    name, aa = _translated("norovirus", "GII")
    acc, cont = ag.check_gii_input(name, aa)
    assert cont > 0.5 and acc


def test_guard_accepts_gii_that_is_not_gii4():
    # the site coordinates are resolved through a GII.4 anchor row, so any GII
    # VP1 is in scope -- only non-GII input is refused
    from Bio import SeqIO

    from calicityper import config as cg

    for rec in SeqIO.parse(str(cg.REF_DIR["orf2"] / "all.orf"), "fasta"):
        gg, ty, acc = rec.id.split("|")[:3]
        if gg == "GII" and not ty.startswith("GII.4"):
            ag.check_gii_input(acc, ag._to_vp1_aa(str(rec.seq)))
            return
    raise AssertionError("no non-GII.4 GII reference to test")


def test_guard_refuses_a_sapovirus_vp1():
    name, aa = _translated("sapovirus", "GI")
    with pytest.raises(ag.AntigenicError) as e:
        ag.check_gii_input(name, aa)
    assert "sapovirus" in str(e.value) and "not norovirus GII" in str(e.value)


def test_guard_refuses_a_norovirus_gi_vp1():
    name, aa = _translated("norovirus", "GI")
    with pytest.raises(ag.AntigenicError) as e:
        ag.check_gii_input(name, aa)
    assert "GI" in str(e.value)


def test_guard_refuses_something_that_is_not_a_vp1():
    with pytest.raises(ag.AntigenicError) as e:
        ag.check_gii_input("junk", "W" * 600)
    assert "matches no bundled VP1 reference" in str(e.value)


# ------------------------------- radicality, mutated-site count, P-dimer export
def test_substitution_score_is_symmetric_and_rejects_non_residues():
    assert ag.substitution_score("A", "A") > ag.substitution_score("A", "W")
    assert ag.substitution_score("K", "R") == ag.substitution_score("R", "K")
    # indel tokens and empty states are not scorable pairs, and must not be 0
    assert ag.substitution_score("ins", "A") is None
    assert ag.substitution_score("A", "del") is None
    assert ag.substitution_score("", "A") is None


def test_site_radicality_excludes_indels_rather_than_scoring_them_zero():
    """An indel at a site must not drag the mean toward zero."""
    only_indel = {294: ("ins", "A"), 295: ("A", "del")}
    assert ag.site_radicality(only_indel) is None
    mixed = {294: ("A", "W"), 295: ("ins", "A")}
    assert ag.site_radicality(mixed) == ag.site_radicality({294: ("A", "W")})


def test_mutated_site_count_counts_sites_not_residues():
    """Two changes inside one site are one mutated site, not two."""
    ref = seq()
    q = list(ref)
    for res in (294, 295):                       # both in Site A
        q[res - 1] = "W"
    ann = ag.annotate("q", "".join(q), "r", ref)
    assert ann.site_changes["A"]["n_changed"] == 2
    assert ann.mutated_site_count == 1
    assert ann.site_radicality["A"] is not None
    assert ann.site_radicality["C"] is None


def test_annotation_carries_the_immunotype_for_its_genotype():
    ann = ag.annotate("q", seq(), "r", seq())
    assert (ann.immunotype, ann.evolutionary_regime) == ("G", "evolving")


def test_type_immunotype_columns_are_norovirus_only():
    """Sapovirus genotype names collide with norovirus ones.

    Sapovirus GI.1 is not norovirus GI.1, and the table is norovirus, so the
    lookup is gated on the genogroup rather than on the genotype string.
    """
    from calicityper import cli

    rows = {
        0: {"query": "a", "genogroup": "GII", "orf2_type": "GII.4"},
        1: {"query": "b", "genogroup": "GII", "orf2_type": "GII.4_Sydney"},
        2: {"query": "c", "genogroup": "GI", "orf2_type": "GI.3"},
        3: {"query": "d", "genogroup": "Caliciviridae_Sapovirus_GI",
            "orf2_type": "GI.1"},
        4: {"query": "e", "genogroup": "GII", "orf2_type": "unassigned"},
        5: {"query": "f", "genogroup": "GII", "orf2_type": "untypeable"},
        6: {"query": "g", "genogroup": "unassigned", "orf2_type": ""},
    }
    cli._apply_immunotype(rows)

    assert (rows[0]["immunotype"], rows[0]["evolutionary_regime"]) == ("G", "evolving")
    # a GII.4 variant label still resolves to its genotype
    assert rows[1]["immunotype"] == "G"
    assert (rows[2]["immunotype"], rows[2]["evolutionary_regime"]) == ("C", "static")
    for i in (3, 4, 5, 6):
        assert rows[i]["immunotype"] == "", rows[i]["query"]
        assert rows[i]["evolutionary_regime"] == "", rows[i]["query"]


def test_p_domain_sequence_is_resolved_through_the_anchor_not_raw_indices():
    anchor = seq()
    row = "-" * 7 + anchor                       # same leading gap in both
    got = ag.p_domain_sequence(row, row)
    start, end = ag.P_DOMAIN
    assert got is not None and len(got) == end - start + 1


def test_p_domain_sequence_refuses_an_anchor_that_does_not_reach_the_span():
    short = seq(ag.P_DOMAIN[1] - 5)
    assert ag.p_domain_sequence(short, short) is None


def test_pdimer_export_only_flagged_deduplicated_and_in_chain_form(tmp_path):
    ref = seq()
    q = list(ref)
    for res in range(1, 40):                     # push past MIN_AA_CHANGES
        q[res - 1] = "W"
    for res in (294, 393, 407, 352):             # variable changes at A/D/E/G
        q[res - 1] = "W"
    row = "".join(q)

    flagged = ag.annotate("q1", row, "r", ref)
    same = ag.annotate("q2", row, "r", ref)      # identical P domain
    assert flagged.candidate_antigenic_shift

    quiet = ag.annotate("q3", ref, "r", ref)
    assert not quiet.candidate_antigenic_shift

    recs = ag.pdimer_records([flagged, same, quiet],
                             {"q1": row, "q2": row, "q3": ref}, ref)
    assert len(recs) == 1, "identical P domains must collapse to one record"
    assert "q1" in recs[0][0] and "q2" in recs[0][0]

    out = tmp_path / "pd.fasta"
    assert ag.write_pdimer_fasta(out, recs) == 1
    body = [l for l in out.read_text().splitlines() if not l.startswith(">")][0]
    assert body.count(":") == 1, "one colon separates the two chains"
    left, right = body.split(":")
    assert left == right, "a homodimer repeats the same chain"
