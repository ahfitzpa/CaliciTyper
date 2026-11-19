"""GII.4 antigenic-site annotation.

The annotation is only as good as the residue-to-column mapping, so most of
these test that mapping and its refusals rather than the arithmetic.
"""
import subprocess
import sys

import pytest

from calicityper import antigenic as ag

REF_LEN = 540          # a plausible full-length GII.4 VP1


def seq(length=REF_LEN, fill="A"):
    return fill * length


def test_site_residues_match_the_published_definitions():
    """Parra et al. 2026 J Virol, quoted in inputs/gii4_antigenic_sites.tsv."""
    assert ag.SITES["A"] == (294, 295, 296, 297, 298, 368, 372, 373)
    assert ag.SITES["C"] == (339, 340, 341, 375, 376, 377, 378)
    assert ag.SITES["D"] == (393, 394, 395, 396, 397, 398)
    assert ag.SITES["E"] == (407, 408, 409, 410, 411, 412, 413, 414)
    assert ag.SITES["G"] == (352, 355, 356, 357, 359, 364)
    # only half of site E is variable
    assert ag.VARIABLE["E"] == (407, 412, 413, 414)
    # the mAb A1431 footprint
    assert ag.CONSERVED["I"] == (402, 403, 504, 505, 506)


def test_residues_map_to_columns_through_the_reference_not_raw_indices():
    """A gap before a residue must shift its column, or every site is wrong."""
    row = "-" * 10 + seq()
    cols = ag.residue_columns(row, [1, 294, 506])
    assert cols[1] == 10          # first residue sits after the gaps
    assert cols[294] == 10 + 293
    assert cols[506] == 10 + 505


def test_a_trimmed_alignment_is_refused():
    """The P2-only alignment case: its residue 1 is not VP1 residue 1.

    gii4_vp1_aa_279.aln starts at residue 279, so reading its first column as
    residue 1 would point every site at the wrong amino acid. Refusing is the
    only safe behaviour, because the file itself carries no offset.
    """
    trimmed = seq(260)            # ~279..540, far short of residue 506
    with pytest.raises(ag.AntigenicError, match="full-length VP1 numbering"):
        ag.residue_columns(trimmed, ag.SITES["A"])


def test_site_c_alone_does_not_call_a_shift():
    """Kendra et al. 2021: site C was enriched among pairs that were NOT
    antigenically different, so it is excluded from the site criterion."""
    ref = seq()
    q = list(ref)
    for r in ag.SITES["C"]:
        q[r - 1] = "W"
    for i in range(100, 160):     # plenty of changes elsewhere, >= 18
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


def test_immunotype_lookup_is_annotation_with_one_evolving_genotype():
    assert ag.immunotype("GII.4") == ("G", "evolving")
    assert ag.immunotype("GII.20") == ("G", "static")
    assert ag.immunotype("GI.3")[0] == ag.immunotype("GI.9")[0] == "C"
    assert ag.immunotype("GII.17")[0] == ag.immunotype("GII.13")[0] == "J"
    assert ag.immunotype("GII.99") == (None, None)
    # a variant suffix must not defeat the lookup
    assert ag.immunotype("GII.4_Sydney") == ("G", "evolving")


def test_cli_refuses_a_non_gii4_genotype(tmp_path):
    """Residue numbers must not be carried to another genotype: VP1 runs
    530-550 residues across human genotypes."""
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
                        str(fa), "--genotype", "GII.4", "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = out.read_text()
    assert text.split("\n")[0].split("\t") == ag.TSV_COLUMNS
    assert "query" in text and "true" in text


# ------------------------------------------------------------------ entropy

def test_entropy_is_zero_for_an_invariant_column_and_two_bits_for_four():
    assert ag.shannon_entropy("AAAA") == 0.0
    assert repr(ag.shannon_entropy("AAAA")) == "0.0", "negative zero reads wrong"
    assert ag.shannon_entropy("ACDE") == 2.0


def test_entropy_ignores_gaps_rather_than_counting_them_as_a_state():
    """A gap is missing data. Counting it as a 21st residue makes a ragged
    column look variable when it is only incomplete."""
    assert ag.shannon_entropy("AAAA----") == 0.0
    assert ag.shannon_entropy("AAAAXXXX") == 0.0
    assert ag.shannon_entropy("----") is None, "nothing to score is not 'conserved'"


def test_entropy_is_a_property_of_the_alignment_not_of_a_pair():
    ref = seq()
    a = list(ref); a[293] = "W"
    b = list(ref); b[293] = "Y"
    rows = [ref, "".join(a), "".join(b)]
    cols = ag.residue_columns(ref, [294, 295])
    ent = ag.site_entropy(rows, cols)
    assert ent[294] > 0, "three different residues at 294"
    assert ent[295] == 0.0, "invariant"
    # the same entropy dict annotates every query, so the column cannot differ
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
                        str(fa), "-o", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    import json
    rec = json.loads((tmp_path / "ag.tsv.provenance.json").read_text())
    assert rec["subcommand"] == "antigenic"
    assert rec["shift_rule"]["min_aa_differences"] == 18
    assert rec["shift_rule"]["site_C_excluded"] is True
    assert rec["entropy"]["conserved_at_or_below"] == 0.3
    # the entropy cut-off came from a 1,572-sequence corpus, not from this run
    assert "not that corpus" in rec["entropy"]["measured_over"]
    assert "Parra" in rec["site_definitions"] and "Kendra" in rec["shift_rule"]["source"]
