"""End-to-end tests for the calicityper command line interface.

These complement tests/test_fast_paths.py (which checks the numeric kernels)
and tests/loo_validation.py (accuracy). The point here is behaviour a user
actually meets: does the CLI run, does it emit the documented columns, and
does it degrade sensibly on input that cannot be typed?

Each test needs BLAST and MAFFT on PATH and is skipped if they are absent.
"""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

from calicityper import config as cg

COLUMNS = ["query", "genogroup", "orf1_type", "orf1_support", "orf1_runners",
           "orf2_type", "orf2_support", "orf2_runners", "note",
           "orf1_sites", "orf1_nearest", "orf1_nearest_identity", "orf1_nearest_type",
           "orf1_accepted",
           "orf2_sites", "orf2_nearest", "orf2_nearest_identity", "orf2_nearest_type",
           "orf2_accepted",
           "orf1_window_cov", "orf2_window_cov", "gii4_window_cov",
           "orf1_putative_novel", "orf1_putative_novel_detail",
           "orf2_putative_novel", "orf2_putative_novel_detail"]

needs_tools = pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("mafft")),
    reason="requires blastn and mafft on PATH",
)


def run_cli(fasta, out, *extra):
    """Invoke the CLI as a subprocess, the way a user does."""
    cmd = [sys.executable, "-m", "calicityper.cli", "type", str(fasta), "-o", str(out), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, f"calicityper failed:\n{p.stdout}\n{p.stderr}"
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_novel(fasta, out, *extra):
    cmd = [sys.executable, "-m", "calicityper.cli", "novel", str(fasta), "-o", str(out), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, f"calicityper novel failed:\n{p.stdout}\n{p.stderr}"
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def reference_genomes(n):
    return list(SeqIO.parse(cg.BLAST_FASTA, "fasta"))[:n]


def test_reference_data_is_installed():
    """The data must ship with the package, not just exist in a source checkout."""
    cg.check_data()
    assert cg.BLAST_FASTA.exists()
    assert (cg.REF_DIR["orf1"] / "GI.aln").exists()
    assert (cg.REF_DIR["orf2"] / "GI.aln").exists()


@needs_tools
def test_types_known_genomes(tmp_path):
    """Two GI reference genomes should be typed as GI with high support."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(2), q, "fasta")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2")

    assert len(rows) == 2
    assert list(rows[0]) == COLUMNS, "output columns changed; README documents these"
    for r in rows:
        assert r["genogroup"] == "GI"
        assert r["orf2_type"].startswith("GI."), r
        assert float(r["orf2_support"]) >= 70.0, r


@needs_tools
def test_non_norovirus_input_is_reported_not_crashed(tmp_path):
    """Random sequence has no BLAST hit: one row, unassigned, explanatory note."""
    q = tmp_path / "junk.fa"
    q.write_text(">junk\n" + "ACGT" * 200 + "\n")
    rows = run_cli(q, tmp_path / "out.tsv")

    assert len(rows) == 1
    assert rows[0]["genogroup"] == "unassigned"
    assert rows[0]["orf1_type"] == "" and rows[0]["orf2_type"] == ""
    assert "no BLAST hit" in rows[0]["note"]


@needs_tools
def test_short_fragment_gets_overlap_note(tmp_path):
    """A 60 nt fragment is below the 100 nt ORF-overlap minimum for typing."""
    g = reference_genomes(1)[0]
    q = tmp_path / "frag.fa"
    q.write_text(f">frag\n{str(g.seq)[1000:1060]}\n")
    rows = run_cli(q, tmp_path / "out.tsv")

    assert len(rows) == 1
    r = rows[0]
    # Either no hit at all, or a hit with too little ORF overlap to type.
    assert r["orf1_type"] == "" or r["orf2_type"] == "" or "overlap" in r["note"]


@needs_tools
def test_every_input_record_yields_exactly_one_row(tmp_path):
    """Row count and order must track the input FASTA, typable or not."""
    recs = reference_genomes(2)
    q = tmp_path / "mixed.fa"
    with open(q, "w") as fh:
        SeqIO.write(recs[0], fh, "fasta")
        fh.write(">junk\n" + "ACGT" * 150 + "\n")
        SeqIO.write(recs[1], fh, "fasta")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2")

    assert len(rows) == 3
    assert [r["query"] for r in rows] == [recs[0].id, "junk", recs[1].id]
    assert rows[1]["genogroup"] == "unassigned"


@needs_tools
def test_merged_fasta_record_is_flagged_not_fatal(tmp_path):
    """A header merged into a sequence line (an interrupted download) must not
    crash MAFFT for the whole chunk: that record is reported untyped, the rest
    are typed. Regression for KJ955493.1 in the RIVM set (job 1110626)."""
    recs = reference_genomes(2)
    s = str(recs[1].seq)
    q = tmp_path / "merged.fa"
    q.write_text(f">{recs[0].id}\n{recs[0].seq}\n"
                 f">broken\n{s[:5000]}>XX000001.1 Norovirus GII, partial cds\n{s[5000:]}\n")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2")

    assert [r["query"] for r in rows] == [recs[0].id, "broken"]
    assert rows[0]["orf1_type"] or rows[0]["orf2_type"]
    assert rows[1]["genogroup"] == "unassigned"
    assert "non-nucleotide" in rows[1]["note"]


@needs_tools
def test_results_are_deterministic_for_a_fixed_seed(tmp_path):
    """The bootstrap is seeded, so two identical runs must agree exactly."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(2), q, "fasta")
    a = run_cli(q, tmp_path / "a.tsv", "--seed", "7", "--bootstraps", "20")
    b = run_cli(q, tmp_path / "b.tsv", "--seed", "7", "--bootstraps", "20")
    assert a == b


@needs_tools
def test_runner_up_types_are_reported(tmp_path):
    """The ranked runner-up types must reach the TSV, best first."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(2), q, "fasta")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2", "--bootstraps", "30")

    for r in rows:
        for region in ("orf1", "orf2"):
            if r[f"{region}_type"] in ("", "unassigned"):
                continue
            parsed = [p.split(":") for p in r[f"{region}_runners"].split("|")]
            assert 1 <= len(parsed) <= 3, r[f"{region}_runners"]
            vals = [float(v) for _, v in parsed]
            assert vals == sorted(vals, reverse=True), "runners must be ranked"
            # the winning type is the best-supported one
            assert parsed[0][0] == r[f"{region}_type"]
            assert abs(vals[0] - float(r[f"{region}_support"])) < 0.05


@needs_tools
def test_report_writes_tsv_md_json_and_html(tmp_path):
    """--report emits all three extra formats without disturbing the TSV."""
    import json

    q = tmp_path / "q.fa"
    with open(q, "w") as fh:
        SeqIO.write(reference_genomes(2), fh, "fasta")
        fh.write(">junk\n" + "ACGT" * 150 + "\n")
    out = tmp_path / "out.tsv"
    prefix = tmp_path / "rep"
    rows = run_cli(q, out, "--threads", "2", "--bootstraps", "30", "--report", str(prefix))

    assert len(rows) == 3
    for ext in ("md", "json", "html"):
        f = prefix.with_suffix("." + ext)
        assert f.exists() and f.stat().st_size > 0, f"{ext} report missing"

    doc = json.loads(prefix.with_suffix(".json").read_text())
    assert doc["tool"] == "calicityper"
    assert doc["n_queries"] == 3
    assert len(doc["results"]) == 3
    assert doc["meta"]["bootstraps"] == 30
    # runner-ups are structured in JSON, not a packed string
    typed = [r for r in doc["results"] if r["orf2_type"] not in ("", "unassigned")]
    assert typed, "expected at least one typed query"
    assert isinstance(typed[0]["orf2_runners"], list)
    assert "support" in typed[0]["orf2_runners"][0]

    page = prefix.with_suffix(".html").read_text()
    assert "calicityper report" in page
    assert "validated for" in page, "HTML must carry the GI/GII validation caveat"


def test_report_module_works_without_bokeh(tmp_path, monkeypatch):
    """Markdown and JSON must not depend on bokeh, and HTML must degrade."""
    from calicityper import report

    monkeypatch.setattr(report, "HAVE_BOKEH", False)
    rows = [dict(query="q1", genogroup="GII", orf1_type="GII.P31", orf1_support=97.0,
                 orf1_runners="GII.P31:97.0|GII.P16:2.0", orf2_type="GII.4_Sydney",
                 orf2_support=99.0, orf2_runners="GII.4_Sydney:99.0", note="")]
    assert report.alignment_figure(["a"], ["ACGT"], "t") is None
    report.to_markdown(rows, tmp_path / "r.md")
    report.to_json(rows, tmp_path / "r.json")
    report.to_html(rows, tmp_path / "r.html", figures=[])
    assert (tmp_path / "r.md").stat().st_size > 0
    assert "GII.4_Sydney" in (tmp_path / "r.html").read_text()


def test_closest_reference_ignores_uncovered_columns():
    """A partial query is compared only where it has bases: the reference that is
    identical over the covered span wins even if it differs elsewhere."""
    from calicityper import report

    q = "----ACGTACGT----"
    refs = ["TTTTACGTACGTTTTT",      # identical over the covered 8 columns
            "ACGTACGAACGTACGT"]      # one difference inside the covered span
    k, p, n, diff = report.closest_reference(q, refs)
    assert (k, p, n, diff) == (0, 0.0, 8, 0)
    assert report.closest_reference("-" * 16, refs) is None


def test_difference_figure_draws_only_differences():
    from calicityper import report
    if not report.HAVE_BOKEH:
        return
    fig = report.difference_figure("q", "--ACGTTCGT-A", "r", "AAACGTACGTCA", "t")
    cells = [r for r in fig.renderers if "rect" in type(r.glyph).__name__.lower()][0]
    # one substitution (col 7) and one query deletion (col 11), each drawn in both rows
    assert sorted(set(cells.data_source.data["col"])) == [7, 11]


def test_complete_vp1_requires_start_stop_and_length():
    from calicityper import putative_novel as nv
    lens = [540, 560]
    orf = "ATG" + "GCT" * 549 + "TAA"                      # 550 aa, complete
    assert len(nv.complete_vp1("CC" + orf + "GG", lens)) == 550
    assert nv.complete_vp1("CC" + orf[:-3] + "GCT", lens) is None      # no stop: truncated
    assert nv.complete_vp1("ATG" + "GCT" * 100 + "TAA", lens) is None  # far too short


def test_classify_follows_per_cluster_bounds():
    from calicityper import putative_novel as nv
    t = {("orf2", "genotype", "GII.4"): 0.10, ("orf2", "genotype", "*pooled*"): 0.13,
         ("orf2", "genogroup", "GII"): 0.44, ("orf2", "genogroup", "*pooled*"): 0.44}
    assert nv.classify("orf2", {"GII.4": 0.05}, None, t)[0] == nv.KNOWN
    assert nv.classify("orf2", {"GII.4": 0.20}, {"GII": 0.30}, t)[0] == nv.NOVEL_TYPE
    assert nv.classify("orf2", {"GII.4": 0.50}, {"GII": 0.50}, t)[0] == nv.NOVEL_GROUP
    # an unlisted (small/tentative) cluster falls back to the pooled bound
    assert nv.classify("orf2", {"GII.NA1": 0.12}, None, t)[0] == nv.KNOWN


def test_putative_novel_columns_and_sapovirus_not_assessed(tmp_path):
    """Every row carries the putative novel columns; non-noroviruses are not assessed."""
    import pandas as pd
    from calicityper import putative_novel as nv
    ex = Path(__file__).resolve().parent.parent / "examples" / "example.fasta"
    out = tmp_path / "t.tsv"
    run_cli(ex, out)
    d = pd.read_csv(out, sep="\t", dtype=str).fillna("").set_index("query")
    for c in ("orf1_putative_novel", "orf2_putative_novel", "orf1_putative_novel_detail", "orf2_putative_novel_detail"):
        assert c in d.columns
    assert d.loc["MN509084.1", "orf2_putative_novel"] == nv.NA
    assert d.loc["MT451991.1", "orf2_putative_novel"] == nv.KNOWN


def test_p_distance_tolerates_ambiguous_residues():
    """An ambiguous codon translates to an ambiguity letter BLOSUM62 has no row
    for ('MTA' -> 'J'), which used to raise and end a whole run. Ambiguous
    positions are skipped rather than counted as mismatches."""
    from Bio.Seq import Seq
    from calicityper import putative_novel as pn

    assert "J" in str(Seq("ATGMTAAAA").translate())          # the real input case
    assert pn.p_distance("MJKAAA", "MLKAAA", "aa") == (0.0, 5)
    assert pn.p_distance("MXKAAA", "MLKAAA", "aa") == (0.0, 5)
    assert pn.p_distance("MLKAAA", "MLKAAA", "aa") == (0.0, 6)
    assert pn.sanitise_aa("MJKBZX*-l") == "MXKXXXXXL"


def test_sapovirus_taxonomy_csvs_match_the_published_table():
    """The CSVs must stay a faithful read of Supplemental Table 1, not a
    hand-transcription: an earlier transcription put six accessions in the wrong
    row. Regenerating must be a no-op."""
    import subprocess
    root = Path(__file__).resolve().parent.parent
    csvs = [root / "validation" / "sapovirus_vp1_taxonomy.csv",
            root / "validation" / "sapovirus_rdrp_taxonomy.csv"]
    if not (root / "validation" / "12864_2024_11048_MOESM1_ESM.xlsx").exists():
        return                                  # source spreadsheet not distributed
    before = [p.read_text() for p in csvs]
    subprocess.run([sys.executable, str(root / "scripts" / "build_sapovirus_taxonomy.py")],
                   check=True, capture_output=True)
    for p, was in zip(csvs, before):
        assert p.read_text() == was, f"{p.name} differs from the published table"


def test_sapovirus_tables_have_the_published_counts():
    """53 VP1 genotypes / 156 strains and 45 P-types / 139 strains, per the paper."""
    root = Path(__file__).resolve().parent.parent
    for name, key, labels, accs in (("sapovirus_vp1_taxonomy.csv", "genotype", 53, 156),
                                    ("sapovirus_rdrp_taxonomy.csv", "ptype", 45, 139)):
        rows = [r for r in csv.DictReader(
            l for l in open(root / "validation" / name) if not l.startswith("#"))]
        assert len(rows) == accs, (name, len(rows))
        assert len({r[key] for r in rows}) == labels, (name, len({r[key] for r in rows}))
        # the one unusable accession stays flagged rather than silently dropped
        if key == "genotype":
            bad = [r for r in rows if r["note"].startswith("UNRESOLVED")]
            assert [r["accession"] for r in bad] == ["KC309412"], bad


def test_rivm_windows_resolve_to_the_published_lengths():
    """One window per gene: region C from the 5' end of VP1, and the 3' end of
    RdRp. Lengths are genogroup-dependent."""
    from calicityper import typing as ty
    for region, group, want in (("orf2", "GII", 252), ("orf2", "GI", 264),
                                ("orf1", "GII", 172), ("orf1", "GIV", 173)):
        cols, n = ty.rivm_window_columns(region, group)
        assert n == want, (region, group, n)
        assert len(cols) == want, (region, group, len(cols))
    # region C is at the 5' end of VP1; the RdRp window at the 3' end
    assert ty.rivm_window_columns("orf2", "GII")[0].min() == 0
    orf1 = ty.rivm_window_columns("orf1", "GII")[0]
    from calicityper import align
    assert orf1.max() >= len(align.refset("orf1", "GII").seqs[0]) - 20


def test_short_input_is_typed_but_not_accepted(tmp_path):
    """150 nt is below the validated range, so the type is reported and the
    acceptance column says why it should not be relied on.

    The old behaviour rejected anything under 172 nt outright, on the grounds
    that no RIVM window could be covered. The truncation harness measured that
    assumption and it was wrong: 100-171 nt typed at 96-100% precision, so the
    sequence is informative, just below the 300-site mark where leave-one-out
    F1 reaches 97%."""
    from calicityper import config as cg
    assert cg.MIN_INPUT_NT == 100 and cg.MIN_INFORMATIVE_SITES == 300
    ref = reference_genomes(1)[0]
    q = tmp_path / "short.fa"
    SeqIO.write([ref[5200:5350]], q, "fasta")       # 150 nt, inside RdRp
    r = run_cli(q, tmp_path / "out.tsv")[0]
    assert r["genogroup"] == "GI", r
    assert r["orf1_type"] not in ("", "unassigned"), r      # the call is still shown
    assert r["orf1_accepted"].startswith("rejected"), r["orf1_accepted"]
    assert "300" in r["orf1_accepted"], r["orf1_accepted"]


def test_input_below_the_alignment_minimum_is_rejected(tmp_path):
    """Under 100 nt there are too few aligned sites for the tree, measured:
    nothing was assignable at 80 nt in the truncation harness."""
    ref = reference_genomes(1)[0]
    q = tmp_path / "tiny.fa"
    SeqIO.write([ref[5200:5280]], q, "fasta")       # 80 nt
    r = run_cli(q, tmp_path / "out.tsv")[0]
    assert r["genogroup"] == "unassigned"
    assert "rejected" in r["note"], r["note"]


def test_window_coverage_is_reported_and_does_not_withhold_the_call(tmp_path):
    """A VP1 fragment missing region C is still typed, with its coverage shown.

    This tool types on every aligned column, not on region C, so sequence
    outside the window carries the signal: in the truncation harness, fragments
    displaced off the window entirely still typed at 99-100% precision. The
    coverage belongs next to the call, not in front of it."""
    ref = reference_genomes(1)[0]
    q = tmp_path / "frag.fa"
    SeqIO.write([ref[6000:6700]], q, "fasta")       # 3' VP1, no region C
    r = run_cli(q, tmp_path / "out.tsv")[0]
    assert r["orf2_window_cov"].startswith("0/"), r["orf2_window_cov"]
    assert r["orf2_type"] not in ("", "unassigned"), r
    assert r["orf2_accepted"] == "accepted", r["orf2_accepted"]


def test_rivm_window_mode_restores_the_gate(tmp_path, monkeypatch):
    """The gate still exists for replicating RIVM's own decisions."""
    monkeypatch.setenv("CALICITYPER_RIVM_WINDOW_MODE", "1")
    ref = reference_genomes(1)[0]
    q = tmp_path / "frag.fa"
    SeqIO.write([ref[6000:6700]], q, "fasta")
    r = run_cli(q, tmp_path / "out.tsv")[0]
    assert r["orf2_type"] in ("", "unassigned"), r
    assert "RIVM window" in r["note"], r["note"]


def test_gii4_variant_window_spans_p2_to_the_vp1_terminus():
    """Region C is the shell domain and carries almost no variant signal, so the
    variant window starts at P2 (aa 279) and runs to the VP1 C-terminus."""
    from calicityper import typing as ty, config as cg
    cols = ty.gii4_variant_window_columns("GII")
    assert cg.GII4_VARIANT_WINDOW_AA == (279, None)
    assert 700 < len(cols) < 900, len(cols)
    region_c = ty.rivm_window_columns("orf2", "GII")[0]
    assert not set(cols) & set(region_c), "variant window must not overlap region C"


def test_variant_window_coverage_is_reported_with_the_variant(tmp_path):
    """A fragment covering region C but little of P2 keeps its variant, with
    the variant-window coverage reported beside it.

    At an 80% floor the harness found 86 of 103 correct variant assignments
    withheld and no error prevented, so the coverage informs the reader rather
    than deciding for them. RIVM window mode still drops the suffix."""
    ex = Path(__file__).resolve().parent.parent / "examples" / "example.fasta"
    rec = [r for r in SeqIO.parse(ex, "fasta") if r.id == "MW305639.1"][0]
    q = tmp_path / "frag.fa"
    SeqIO.write([rec[5060:6400]], q, "fasta")        # region C + part of P2
    r = run_cli(q, tmp_path / "out.tsv")[0]
    assert r["orf2_type"].startswith("GII.4"), r["orf2_type"]
    cov, width = r["gii4_window_cov"].split("/")
    assert int(cov) / int(width) < 0.8
    assert "variant" in r["note"], r["note"]


def test_setup_verifies_the_bundled_reference_set():
    """Every artefact the typing path reads is checksummed in the manifest."""
    from calicityper import setup_data, config as cg
    m = setup_data.load_manifest()
    assert len(m["sha256_16"]) >= 30, len(m["sha256_16"])
    ok, problems = setup_data.verify(cg.DATA, m, quiet=True)
    assert ok, problems[:5]


def test_setup_refuses_an_archive_that_does_not_match(tmp_path):
    """A reference set that differs from the manifest must not be installed:
    the typing windows are alignment column ranges, so a different set changes
    calls without any other symptom."""
    import shutil, tarfile
    from calicityper import setup_data, config as cg
    src = tmp_path / "src"
    shutil.copytree(cg.DATA, src)
    with open(src / "reference_sequences_VP1" / "all.orf", "a") as fh:
        fh.write("\n>tampered\nACGT\n")
    arch = tmp_path / "bad.tar.gz"
    with tarfile.open(arch, "w:gz") as tar:
        for item in sorted(src.iterdir()):
            tar.add(item, arcname=item.name)
    target = tmp_path / "out"
    rc = setup_data.install(arch, target, setup_data.load_manifest())
    assert rc == 1
    assert not target.exists(), "nothing may be installed when verification fails"


def test_setup_merges_into_a_directory_that_already_has_the_manifest(tmp_path):
    """The package data directory always holds the manifest, accession list and
    bounds table. Treating "not empty" as a reason to refuse made
    `calicityper setup` exit 1 on every fresh install, which broke CI."""
    import shutil, tarfile
    from calicityper import setup_data, config as cg

    # a target that looks like a fresh checkout: small files, no reference set
    target = tmp_path / "data"
    target.mkdir()
    for name in ("reference_manifest.json", "reference_accessions.tsv",
                 "putative_novel_bounds.tsv"):
        if (cg.DATA / name).exists():
            shutil.copy(cg.DATA / name, target / name)

    arch = tmp_path / "ref.tar.gz"
    with tarfile.open(arch, "w:gz") as tar:
        for d in ("database", "reference_sequences_RdRp", "reference_sequences_VP1"):
            if (cg.DATA / d).exists():
                tar.add(cg.DATA / d, arcname=d)

    assert setup_data.install(arch, target, setup_data.load_manifest()) == 0
    ok, problems = setup_data.verify(target, quiet=True)
    assert ok, problems[:3]
    # the small files must survive the merge
    assert (target / "reference_manifest.json").exists()
    assert (target / "reference_accessions.tsv").exists()


@needs_tools
def test_novel_command_screens_complete_sequence_only(tmp_path):
    """`calicityper novel` applies Chhabra et al. 2019 to what the paper
    defines it on, and declines anything else with the reason.

    A complete reference genome is assessable in both regions and lands inside
    its own cluster's 2x SD bound; a region C fragment has neither a complete
    VP1 reading frame nor the RdRp typing region, so it is `not_assessable`
    rather than screened against a partial sequence.
    """
    ref = reference_genomes(1)[0]
    q = tmp_path / "q.fa"
    frag = ref[5400:5700]
    frag.id, frag.description = "regionC_300nt", ""
    SeqIO.write([ref, frag], q, "fasta")
    rows = run_novel(q, tmp_path / "novel.tsv")
    by = {(r["query"], r["region"]): r for r in rows}
    whole = by[(ref.id, "orf2")]
    assert whole["flag"] == "known", whole
    assert float(whole["distance"]) < float(whole["bound"]), whole
    assert "VP1 aa" in whole["detail"], whole["detail"]
    assert "amino acid" in by[(ref.id, "orf2")]["basis"]
    assert "nucleotide" in by[(ref.id, "orf1")]["basis"]      # P-types are nt, not aa
    for region in ("orf1", "orf2"):
        assert by[("regionC_300nt", region)]["flag"] == "not_assessable"


@needs_tools
def test_novel_command_refuses_non_norovirus(tmp_path):
    """The 2x SD bounds are norovirus bounds; sapovirus has its own published
    cut-offs, so it is declined rather than measured against the wrong ones."""
    sapo = [r for r in SeqIO.parse(cg.BLAST_FASTA, "fasta") if "Sapovirus" in r.id][0]
    sapo.id, sapo.description = "sapovirus", ""
    q = tmp_path / "s.fa"
    SeqIO.write([sapo], q, "fasta")
    rows = run_novel(q, tmp_path / "novel.tsv")
    assert rows and all(r["flag"] == "not_assessable" for r in rows), rows
    assert all("sapovirus" in r["detail"] for r in rows), rows


@needs_tools
def test_novel_screen_runs_only_on_unaccepted_calls(tmp_path):
    """In `type`, a bootstrap-accepted call is a known type by definition, so
    the expensive screen is spent only where the bootstrap did not accept."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(1), q, "fasta")
    r = run_cli(q, tmp_path / "out.tsv")[0]
    assert r["orf2_accepted"] == "accepted", r
    assert r["orf2_putative_novel"] == "known", r
    assert "not screened" in r["orf2_putative_novel_detail"], r


def test_patristic_recovers_an_additive_tree_exactly():
    """The patristic transform must be a tree distance, not an approximation of
    one: on a matrix that IS additive, it returns the matrix unchanged."""
    import numpy as np
    from calicityper.patristic import patristic_matrix
    M = np.zeros((4, 4))
    M[0, 1] = M[1, 0] = 1 + 2                      # cherry (A,B)
    M[2, 3] = M[3, 2] = 3 + 4                      # cherry (C,D)
    for x, lx in ((0, 1), (1, 2)):
        for y, ly in ((2, 3), (3, 4)):
            M[x, y] = M[y, x] = lx + 5 + ly        # 5 = internal edge
    assert np.allclose(patristic_matrix(M), M)


def test_bounds_are_kept_per_distance_measure():
    """A bound is only meaningful in its own units, so the table carries both
    and reading it with the wrong measure is an error, not a silent mismatch."""
    import pytest as _pytest
    from calicityper import putative_novel as nv
    for measure in nv.MEASURES:                       # ml, nj, p
        t = nv.bounds(measure)
        assert ("orf2", "genotype", "*pooled*") in t, measure
    # the measures must actually differ, or the table is not doing its job
    assert nv.bounds("ml")[("orf1", "genogroup", "*pooled*")] != \
           nv.bounds("p")[("orf1", "genogroup", "*pooled*")]
    with _pytest.raises(RuntimeError):
        nv.bounds("nonexistent-measure")


def test_maximum_likelihood_patristic_is_the_default_measure():
    """Chhabra et al. apply the 2x SD criterion to patristic distances read off
    a maximum likelihood tree, so that is the default: IQ-TREE here, PhyML
    there. 'patristic' is accepted as a name for it."""
    from calicityper import putative_novel as nv
    assert nv.DISTANCE == "ml"
    assert nv.MEASURES == ("ml", "nj", "p")
