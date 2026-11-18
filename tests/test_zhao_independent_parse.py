"""Our parse of Supplemental Table 1 against an independently written one.

scripts/build_sapovirus_taxonomy.py reads the spreadsheet with the standard
library; scripts/third_party/convert_sapovirus_table.py reads it with openpyxl
and was written separately, to a different output shape, by another author.
Two parsers agreeing is much better evidence that the published table has been
read correctly than either one passing its own tests, because the failure mode
that matters here -- a merged label column silently attaching a strain to the
wrong type -- would not raise anything in either script on its own.

Comparing them, every difference must fall into one of three groups:

  relabelled   a deviation this project makes deliberately, recorded in
               RELABEL with its evidence
  duplicate    the supplement lists the strain twice under different types;
               the independent parse reports both joined with ';' and keeps
               neither, we take the later listing
  typo         the published accession does not resolve, or resolves to
               something that is not a calicivirus; recorded in PAPER_TYPOS

Anything outside those three is a parsing disagreement and fails.
"""
import csv
import importlib.util
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
V = ROOT / "validation"
PARSE = V / "zhao_supplement_independent_parse.csv"

needs_files = pytest.mark.skipif(
    not (PARSE.exists() and (V / "sapovirus_vp1_taxonomy.csv").exists()),
    reason="requires the taxonomy CSVs and the independent parse")


def _taxonomy():
    def read(p):
        return [r for r in csv.DictReader(l for l in open(p) if not l.startswith("#"))
                if not (r.get("note") or "").startswith(("UNRESOLVED", "RETIRED"))]
    vp1 = {r["accession"].split(".")[0]: r["genotype"]
           for r in read(V / "sapovirus_vp1_taxonomy.csv")}
    rdrp = {r["accession"].split(".")[0]: r["ptype"]
            for r in read(V / "sapovirus_rdrp_taxonomy.csv")}
    return vp1, rdrp


def _builder():
    spec = importlib.util.spec_from_file_location(
        "build_sapovirus_taxonomy", ROOT / "scripts" / "build_sapovirus_taxonomy.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def canon(label):
    """GII.NA1 and GII.NA.1 are the same name punctuated two ways."""
    return re.sub(r"\.NA\.?(\d)", r".NA.\1", (label or "").strip())


def _differences():
    vp1, rdrp = _taxonomy()
    b = _builder()
    relabelled = {(acc, lab) for lab, acc in b.RELABEL}
    out = []
    for r in csv.DictReader(open(PARSE)):
        acc = r["accession"]
        for col, mine, region in (("vp1_genotype", vp1.get(acc), "VP1"),
                                  ("p_type", rdrp.get(acc), "RdRp")):
            want = r[col].strip()
            if not want:
                continue
            if mine is None:
                out.append(dict(acc=acc, region=region, theirs=want, ours=None,
                                kind="absent"))
            elif canon(mine) != canon(want):
                kind = ("duplicate" if ";" in want
                        else "relabelled" if (acc, want) in relabelled
                        else "unexplained")
                out.append(dict(acc=acc, region=region, theirs=want, ours=mine,
                                kind=kind))
    return out, b


@needs_files
def test_no_unexplained_parsing_disagreement():
    diffs, _ = _differences()
    bad = [d for d in diffs if d["kind"] == "unexplained"]
    assert not bad, (
        "the two parsers read Supplemental Table 1 differently, and the "
        "difference is not a recorded relabel, duplicate listing or accession "
        "typo:\n" + "\n".join(
            f"  {d['acc']} {d['region']}: independent parse {d['theirs']}, ours {d['ours']}"
            for d in bad))


@needs_files
def test_strains_we_drop_are_exactly_the_published_accession_errors():
    diffs, b = _differences()
    absent = {d["acc"] for d in diffs if d["kind"] == "absent"}
    # PAPER_TYPOS replaces the accession, so the published string is gone from
    # our table; KC309412 has no resolution and is excluded from the refsets.
    expected = set(b.PAPER_TYPOS) | {"KC309412"}
    assert absent == expected, (
        f"strains in the published table but not ours: {sorted(absent)};\n"
        f"expected exactly the unresolvable published accessions: {sorted(expected)}")


@needs_files
def test_duplicate_listings_resolve_to_the_later_one():
    """The supplement lists nine strains twice, under an old broad genogroup and
    again under one it coins later in the same table. Taking the later listing
    every time is a rule, not nine separate judgements, and it matches the
    principle used elsewhere: a type the paper coins must not be left empty."""
    diffs, _ = _differences()
    dups = [d for d in diffs if d["kind"] == "duplicate"]
    assert dups, "expected the independent parse to report duplicate listings"
    wrong = [d for d in dups if canon(d["ours"]) != canon(d["theirs"].split(";")[-1])]
    assert not wrong, (
        "a duplicate listing was not resolved to the later of the two:\n" + "\n".join(
            f"  {d['acc']} {d['region']}: listed {d['theirs']}, we keep {d['ours']}"
            for d in wrong))
