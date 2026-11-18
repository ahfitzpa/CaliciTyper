"""The installed reference set must be internally consistent.

These check the data as shipped, not the code that types against it. They
exist because a reference set can be wrong in ways that every functional test
still passes: the per-genogroup files and the pooled all.orf are written in
the same run but by different code paths, and nothing previously compared
them. A strain relabelled out of its last genogroup left that genogroup's
files on disk, inside the manifest and inside the published archive, holding
an assignment the taxonomy had already retracted.
"""
import collections
import json
from pathlib import Path

import pytest
from Bio import SeqIO

from calicityper import config as cg

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "calicityper" / "data" / "reference_manifest.json"

REGIONS = {"orf2": "VP1", "orf1": "RdRp"}


def _members(d):
    """{stem: [(label, accession)]} for every per-genogroup .orf, plus all.orf."""
    per, pooled = collections.defaultdict(list), []
    for orf in sorted(Path(d).glob("*.orf")):
        for r in SeqIO.parse(orf, "fasta"):
            group, label, acc = r.id.split("|")
            (pooled if orf.stem == "all" else per[orf.stem]).append((label, acc))
    return per, pooled


def _skip_unless_installed(region):
    d = Path(cg.SAPOVIRUS_REF_DIR[region])
    if not (d / "all.orf").exists():
        pytest.skip("requires an installed sapovirus reference set")
    return d


@pytest.mark.parametrize("region", sorted(REGIONS))
def test_pooled_set_matches_the_per_genogroup_sets(region):
    """all.orf is rewritten in full each build; the per-genogroup files are not.

    Any difference means a per-genogroup file survived a build that no longer
    produces it, so the set contains two answers for the same strain.
    """
    d = _skip_unless_installed(region)
    per, pooled = _members(d)
    flat = sorted(x for members in per.values() for x in members)
    assert flat == sorted(pooled), (
        f"{REGIONS[region]}: per-genogroup .orf files and all.orf disagree.\n"
        f"  only in the per-genogroup files: {sorted(set(flat) - set(pooled))}\n"
        f"  only in all.orf: {sorted(set(pooled) - set(flat))}\n"
        "Rebuild with scripts/build_sapovirus_refsets.py, which now prunes "
        "files for groups that no longer have members.")


@pytest.mark.parametrize("region", sorted(REGIONS))
def test_no_accession_carries_two_labels(region):
    d = _skip_unless_installed(region)
    per, _ = _members(d)
    labels = collections.defaultdict(set)
    for stem, members in per.items():
        for label, acc in members:
            labels[acc].add((stem, label))
    clashes = {a: sorted(v) for a, v in labels.items() if len({l for _, l in v}) > 1}
    assert not clashes, (
        f"{REGIONS[region]}: an accession is typed two different ways by the "
        f"reference set itself: {clashes}")


@pytest.mark.parametrize("region", sorted(REGIONS))
def test_every_reference_file_is_in_the_manifest(region):
    """A file on disk but not in the manifest is unverified; the reverse is a
    broken install. `calicityper setup` checks the manifest against what it
    extracted, which cannot catch a file that was in the archive all along."""
    d = _skip_unless_installed(region)
    if not MANIFEST.exists():
        pytest.skip("no reference manifest")
    listed = set(json.loads(MANIFEST.read_text())["sha256_16"])
    for f in sorted(d.iterdir()):
        if f.suffix not in (".orf", ".fasta", ".aln"):
            continue
        rel = str(f.relative_to(ROOT)) if f.is_relative_to(ROOT) else None
        rel = rel or f"calicityper/data/{d.name}/{f.name}"
        assert rel in listed, f"{rel} is installed but not checksummed by the manifest"
