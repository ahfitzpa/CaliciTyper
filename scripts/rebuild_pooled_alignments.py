#!/usr/bin/env python3
"""Rebuild every pooled `all.aln` deterministically.

`mafft --auto` is not deterministic above roughly a hundred sequences when it
is given more than one thread: aligning the same `all.orf` twice at 16 threads
produces two different alignments. All four pooled reference alignments are
over that size, so their manifest checksums changed on every rebuild from
byte-identical input, which defeats the point of checksumming the set. The
per-genogroup alignments are well under the threshold and are deterministic at
any thread count; they are not touched here.

At one thread `--auto` is stable, which is what this script relies on, and it
proves it rather than asserting it: each alignment is built twice and the two
are compared before anything is written.

The join key between an alignment and its `.orf` is the bare accession, the
last `|`-separated field, so the two genera's differing id conventions --
norovirus keeps `GI|GI.1|M87661.2`, sapovirus writes `HM030923` -- both work.
That is checked too: a rebuild that fails to round-trip every accession is
discarded rather than installed.

    python scripts/rebuild_pooled_alignments.py            # report
    python scripts/rebuild_pooled_alignments.py --apply
"""
from __future__ import annotations

import argparse
import hashlib
import io
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import config as cg           # noqa: E402

# one thread, deliberately: see the module docstring
MAFFT = ["mafft", "--quiet", "--auto", "--thread", "1"]


def targets():
    for label, d in (("norovirus VP1", cg.REF_DIR["orf2"]),
                     ("norovirus RdRp", cg.REF_DIR["orf1"]),
                     ("sapovirus VP1", cg.SAPOVIRUS_REF_DIR["orf2"]),
                     ("sapovirus RdRp", cg.SAPOVIRUS_REF_DIR["orf1"])):
        d = Path(d)
        if (d / "all.orf").exists():
            yield label, d


def align(orf):
    r = subprocess.run(MAFFT + [str(orf)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"mafft failed on {orf}:\n{r.stderr[-400:]}")
    return r.stdout


def to_clustal(fasta_text, id_style):
    """Clustal text, with ids written the way this directory already writes them.

    Biopython's clustal writer truncates an identifier at 30 characters. That
    is harmless for sapovirus, whose longest is nine, but norovirus keeps the
    full `genogroup|type|accession` and four of its VP1 ids are longer than
    that. They are already truncated in the shipped file -- AY032605.1,
    AY502023.1, GU445325.2 and JN400623.1 lose their final character, stop
    matching `all.orf`, and are silently dropped from the pooled reference set
    by align.RefSet. So the full-id style gets a writer that does not
    truncate; the accession style keeps Biopython's, so that an alignment
    which has not changed stays byte-identical.
    """
    aln = AlignIO.read(io.StringIO(fasta_text), "fasta")
    if id_style != "full":
        for rec in aln:
            rec.id = rec.id.split("|")[-1]
            rec.name, rec.description = rec.id, ""
        out = io.StringIO()
        AlignIO.write(aln, out, "clustal")
        return out.getvalue()

    ids = [rec.id for rec in aln]
    width = max(36, max(len(i) for i in ids) + 6)
    seqs = [str(rec.seq) for rec in aln]
    lines = ["CLUSTAL X (1.81) multiple sequence alignment", "", ""]
    for start in range(0, len(seqs[0]), 50):
        for i, s_ in zip(ids, seqs):
            lines.append(f"{i.ljust(width)}{s_[start:start + 50]}")
        lines.append("")
        lines.append("")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--apply", action="store_true", help="write the alignments")
    a = ap.parse_args()

    changed = []
    for label, d in targets():
        orf, dest = d / "all.orf", d / "all.aln"
        n = sum(1 for _ in SeqIO.parse(orf, "fasta"))
        # how this directory already names sequences in its alignment
        first = AlignIO.read(dest, "clustal")[0].id if dest.exists() else ""
        style = "full" if "|" in first else "accession"

        one, two = align(orf), align(orf)
        if one != two:
            raise SystemExit(
                f"{label}: mafft is still not deterministic at one thread; "
                f"refusing to install an alignment that cannot be reproduced")
        text = to_clustal(one, style)

        # the alignment is joined to the .orf on the bare accession; a rebuild
        # that loses or mangles one would silently drop references
        want = {r.id.split("|")[-1] for r in SeqIO.parse(orf, "fasta")}
        got = {r.id.split("|")[-1] for r in AlignIO.read(io.StringIO(text), "clustal")}
        if want != got:
            raise SystemExit(
                f"{label}: {len(want - got)} accession(s) lost and "
                f"{len(got - want)} gained on rewrite, e.g. "
                f"{sorted(want ^ got)[:3]}; not installing")

        old = dest.read_bytes() if dest.exists() else b""
        new = text.encode()
        same = old == new
        print(f"  {label:16s} {n:4d} seqs  ids={style:10s} "
              f"{'unchanged' if same else 'CHANGES'}  "
              f"{hashlib.sha256(old).hexdigest()[:12]} -> "
              f"{hashlib.sha256(new).hexdigest()[:12]}")
        if not same:
            changed.append((label, dest, new))

    if not changed:
        print("\nevery pooled alignment is already the deterministic one")
        return 0
    if not a.apply:
        print(f"\n{len(changed)} would change; pass --apply to write")
        return 0
    for label, dest, new in changed:
        dest.write_bytes(new)
    print(f"\nwrote {len(changed)} alignment(s). Now rebuild the manifest and "
          f"the archive:\n  python scripts/build_reference_manifest.py --apply\n"
          f"  python scripts/build_reference_archive.py --out-dir zenodo --apply")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
