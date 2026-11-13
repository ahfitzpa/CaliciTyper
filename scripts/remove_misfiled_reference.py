#!/usr/bin/env python3
"""Remove a reference that is filed under the wrong genogroup.

LC342058.1 is a GI norovirus (RIVM: GI.PNA4 / GI.Pg, bootstrap 100, ORF2 GI.6)
filed into the GII RdRp reference set as GII.PNA4. The cause is one cell of the
published table: the strain is listed as `GI.PNA4|No24-Saitama-98-26` while its
assigned-P-type column reads GII.PNA4. The consequence is measurable -- the two
records labelled GII.PNA4 are 0.385 apart over the 762 nt window, while
AB074893.1, the genuine GII.PNA4, is 0.122 from the GII.P11 references -- so
the GII.PNA4 clade never forms and no query of that type can reach bootstrap
support. AB074893.1 scored 0% support against RIVM's 100%.

Removing a row from a CLUSTAL alignment does not move any column, so no other
reference's typing window changes; what does change is reference occupancy,
which is computed over the references present, so the benchmark is re-run after
this.

    python scripts/remove_misfiled_reference.py --accession LC342058.1 \
        --region orf1 --group GII [--apply]
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import config as cg          # noqa: E402


def drop_fasta(path, acc, apply):
    """Remove the record whose header holds `acc`. Returns True if it was there."""
    text = path.read_text().splitlines()
    out, hit, skipping = [], False, False
    for line in text:
        if line.startswith(">"):
            skipping = acc in line
            hit = hit or skipping
        if not skipping:
            out.append(line)
    if hit and apply:
        path.write_text("\n".join(out) + "\n")
    return hit


def drop_clustal(path, acc, apply):
    """Remove the accession's lines from a CLUSTAL alignment, keeping the rest
    byte for byte: every other sequence's columns stay exactly where they are."""
    lines = path.read_text().splitlines()
    out = [ln for ln in lines if not ln.startswith(acc) and acc not in ln.split(" ")[0]]
    hit = len(out) != len(lines)
    if hit and apply:
        path.write_text("\n".join(out) + "\n")
    return hit


def drop_csv(path, acc, apply):
    lines = path.read_text().splitlines()
    out = [ln for ln in lines if not ln.split(",")[0].startswith(acc)]
    hit = len(out) != len(lines)
    if hit and apply:
        path.write_text("\n".join(out) + "\n")
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--region", default="orf1", choices=["orf1", "orf2"])
    ap.add_argument("--group", required=True)
    ap.add_argument("--apply", action="store_true", help="write the changes (otherwise report)")
    a = ap.parse_args()

    d = cg.REF_DIR[a.region]
    jobs = [(d / f"{a.group}.orf", drop_fasta), (d / f"{a.group}.fasta", drop_fasta),
            (d / "all.orf", drop_fasta), (d / f"{a.group}.aln", drop_clustal),
            (d / "all.aln", drop_clustal), (d / f"{a.group}_annotation.csv", drop_csv)]
    # .dnd guide trees are legacy files the package never reads; left alone, as
    # scripts/fix_rdrp_references.py leaves them
    for path, fn in jobs:
        if not path.exists():
            print(f"  {path.name:24} absent")
            continue
        if a.apply:
            shutil.copy(path, path.with_suffix(path.suffix + ".bak"))
        hit = fn(path, a.accession, a.apply)
        print(f"  {path.name:24} {'removed' if hit and a.apply else 'found' if hit else 'not present'}")
    if not a.apply:
        print("\nreport only; pass --apply to write (a .bak is kept for each file touched)")
        return 0
    print(f"\n{a.accession} removed from the {a.group} {a.region} reference set.\n"
          f"Next: rebuild the manifest, then re-run the benchmark -- reference occupancy "
          f"is computed over the references present, so window columns can shift by one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
