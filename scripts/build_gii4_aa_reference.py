#!/usr/bin/env python3
"""Build GII.4 translated-VP1 amino-acid reference alignments (full VP1 + aa 279 P2 window).

Follows the tool's maintenance-script convention (cf. scripts/build_putative_novel_bounds.py):
reads packaged ``calicityper/data/gii4_variant_references.tsv`` and the installed reference
set (``CALICITYPER_DATA/reference_sequences_VP1/GII.orf`` or ``GII.fasta``), orders the
canonical 540-aa ``JX459908.1|Sydney`` anchor first, and writes:
  - gii4_vp1_aa_full.fasta / .aln  (VP1 aa 1..C-term; used by `calicityper antigenic`)
  - gii4_vp1_aa_279.fasta  / .aln  (VP1 aa 279..C-term; used by F1 variant-window LOO)
  - gii4_variant_index.tsv
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import antigenic as ag  # noqa: E402
from calicityper import config as cg     # noqa: E402

P2_START_AA = 279


def _write_fasta(path: Path, recs: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for hdr, seq in recs:
            fh.write(f">{hdr}\n{seq}\n")
    print(f"wrote {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=None, help="Override CALICITYPER_DATA")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output directory [default: <CALICITYPER_DATA>/reference_sequences_VP1]")
    ap.add_argument("--p2-start-aa", type=int, default=P2_START_AA)
    args = ap.parse_args(argv)

    data_dir = args.data_dir or cg.DATA
    out_dir = args.out_dir or (Path(data_dir) / "reference_sequences_VP1")
    full_refs = ag.load_gii4_variant_references(data_dir=data_dir)
    p2_refs = [(h, s[args.p2_start_aa - 1:]) for h, s in full_refs if len(s) >= args.p2_start_aa]

    full_fa = out_dir / "gii4_vp1_aa_full.fasta"
    p2_fa = out_dir / f"gii4_vp1_aa_{args.p2_start_aa}.fasta"
    idx_tsv = out_dir / "gii4_variant_index.tsv"

    _write_fasta(full_fa, full_refs)
    _write_fasta(p2_fa, p2_refs)

    with idx_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["accession", "variant", "full_vp1_aa_length", "aa_window_length", "is_540aa_anchor"],
            delimiter="\t",
        )
        w.writeheader()
        for hdr, seq in full_refs:
            acc, var = hdr.split("|", 1)
            w.writerow({
                "accession": acc,
                "variant": var,
                "full_vp1_aa_length": len(seq),
                "aa_window_length": len(seq[args.p2_start_aa - 1:]),
                "is_540aa_anchor": str(len(seq) == 540).lower(),
            })
    print(f"wrote {idx_tsv}")

    mafft = shutil.which("mafft")
    if mafft:
        for fa in (full_fa, p2_fa):
            aln = fa.with_suffix(".aln")
            with aln.open("w", encoding="utf-8") as out:
                subprocess.run([mafft, "--quiet", "--amino", "--auto", str(fa)], stdout=out, check=True)
            print(f"wrote {aln}")
    else:
        print("MAFFT not on PATH — wrote FASTA files only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
