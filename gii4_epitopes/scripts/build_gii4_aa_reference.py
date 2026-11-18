#!/usr/bin/env python3
"""F1: build a GII.4-only translated-VP1 amino-acid reference alignment.

Reads the verified reference set (CALICITYPER_DATA), selects the GII.4 VP1
nucleotide references listed in inputs/gii4_variant_references.tsv, extracts
their VP1 coding sequence using the bundled GII annotation, translates it,
trims to the packaged variant window (VP1 aa 279 to C-terminus), and aligns
with MAFFT (--amino --auto).

Outputs (results/, untracked):
  gii4_vp1_aa_279.fasta     unaligned translated window, per reference
  gii4_vp1_aa_279.aln       MAFFT alignment (if MAFFT available)
  gii4_variant_index.tsv    accession -> variant label + window length

Why AA-only: synonymous substitutions accumulate without changing the capsid;
a nucleotide NJ tree of GII.4 variants manufactures phantom clades (Parra
2017, GII.12 example). The GII.4-only alignment also avoids P2 insertions at
residues 393-395 being shifted by gaps from GII.3/GII.17 references.

This does NOT change how `calicityper type` calls variants.
"""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

from common import INPUTS, RESULTS, check_data_dir, read_fasta, read_tsv, translate_nt, write_fasta, write_tsv

P2_START_AA = 279   # from inputs/variant_thresholds.tsv; change there and keep in sync


def load_vp1_coordinates(annotation_path):
    coords = {}
    with open(annotation_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if "vp1_start" in r and "vp1_end" in r:
                coords[r["id"].split(".")[0]] = (int(r["vp1_start"]), int(r["vp1_end"]))
    return coords


def main():
    data = check_data_dir()
    vp1_dir = data / "reference_sequences_VP1"
    gii_fasta = vp1_dir / "GII.fasta"
    annot = vp1_dir / "GII_annotation.csv"
    for p in (gii_fasta, annot):
        if not p.exists():
            raise SystemExit(f"missing {p} in reference set {data}")

    refs = read_tsv(INPUTS / "gii4_variant_references.tsv",
                    required=["accession", "variant"])
    coords = load_vp1_coordinates(annot)
    nt = read_fasta(gii_fasta)
    aa, index = {}, []
    for r in refs:
        acc, variant = r["accession"], r["variant"]
        base = acc.split(".")[0]
        seq = nt.get(acc) or nt.get(base)
        if seq is None:
            print(f"WARN no nucleotide sequence for {acc}", file=sys.stderr)
            continue
        if base in coords:
            s, e = coords[base]
            seq = seq[s - 1:e]
        if len(seq) < 3 or (len(seq) - len(seq) % 3) < (P2_START_AA - 1) * 3:
            print(f"WARN {acc}: VP1 too short for aa {P2_START_AA}", file=sys.stderr)
            continue
        full_aa = translate_nt(seq)
        window = full_aa[P2_START_AA - 1:]
        aa[f"{acc}|{variant}"] = window
        index.append({"accession": acc, "variant": variant,
                      "aa_window_length": len(window),
                      "basis": "GII.4-only translated VP1 aa 279-Cterm"})

    write_fasta(RESULTS / "gii4_vp1_aa_279.fasta", aa)
    write_tsv(RESULTS / "gii4_variant_index.tsv", index,
              ["accession", "variant", "aa_window_length", "basis"])

    mafft = shutil.which("mafft")
    if mafft:
        aln = RESULTS / "gii4_vp1_aa_279.aln"
        with open(aln, "w", encoding="utf-8") as out:
            subprocess.run([mafft, "--quiet", "--amino", "--auto",
                            str(RESULTS / "gii4_vp1_aa_279.fasta")],
                           stdout=out, check=True)
        print(f"wrote {aln}")
    else:
        print("MAFFT not on PATH — wrote FASTA only; run `mafft --amino --auto` on the HPC.")
    print(f"Built {len(aa)} GII.4 AA references across "
          f"{len(set(r['variant'] for r in index))} variants. Review before use.")


if __name__ == "__main__":
    main()
