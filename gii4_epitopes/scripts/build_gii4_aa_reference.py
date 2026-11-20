#!/usr/bin/env python3
"""F1: build GII.4-only translated-VP1 amino-acid reference alignments (full VP1 + aa 279 P2 window).

Updated from technical audit:
- Supports shipped ``reference_sequences_VP1/GII.orf`` (pipe-delimited headers) automatically
  in addition to ``GII.fasta`` + ``GII_annotation.csv`` / ``database/refseq_annotation.csv``.
- Emits BOTH ``gii4_vp1_aa_full.{fasta,aln}`` (full-length VP1 aa 1..C-term, required by
  ``calicityper antigenic``) AND ``gii4_vp1_aa_279.{fasta,aln}`` (aa 279..C-term, required
  by ``validate_variant_thresholds.py``), with canonical 540-aa ``JX459908.1|Sydney`` first.
"""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

from common import INPUTS, RESULTS, check_data_dir, read_tsv, translate_nt, write_fasta, write_tsv

P2_START_AA = 279
CANONICAL_ANCHOR = "JX459908.1"


def load_vp1_coordinates(candidates):
    coords = {}
    for annotation_path in candidates:
        if not annotation_path.exists():
            continue
        with open(annotation_path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rid = (r.get("id") or r.get("accession") or "").split(".")[0]
                if rid and r.get("vp1_start") and r.get("vp1_end"):
                    coords[rid] = (int(r["vp1_start"]), int(r["vp1_end"]))
    return coords


def read_fasta_by_accession(path):
    seqs = {}
    name, buf = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    s = "".join(buf)
                    for tok in name.split("|"):
                        seqs.setdefault(tok, s)
                        seqs.setdefault(tok.split(".")[0], s)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.upper())
        if name is not None:
            s = "".join(buf)
            for tok in name.split("|"):
                seqs.setdefault(tok, s)
                seqs.setdefault(tok.split(".")[0], s)
    return seqs


def main():
    data = check_data_dir()
    vp1_dir = data / "reference_sequences_VP1"
    nt_source = next((p for p in (vp1_dir / "GII.orf", vp1_dir / "GII.fasta") if p.exists()), None)
    if nt_source is None:
        raise SystemExit(f"missing GII.orf / GII.fasta in {vp1_dir}")

    is_spliced_orf = nt_source.suffix == ".orf"
    coords = {} if is_spliced_orf else load_vp1_coordinates([
        vp1_dir / "GII_annotation.csv",
        data / "database" / "refseq_annotation.csv",
    ])

    refs = read_tsv(INPUTS / "gii4_variant_references.tsv", required=["accession", "variant"])
    refs.sort(key=lambda r: (0 if r["accession"].strip() == CANONICAL_ANCHOR else 1))
    nt = read_fasta_by_accession(nt_source)

    aa_full, aa_279, index = {}, {}, []
    for r in refs:
        acc, variant = r["accession"].strip(), r["variant"].strip()
        base = acc.split(".")[0]
        seq = nt.get(acc) or nt.get(base)
        if seq is None:
            print(f"WARN no nucleotide sequence for {acc}", file=sys.stderr)
            continue
        if not is_spliced_orf and base in coords:
            s, e = coords[base]
            seq = seq[s - 1:e]
        usable = len(seq) - (len(seq) % 3)
        if usable < P2_START_AA * 3:
            print(f"WARN {acc}: VP1 too short for aa {P2_START_AA}", file=sys.stderr)
            continue
        full_aa = translate_nt(seq[:usable]).rstrip("*")
        window = full_aa[P2_START_AA - 1:]
        key = f"{acc}|{variant}"
        aa_full[key] = full_aa
        aa_279[key] = window
        index.append({
            "accession": acc,
            "variant": variant,
            "full_vp1_aa_length": len(full_aa),
            "aa_window_length": len(window),
            "basis": "GII.4-only translated VP1 (full + aa 279-Cterm)",
        })

    write_fasta(RESULTS / "gii4_vp1_aa_full.fasta", aa_full)
    write_fasta(RESULTS / "gii4_vp1_aa_279.fasta", aa_279)
    write_tsv(
        RESULTS / "gii4_variant_index.tsv",
        index,
        ["accession", "variant", "full_vp1_aa_length", "aa_window_length", "basis"],
    )

    mafft = shutil.which("mafft")
    if mafft:
        for stem in ("gii4_vp1_aa_full", "gii4_vp1_aa_279"):
            aln = RESULTS / f"{stem}.aln"
            with open(aln, "w", encoding="utf-8") as out:
                subprocess.run([mafft, "--quiet", "--amino", "--auto", str(RESULTS / f"{stem}.fasta")],
                               stdout=out, check=True)
            print(f"wrote {aln}")
    else:
        print("MAFFT not on PATH — wrote FASTA files only; run `mafft --amino --auto` on the HPC.")
    print(f"Built {len(aa_full)} GII.4 AA references across {len({r['variant'] for r in index})} variants.")


if __name__ == "__main__":
    main()
