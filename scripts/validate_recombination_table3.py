#!/usr/bin/env python3
"""Recombination validation against Table 3 of Zhao et al. 2024 (real strains).

The recombination path has only ever been tested on synthetic alignments with
planted breakpoints. Table 3 of the sapovirus classification lists thirteen
real strains whose VP1 and RdRp types disagree, and reports what dedicated
recombination software made of them:

  * twelve were confirmed as recombinants, "particularly between genotypes and
    specifically at the junction of the RdRp and VP1 regions";
  * the sea lion strain JN420370 was NOT: "Results from both RDP4 and Simplot
    software demonstrated no evidence of recombination in the sequence from sea
    lion", its VP1/RdRp discordance reflecting deep divergence instead.

That makes the table a validation set with a built-in negative control, and one
that needs no new sequence data. For each strain this script fetches the
genome, assembles a reference panel from the sapovirus whole genomes in the
BLAST database, runs `calicityper recomb --virus sapovirus`, and records the
status, the breakpoint interval and the region the tool names.

One strain per invocation (--index, for a SLURM array), then --summarise.

    python scripts/validate_recombination_table3.py --index 0 --out-dir results/recomb_table3
    python scripts/validate_recombination_table3.py --summarise --out-dir results/recomb_table3
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Table 3 of Zhao et al. 2024. expect_recombination follows that paper's own
# reported result, not our expectation of it.
TABLE3 = [
    # accession,   host,       genotype,  ptype,        expect_recombination
    ("JN420370",  "sea lion", "GXXX.1",  "GXIII.P1",   False),   # negative control
    ("GQ261222",  "human",    "GI.2",    "GI.P1",      True),
    ("MN794218",  "human",    "GII.4",   "GII.PNA1",   True),
    ("MK291482",  "human",    "GII.4",   "GII.P1",     True),
    ("MZ488271",  "human",    "GII.4",   "GII.P1",     True),
    ("MZ488272",  "human",    "GII.4",   "GII.P1",     True),
    ("MG012446",  "human",    "GII.4",   "GII.P1",     True),
    ("KP067444",  "human",    "GII.4",   "GII.P1",     True),
    ("MK296750",  "human",    "GII.4",   "GII.P1",     True),
    ("MH541043",  "human",    "GII.3",   "GII.P9",     True),
    ("MH541048",  "human",    "GII.3",   "GII.P9",     True),
    ("AB455793",  "human",    "GII.3",   "GII.P9",     True),
    ("LC504410",  "human",    "GII.3",   "GII.PNA2",   True),
]
MIN_PANEL = 8            # pad to at least this many candidate parents
PANEL_SIZE = 10          # reference genomes per run; OpenRDP cost grows steeply with panel size
MIN_QUERY_NT = 4000      # a genome that does not span RdRp and VP1 cannot show this junction


def fetch(acc, path):
    if path.exists() and path.stat().st_size > 1000:
        return
    url = f"{EUTILS}?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "calicityper-benchmark"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                txt = r.read().decode()
            if txt.startswith(">"):
                path.write_text(txt)
                return
        except Exception as e:                                  # NCBI drops connections
            print(f"  fetch {acc} attempt {attempt+1} failed: {e}", flush=True)
        time.sleep(5 * (attempt + 1))
    raise SystemExit(f"could not fetch {acc}")


def panel_for(query_fa, workdir, n=PANEL_SIZE):
    """The n sapovirus reference genomes most similar to the query, by BLAST."""
    from calicityper import config as cg
    out = workdir / "panel_hits.tsv"
    # -task blastn, not the megablast default: megablast is tuned for near-identical
    # sequences and finds nothing for a divergent sapovirus (JN420370 matched only
    # itself, MZ488271 nothing at all), leaving an empty reference panel.
    cmd = [cg.exe("blastn"), "-query", str(query_fa), "-subject", str(cg.BLAST_FASTA),
           "-task", "blastn", "-evalue", "1e-5",
           "-outfmt", "6 sseqid bitscore", "-max_hsps", "1"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out.write_text(p.stdout)
    seen, keep = set(), []
    for line in p.stdout.splitlines():
        sid = line.split("\t")[0]
        if "Sapovirus" in sid and sid not in seen:
            seen.add(sid)
            keep.append(sid)
    return keep[:n]


def run_one(idx, out_dir, threads):
    from Bio import SeqIO
    from calicityper import config as cg
    acc, host, genotype, ptype, expect = TABLE3[idx]
    w = Path(out_dir) / acc
    w.mkdir(parents=True, exist_ok=True)
    print(f"=== {acc} ({host}, {genotype}[{ptype}]), expect recombination: {expect}", flush=True)

    qfa = w / f"{acc}.fasta"
    fetch(acc, qfa)
    rec = next(SeqIO.parse(qfa, "fasta"))
    print(f"    query length {len(rec.seq)} nt", flush=True)
    if len(rec.seq) < MIN_QUERY_NT:
        json.dump(dict(accession=acc, skipped=f"query only {len(rec.seq)} nt"),
                  open(w / "result.json", "w"), indent=1)
        print("    SKIPPED: too short to span the RdRp/VP1 junction", flush=True)
        return 0

    wanted = panel_for(qfa, w)
    refs = {r.id: r for r in SeqIO.parse(cg.BLAST_FASTA, "fasta")}
    panel = [refs[i] for i in wanted if i in refs]
    # the query may itself be a reference genome; it must not be its own parent
    panel = [r for r in panel if acc not in r.id]
    # A recombination scan needs candidate parents. If BLAST yields too few --
    # a very divergent query -- fall back to padding with other sapovirus
    # reference genomes rather than scanning against nothing.
    if len(panel) < MIN_PANEL:
        have = {r.id for r in panel}
        for rid, rec_ in refs.items():
            if len(panel) >= MIN_PANEL:
                break
            if "Sapovirus" in rid and rid not in have and acc not in rid:
                panel.append(rec_)
        print(f"    padded panel to {len(panel)} (BLAST alone gave too few)", flush=True)
    if len(panel) < 3:
        json.dump(dict(accession=acc, skipped=f"only {len(panel)} reference genomes available"),
                  open(w / "result.json", "w"), indent=1)
        print("    SKIPPED: no usable reference panel", flush=True)
        return 0
    rec.id, rec.description = acc, ""
    combined = w / "panel.fasta"
    SeqIO.write(panel + [rec], combined, "fasta")
    print(f"    panel: {len(panel)} reference genomes", flush=True)

    o = w / "recomb"
    cmd = [sys.executable, "-m", "calicityper.cli", "recomb", "--input", str(combined),
           "--query", acc, "--virus", "sapovirus", "-o", str(o), "--threads", str(threads)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    (w / "recomb.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    jf = o / "recombination.json"
    if p.returncode != 0 or not jf.exists():
        print(f"    FAILED (exit {p.returncode}); see recomb.log", flush=True)
        json.dump(dict(accession=acc, failed=True, returncode=p.returncode),
                  open(w / "result.json", "w"), indent=1)
        return 1

    d = json.load(open(jf))
    res = d["results"] if isinstance(d, dict) and "results" in d else d
    r0 = (res if isinstance(res, list) else [res])[0]
    block = r0.get("recombination", r0)
    events = []
    for e in block.get("events", []):
        bp = (e.get("breakpoints") or [{}])[0]
        events.append(dict(event=e.get("event_id"), start=bp.get("alignment_start"),
                           end=bp.get("alignment_end"), region=e.get("biological_region"),
                           methods=e.get("supporting_methods", []),
                           plausibility=e.get("biological_plausibility")))
    out = dict(accession=acc, host=host, genotype=genotype, ptype=ptype,
               expect_recombination=expect, status=block.get("status"),
               score=block.get("score"), confidence=block.get("confidence"),
               n_events=len(events), events=events, query_nt=len(rec.seq),
               panel=len(panel))
    json.dump(out, open(w / "result.json", "w"), indent=1)
    print(f"    status={out['status']} score={out['score']} events={len(events)}", flush=True)
    for e in events[:4]:
        print(f"      {e['event']}: {e['start']}-{e['end']} | {e['region']} | "
              f"{','.join(e['methods'])}", flush=True)
    return 0


DETECTED = {"likely", "strong"}


def summarise(out_dir):
    rows = []
    for acc, host, genotype, ptype, expect in TABLE3:
        f = Path(out_dir) / acc / "result.json"
        if not f.exists():
            rows.append(dict(accession=acc, status="NOT RUN", expect=expect))
            continue
        d = json.load(open(f))
        if d.get("skipped") or d.get("failed"):
            rows.append(dict(accession=acc, status=d.get("skipped") or "FAILED", expect=expect))
            continue
        regions = {e.get("region") for e in d.get("events", []) if e.get("region")}
        rows.append(dict(accession=acc, host=host, dual=f"{genotype}[{ptype}]", expect=expect,
                         status=d.get("status"), score=d.get("score"),
                         n_events=d.get("n_events"), regions="; ".join(sorted(r for r in regions if r))))
    print(f"\n{'accession':11s} {'dual type':22s} {'expected':9s} {'status':9s} {'score':>6s} {'ev':>3s}  regions")
    agree = tested = 0
    for r in rows:
        exp = "recomb" if r.get("expect") else "NONE"
        st = str(r.get("status"))
        if st in DETECTED or st in ("none", "weak", "possible"):
            tested += 1
            called = st in DETECTED
            agree += (called == bool(r.get("expect")))
        sc = f"{r['score']:.3f}" if isinstance(r.get("score"), (int, float)) else "-"
        print(f"{r['accession']:11s} {str(r.get('dual','')):22s} {exp:9s} {st:9s} {sc:>6s} "
              f"{str(r.get('n_events','-')):>3s}  {r.get('regions','')}")
    if tested:
        print(f"\nagreement with the published result: {agree}/{tested}")
    neg = [r for r in rows if r.get("expect") is False and str(r.get("status")) in
           DETECTED | {"none", "weak", "possible"}]
    for r in neg:
        ok = str(r["status"]) not in DETECTED
        print(f"negative control {r['accession']}: status {r['status']} -- "
              f"{'as published (no recombination)' if ok else 'DISAGREES with the published result'}")
    with open(Path(out_dir) / "table3_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "host", "dual", "expect", "status",
                                           "score", "n_events", "regions"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})
    print(f"\nwrote {out_dir}/table3_summary.csv")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--index", type=int, help="row of Table 3 to run (0-12)")
    ap.add_argument("--summarise", action="store_true")
    ap.add_argument("--out-dir", default="results/recomb_table3")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    if a.summarise:
        return summarise(a.out_dir)
    if a.index is None or not 0 <= a.index < len(TABLE3):
        raise SystemExit(f"--index must be 0..{len(TABLE3)-1}")
    return run_one(a.index, a.out_dir, a.threads)


if __name__ == "__main__":
    raise SystemExit(main())
