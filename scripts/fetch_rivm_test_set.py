#!/usr/bin/env python3
"""Fetch the GenBank sequences that the RIVM Typing Tool called, for comparison.

Reads the versioned accessions from validation/rivm_accessions.tsv (shipped in
the repository) and downloads exactly those versions from NCBI into one FASTA.
Versions matter: GenBank records are revised, and a bare accession returns the
newest version, which is not necessarily the sequence RIVM typed. RIVM's own
calls (validation/results_rivm.csv, not redistributed) are only needed to
score the run afterwards, not to fetch it.

Why not entrez-direct: `efetch` POSTs batches and NCBI drops those connections
under throttling ("SSL_read: unexpected eof"), which fails silently per batch.
Plain GET requests of 50 accessions with a polite rate are reliable, and this
script is resumable and verifies what it got.

Usage:
    python scripts/fetch_rivm_test_set.py [out.fasta] [--batch 50] [--api-key KEY]

An NCBI API key (free, from https://www.ncbi.nlm.nih.gov/account/) raises the
rate limit from 3 to 10 requests/second. Set NCBI_API_KEY or pass --api-key.
"""
import argparse
import csv
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ACCESSIONS = ROOT / "validation" / "rivm_accessions.tsv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "calicityper"


def wanted_accessions(every=False):
    """Accessions to fetch.

    By default only the records RIVM assigned a capsid or polymerase type
    (4,392), which is what the genotype comparison needs. With `every`, all
    5,986 records, including the non-norovirus caliciviruses -- needed to score
    genogroup DETECTION, since that must also reject sapoviruses, lagoviruses
    and vesiviruses correctly.
    """
    with open(ACCESSIONS) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not every:
        rows = [r for r in rows if r["genotype_set"] == "yes"]
    return sorted({r["accession"] for r in rows})


def drop_unterminated_tail(path):
    """If an interrupted run left the last record without a final newline, cut it off.

    That record may be truncated, and already_have() would count it as fetched;
    the next appended header would also land on the same line, merging two
    records (KJ955493.1 / KM025143.1, job 1097930, 2026-11-09).
    """
    if not path.exists() or path.stat().st_size == 0:
        return None
    data = path.read_bytes()
    if data.endswith(b"\n"):
        return None
    start = data.rfind(b"\n>") + 1          # 0 if the only record is the broken one
    acc = data[start + 1:].split(None, 1)[0].decode(errors="replace") if data[start:start + 1] == b">" else "?"
    with open(path, "r+b") as fh:
        fh.truncate(start)
    return acc


def already_have(path):
    """Versioned accessions already present in a partially written FASTA, for resuming."""
    if not path.exists():
        return set()
    have = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                have.add(line[1:].split()[0])
    return have


def fetch(ids, email, api_key, tries=5):
    q = {"db": "nuccore", "rettype": "fasta", "retmode": "text",
         "tool": TOOL, "email": email, "id": ",".join(ids)}
    if api_key:
        q["api_key"] = api_key
    url = f"{EUTILS}?{urllib.parse.urlencode(q)}"
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                text = r.read().decode("utf-8", "replace")
            if text.count(">") == 0:
                raise ValueError("empty response")
            return text
        except (urllib.error.URLError, urllib.error.HTTPError, ValueError, OSError) as e:
            wait = min(60, 2 ** attempt)
            print(f"    retry {attempt}/{tries} after {type(e).__name__}: {e} (sleep {wait}s)",
                  file=sys.stderr, flush=True)
            time.sleep(wait)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?", default="rivm_test_set.fasta")
    ap.add_argument("--batch", type=int, default=50, help="accessions per request (keeps the URL short)")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", "calicityper@example.org"))
    ap.add_argument("--all", action="store_true", dest="every",
                    help="fetch all 5,986 RIVM records, not just the 4,392 with a type call "
                         "(needed to score genogroup detection on non-norovirus caliciviruses)")
    a = ap.parse_args()

    out = Path(a.out)
    accs = wanted_accessions(a.every)
    cut = drop_unterminated_tail(out)
    if cut:
        print(f"dropped unterminated last record {cut} (interrupted run); refetching it")
    have = already_have(out)
    todo = [x for x in accs if x not in have]
    delay = 0.12 if a.api_key else 0.40           # stay under 10/s with a key, 3/s without

    print(f"RIVM-called accessions : {len(accs)}")
    print(f"already in {out.name:22s}: {len(have)}")
    print(f"to fetch               : {len(todo)}  "
          f"({-(-len(todo) // a.batch)} requests, api_key={'yes' if a.api_key else 'no'})")
    if not todo:
        print("nothing to do.")
        return 0

    failed = []
    with open(out, "a") as fh:
        for i in range(0, len(todo), a.batch):
            chunk = todo[i:i + a.batch]
            text = fetch(chunk, a.email, a.api_key)
            if text is None:
                failed += chunk
                print(f"  [{i // a.batch + 1}] FAILED {len(chunk)} accessions", flush=True)
                continue
            fh.write(text if text.endswith("\n") else text + "\n")
            fh.flush()
            got = text.count(">")
            print(f"  [{i // a.batch + 1:3d}] {got:3d}/{len(chunk)} seqs  "
                  f"({i + len(chunk)}/{len(todo)})", flush=True)
            time.sleep(delay)

    final = already_have(out)
    missing = [x for x in accs if x not in final]
    print(f"\nwrote {out}  sequences={len(final)}")
    if missing:
        Path("missing_accessions.txt").write_text("\n".join(missing) + "\n")
        print(f"MISSING {len(missing)} accessions (listed in missing_accessions.txt); "
              f"rerun this script to retry them")
    if failed:
        print(f"{len(failed)} accessions were in batches that failed outright")
    return 0


if __name__ == "__main__":
    sys.exit(main())
