# External tools: BLAST+, MAFFT and OpenRDP

CaliciTyper shells out to three external programs. This is how to install them
reproducibly, and why OpenRDP is not simply a dependency of the Python package.

| Tool | Used for | Found via |
|---|---|---|
| BLAST+ (`blastn`, `makeblastdb`) | genogroup assignment | `$CALICITYPER_BLASTN`, `$CALICITYPER_MAKEBLASTDB`, else `PATH` |
| MAFFT | adding queries to reference alignments | `$CALICITYPER_MAFFT`, else `PATH` |
| OpenRDP | recombination scanners | `$CALICITYPER_OPENRDP`, else `PATH` |

BLAST+ and MAFFT are required for typing. OpenRDP is needed only for
`calicityper recomb`.

## Why OpenRDP is not in the main environment

It cannot be. OpenRDP 0.1.0 requires **numpy < 2**; the current environment
runs numpy 2.5. Installing OpenRDP alongside CaliciTyper forces a numpy
downgrade across the whole environment, which changes the numerical stack the
typing path runs on. That is a bad trade for one optional subcommand.

So OpenRDP lives in its own environment, and CaliciTyper is told where to find
it. The consequence to be aware of: **when OpenRDP is absent, 13 tests -- every
test that exercises a real scanner -- are silently deselected.** A green test
run does not mean the recombination path was tested. Check for
`deselected` in the pytest summary.

## Recommended: one Apptainer image, version-pinned

Best option on an HPC system, and not conda. A single SIF pins all three tools
and their numerical stacks at exact versions, is a single file you can
checksum and archive beside your results, and sidesteps the numpy conflict by
keeping OpenRDP's interpreter separate from CaliciTyper's.

Apptainer is available on Sonic as a module:

```bash
module load apptainer/1.3.4-gcc-11.5.0-ojp6nts
```

`containers/calicityper-tools.def` builds the image:

```bash
apptainer build calicityper-tools.sif containers/calicityper-tools.def
sha256sum calicityper-tools.sif          # record this beside your results
```

Point CaliciTyper at the containerised binaries with small wrappers:

```bash
SIF=/path/to/calicityper-tools.sif
mkdir -p ~/bin/ct-tools
for t in blastn makeblastdb mafft openrdp; do
  printf '#!/bin/sh\nexec apptainer exec %s %s "$@"\n' "$SIF" "$t" > ~/bin/ct-tools/$t
  chmod +x ~/bin/ct-tools/$t
done
export CALICITYPER_BLASTN=~/bin/ct-tools/blastn
export CALICITYPER_MAKEBLASTDB=~/bin/ct-tools/makeblastdb
export CALICITYPER_MAFFT=~/bin/ct-tools/mafft
export CALICITYPER_OPENRDP=~/bin/ct-tools/openrdp
```

Bind any path outside your home directory that the run must read or write,
for example `apptainer exec -B /scratch "$SIF" ...`. On a cluster, build the
image once on a machine where you have root (or use `--fakeroot`), then copy
the SIF.

## Alternative: a separate virtual environment

What is in use now, and perfectly adequate. Not conda, no container runtime.

```bash
python3.11 -m venv ~/scratch/openrdp_env
~/scratch/openrdp_env/bin/pip install 'numpy<2' openrdp==0.1.0
export CALICITYPER_OPENRDP=~/scratch/openrdp_env/bin/openrdp
```

Set that variable in every SLURM script that calls `recomb`, not only in your
interactive shell. Omitting it does not fail loudly: the run completes and
reports `not_assessed`, which is honest but easy to miss in a batch of jobs.

`pipx install openrdp` or `uv tool install openrdp` give the same isolation
with less ceremony, but neither pins BLAST+ or MAFFT, so neither helps with
version control of the full toolchain.

## Pinned versions

The versions these results were produced with:

| Tool | Version |
|---|---|
| BLAST+ | 2.17.0+ |
| MAFFT | 7.526 (2024/Apr/26) |
| OpenRDP | 0.1.0 (numpy 1.26.4) |
| CaliciTyper | 0.8.0 (numpy 2.5.5) |

Record these, or the SIF checksum, alongside any results you report. BLAST
genogroup assignment also depends on the database search space, which
CaliciTyper pins explicitly so that a call cannot change as the reference
database grows.

## Checking what is available

```bash
calicityper recomb --dry-run --input your.fasta --query ID
```

reports which tools were found and which would be skipped, without running
anything.
