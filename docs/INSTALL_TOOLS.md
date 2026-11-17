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

UCD's own Apptainer guidance (ServiceNow KB0011720) gives the build command as:

```bash
export APPTAINER_CACHEDIR=/scratch/$USER/APPTAINER_CACHE
export APPTAINER_TMPDIR=/tmp
apptainer build --fakeroot --sandbox --bind="$TMPDIR:/tmp" out.sif in.def
```

Two notes on applying that here.

**`--fakeroot` does not work for this account.** `fakeroot` is not installed on
Sonic and the account has no `/etc/subuid` or `/etc/subgid` range, so Apptainer
falls back to a root-mapped namespace. That is worth raising with the Sonic
administrators, since their own documentation prescribes the flag: ask for a
subuid range, or for `fakeroot` to be installed.

**It is not needed for this image.** `--fakeroot` matters when `%post` must act
as root, and the usual reason is a package manager. In a root-mapped namespace
`apt-get` cannot drop privileges to the `_apt` user and dies with
`E: setgroups 65534 failed - setgroups (1: Operation not permitted)`, which is
what the first attempt at this image hit. Everything else in `%post` works
normally. `containers/calicityper-tools.def` therefore uses no package manager
at all: it starts from an image that already carries Python, and fetches BLAST+
as its prebuilt tarball and MAFFT as its portable binary package, both unpacked
with Python's standard library. No compilation, no apt, no privilege.

**`--sandbox` is the fallback if the SIF step fails.** Building a SIF ends with
`mksquashfs`, which on the login node fails with `Failed to create thread`
against the process limit there. `--sandbox` writes a writable directory tree
instead and never calls `mksquashfs`; the tools are usable directly from it, and
it converts to a SIF afterwards:

```bash
apptainer build --sandbox tools.sandbox containers/calicityper-tools.def
apptainer exec tools.sandbox blastn -version      # usable as-is
apptainer build tools.sif tools.sandbox           # convert when convenient
```

Build on a compute node rather than the login node either way;
`calicityper_runs/scripts/build_tools_sif.sh` does this, with the sandbox
fallback wired in.

**Prebuilt images need no privilege at all.** `apptainer pull` is the simplest
route where a suitable image exists:

```bash
apptainer pull blast.sif docker://quay.io/biocontainers/blast:2.17.0--h66d330f_0
apptainer pull mafft.sif docker://quay.io/biocontainers/mafft:7.526--h4bc722e_0
```

There is no OpenRDP biocontainer, so OpenRDP still needs the virtual
environment below. Check the tags before relying on them; they are not pinned
here because they have not been tested on this system.

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
