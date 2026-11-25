# External tools: BLAST+, MAFFT, IQ-TREE, ClustalW and OpenRDP

CaliciTyper invokes external bioinformatics tools depending on the subcommand. This guide covers reproducible installation and explains why OpenRDP is kept outside the main Python environment.

| Tool | Used for | Required by | Found via |
|---|---|---|---|
| BLAST+ (`blastn`, `makeblastdb`) | genogroup assignment | `type`, `novel` | `$CALICITYPER_BLASTN`, `$CALICITYPER_MAKEBLASTDB`, else `PATH` |
| MAFFT | profile alignment onto fixed reference columns | `type`, `novel`, `recomb --input`, `antigenic --unaligned` | `$CALICITYPER_MAFFT`, else `PATH` |
| IQ-TREE (`iqtree`) | ML patristic trees for `mean + 2 SD` novelty bounds (NJ fallback if absent) | `novel` (optional) | `$CALICITYPER_IQTREE`, else `PATH` |
| ClustalW (`clustalw2`) | legacy `--backend clustalw` alignment path | `type --backend clustalw` (optional) | `$CALICITYPER_CLUSTALW`, else `PATH` |
| OpenRDP (`openrdp`) | recombination scanners | `recomb` | `$CALICITYPER_OPENRDP`, else `PATH` |

## Why OpenRDP is not in the main environment

OpenRDP requires `numpy<2`, whereas CaliciTyper runs on `numpy>=2.0,<3`. Installing OpenRDP alongside CaliciTyper forces a NumPy downgrade across the whole typing stack. Keeping OpenRDP in an isolated virtual environment (or inside the Apptainer container's separate venv) avoids that conflict.

Note for test runs: when `openrdp` is absent from `PATH` and `CALICITYPER_OPENRDP` is unset, `pytest tests/` deselects the 13 tests that exercise live OpenRDP scanners. Check for `deselected` in the `pytest` summary if you need to verify the recombination path.

## Recommended for HPC: Apptainer containers

Two layered Apptainer definitions live under `containers/` and build without `apt-get` or `--fakeroot` (working cleanly inside unprivileged root-mapped user namespaces on HPC clusters):

1. `containers/calicityper-tools.def` — starts from a digest-pinned `python:3.11` base image, installs the locked conda stack (`BLAST+`, `MAFFT`, `ClustalW`, `IQ-TREE`, Python `3.12`, NumPy `>=2`, Numba, Biopython, Pandas, Bokeh) via `micromamba` from `environment.lock.yml`, and installs OpenRDP (`PoonLab/OpenRDP` commit `b42d9e5`, `numpy 1.26.4`) into a separate `/opt/openrdp` virtual environment.
2. `containers/calicityper.def` — layers `calicityper` offline on top of `calicityper-tools.sif`, pre-warms Numba caches, and sets `calicityper` as the container entrypoint.

### Build both images (from the repository root)

```bash
# Optional: pre-fetch micromamba and OpenRDP into containers/cache/ for offline builds
bash containers/fetch_cache.sh

apptainer build calicityper-tools.sif containers/calicityper-tools.def
apptainer build calicityper.sif containers/calicityper.def
apptainer test calicityper.sif
sha256sum calicityper.sif    # record alongside reported results
```

On SLURM clusters, build on a compute node rather than a login node (where `mksquashfs` may hit per-user thread limits). A SLURM wrapper is provided in `containers/build_images.sh`. If `mksquashfs` fails in your environment, build a writable `--sandbox` directory first and convert it afterwards:

```bash
apptainer build --sandbox calicityper-tools.sandbox containers/calicityper-tools.def
apptainer build calicityper-tools.sif calicityper-tools.sandbox
```

### Run `calicityper.sif`

Reference data (`doi:10.5281/zenodo.24354774`) is kept outside the image and bind-mounted to `/data` (or pointed to via `CALICITYPER_DATA`):

```bash
REAL_DIR=$(pwd -P)
mkdir -p "$REAL_DIR/refdata"

# Download and verify the manifest-matched reference archive once
apptainer exec -B "$REAL_DIR/refdata:/data" calicityper.sif \
    calicityper setup --data-dir /data

# Run typing (resolve symlinks with pwd -P on HPC scratch mounts)
apptainer run -B "$REAL_DIR/refdata:/data" -B "$REAL_DIR:/work" calicityper.sif \
    type /work/examples/example.fasta -o /work/example_typed.tsv --threads 4
```

## Alternative: standalone OpenRDP virtual environment

If you run CaliciTyper in a standard Conda environment (`environment.yml` or `environment.lock.yml`) and want `calicityper recomb`, install OpenRDP in a separate Python 3.11 virtual environment:

```bash
python3.11 -m venv ~/openrdp_env
~/openrdp_env/bin/pip install 'numpy<2' \
    https://github.com/PoonLab/OpenRDP/archive/b42d9e5fb84c14e817e08aab2668aeb2c6a94f83.tar.gz
export CALICITYPER_OPENRDP=~/openrdp_env/bin/openrdp
```

Export `CALICITYPER_OPENRDP` in batch/SLURM scripts as well as interactive shells. When OpenRDP is missing, `calicityper recomb` records skipped scanners and reports `not_assessed` rather than crashing.

## Tested pinned versions

The container image and manuscript benchmarks were verified with:

| Tool | Version |
|---|---|
| CaliciTyper | 0.8.1 |
| BLAST+ | 2.17.0+ |
| MAFFT | 7.526 (2024/Apr/26) |
| IQ-TREE | 3.1.4 |
| ClustalW | 2.1 |
| OpenRDP | 0.1.0 (`PoonLab/OpenRDP@b42d9e5`, `numpy 1.26.4`) |

Record these versions (or the `calicityper.sif` SHA-256 digest and each run's `*.provenance.json` sidecar) alongside reported results.

## Checking tool availability

Verify which recombination scanners and alignment tools are detected without running an analysis:

```bash
calicityper recomb --dry-run --input examples/example.fasta --query MW305639.1
```
