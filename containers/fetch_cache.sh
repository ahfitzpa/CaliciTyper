#!/usr/bin/env bash
# Populate containers/cache/ so an Apptainer build needs no network.
#
# Each file is optional: calicityper-tools.def uses the cache when a file is
# there and downloads when it is not. Run this while the upstream hosts are
# reachable, and the image can then be rebuilt whenever one of them is not.
# mafft.cbrc.jp returned 503 to everything for days, which is the case this
# exists for.
set -uo pipefail
cd "$(dirname "$0")/cache" || exit 1

get () {  # name url
    if [ -s "$1" ]; then echo "  have   $1"; return 0; fi
    if curl -fsSL --max-time 600 -o "$1.part" "$2"; then
        mv "$1.part" "$1"; echo "  got    $1 ($(du -h "$1" | cut -f1))"
    else
        rm -f "$1.part"; echo "  FAILED $1  <- $2"
    fi
}

get blast-2.17.0-x64-linux.tar.gz \
    "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz"

# Preferred: the author's prebuilt portable package. If its host is down the
# build falls back to the source below, which is the same version from the
# author's own GitLab repository and is compiled in the image.
get mafft-7.526-linux.tgz \
    "https://mafft.cbrc.jp/alignment/software/mafft-7.526-linux.tgz"
get mafft-v7.526-src.tar.gz \
    "https://gitlab.com/sysimm/mafft/-/archive/v7.526/mafft-v7.526.tar.gz"

# b42d9e5 (master), not the v0.1.0 tag: the tag does not cap numpy and so
# resolves numpy 2.x, defeating the point of isolating OpenRDP.
get openrdp-b42d9e5.tar.gz \
    "https://github.com/PoonLab/OpenRDP/archive/b42d9e5fb84c14e817e08aab2668aeb2c6a94f83.tar.gz"

echo
sha256sum ./*.tar.gz ./*.tgz 2>/dev/null | tee SHA256SUMS
