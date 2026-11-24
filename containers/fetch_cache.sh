#!/usr/bin/env bash
# Optional offline cache for containers/calicityper-tools.def (falls back to download if absent).
set -uo pipefail
CACHE_DIR="$(dirname "$0")/cache"
mkdir -p "$CACHE_DIR"
cd "$CACHE_DIR" || exit 1

get () {
    if [ -s "$1" ]; then echo "  have   $1"; return 0; fi
    if curl -fsSL --max-time 600 -o "$1.part" "$2"; then
        mv "$1.part" "$1"; echo "  got    $1 ($(du -h "$1" | cut -f1))"
    else
        rm -f "$1.part"; echo "  FAILED $1  <- $2"
    fi
}

get micromamba-2.0.5-linux-64 \
    "https://github.com/mamba-org/micromamba-releases/releases/download/2.0.5-0/micromamba-linux-64"

# Commit b42d9e5 caps numpy<2; the v0.1.0 tag does not.
get openrdp-b42d9e5.tar.gz \
    "https://github.com/PoonLab/OpenRDP/archive/b42d9e5fb84c14e817e08aab2668aeb2c6a94f83.tar.gz"

echo
sha256sum ./micromamba-2.0.5-linux-64 ./openrdp-b42d9e5.tar.gz 2>/dev/null | tee SHA256SUMS
