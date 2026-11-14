#!/usr/bin/env bash
# Thin wrapper kept for backwards compatibility: the real implementation is
# fetch_rivm_test_set.py (resumable, and it does not use entrez-direct, whose
# batched POSTs are dropped by NCBI under throttling).
set -euo pipefail
exec python3 "$(dirname "$0")/fetch_rivm_test_set.py" "$@"
