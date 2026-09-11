#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

SOURCE_RUN="$ROOT_DIR/build/test_33_fast_jax_stability/fast_p16_steps10_dt0p05/out"
if [[ ! -d "$SOURCE_RUN" ]]; then
  make test_33_fast_jax_stability
fi

JAMIL_FAST_PARITY_SOURCE_RUN="$SOURCE_RUN" \
JAMIL_FAST_PARITY_ROOT="$ROOT_DIR/build/test_34_fast_offequilibrium_parity" \
JAMIL_FAST_PARITY_STEPS="0 10" \
JAMIL_FAST_PARITY_BEADS="1" \
JAMIL_FAST_PARITY_MAX_FORCE_ERROR="1e-5" \
projects/jamil_replication/run_fast_offequilibrium_parity.sh

grep -q "cases: \`2\`" projects/jamil_replication/fast_offequilibrium_parity_summary.out
test -s build/test_34_fast_offequilibrium_parity/offequilibrium_oracle_cases.jsonl
test -s build/test_34_fast_offequilibrium_parity/offequilibrium_force_report.json

echo "test_34_fast_offequilibrium_parity: PASS"
