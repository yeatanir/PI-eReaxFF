#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

if [[ ! -x ".venv-fast/bin/python" ]]; then
  echo "test_33_fast_jax_stability: missing .venv-fast/bin/python; run the JAX install step first" >&2
  exit 77
fi

JAMIL_FAST_STABILITY_ROOT="$ROOT_DIR/build/test_33_fast_jax_stability" \
JAMIL_FAST_STABILITY_DT_CASES="0.05" \
JAMIL_FAST_STABILITY_STEP_CASES="10" \
projects/jamil_replication/run_fast_stability_check.sh

grep -q "| 16 | 10 | 0.050" projects/jamil_replication/fast_stability_summary.out
grep -q "| PASS |" projects/jamil_replication/fast_stability_summary.out
test -s build/test_33_fast_jax_stability/fast_p16_steps10_dt0p05/out/fast_jax_work/worker_timing.csv

echo "test_33_fast_jax_stability: PASS"
