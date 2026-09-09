#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
WORK_DIR="$ROOT_DIR/build/test_29_oracle_dispatch_modes"
BIN="$ROOT_DIR/build/test_29_oracle_dispatch_modes_bin"
SUMMARY="$SUMMARY_DIR/test_29_oracle_dispatch_modes.summary"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

if [[ ! -x "$BIN" ]]; then
  {
    echo "test_29_oracle_dispatch_modes: FAIL"
    echo "Missing executable: $BIN"
  } > "$SUMMARY"
  cat "$SUMMARY"
  exit 1
fi

export OMP_NUM_THREADS=1
export EREAXFF_ADAPTER_WORK_ROOT="$WORK_DIR"
export EREAXFF_ORACLE_SMOKE_BINARY="$ROOT_DIR/build/oracle_single_point_force_export"

"$BIN" > "$WORK_DIR/stdout.log"

python3 - "$WORK_DIR" <<'PY'
from __future__ import annotations

import pathlib
import sys

root = pathlib.Path(sys.argv[1])
expected = {
    "serial": {"batches": 2, "jobs": 2, "batch_sizes": [1, 1]},
    "bounded_pool": {"batches": 1, "jobs": 2, "batch_sizes": [2]},
    "openmp_all": {"batches": 1, "jobs": 2, "batch_sizes": [2]},
}

for mode, spec in expected.items():
    log = root / mode / "oracle_timing.log"
    if not log.exists():
        raise SystemExit(f"FAIL: missing timing log for {mode}")
    lines = [line.strip() for line in log.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    batch_lines = [line for line in lines if line.startswith("BATCH")]
    call_lines = [line for line in lines if line.startswith("CALL")]
    if len(batch_lines) != spec["batches"]:
        raise SystemExit(f"FAIL: expected {spec['batches']} batch records for {mode}, got {len(batch_lines)}")
    if len(call_lines) != spec["jobs"]:
        raise SystemExit(f"FAIL: expected {spec['jobs']} call records for {mode}, got {len(call_lines)}")
    batch_sizes = []
    beads = []
    batches = []
    for line in batch_lines:
        parts = line.split()
        if len(parts) < 13:
            raise SystemExit(f"FAIL: malformed BATCH record for {mode}: {line}")
        batch_sizes.append(int(parts[5]))
    if batch_sizes != spec["batch_sizes"]:
        raise SystemExit(f"FAIL: unexpected batch sizes for {mode}: {batch_sizes}")
    for line in call_lines:
        parts = line.split()
        if len(parts) < 13:
            raise SystemExit(f"FAIL: malformed CALL record for {mode}: {line}")
        batches.append(int(parts[2]))
        beads.append(int(parts[4]))
    if beads != list(range(1, spec["jobs"] + 1)):
        raise SystemExit(f"FAIL: unexpected bead order for {mode}: {beads}")
    if mode == "serial":
        if batches != [1, 2]:
            raise SystemExit(f"FAIL: unexpected batch numbers for {mode}: {batches}")
    else:
        if any(batch != 1 for batch in batches):
            raise SystemExit(f"FAIL: unexpected batch numbers for {mode}: {batches}")

stdout = (root / "stdout.log").read_text(encoding="utf-8")
if "test_29_oracle_dispatch_modes: PASS" not in stdout:
    raise SystemExit("FAIL: Fortran test did not report PASS")
PY

{
  echo "test_29_oracle_dispatch_modes: PASS"
  echo "  modes = serial bounded_pool openmp_all"
  echo "  p = 2"
  echo "  timing_log = oracle_timing.log"
} | tee "$SUMMARY"
