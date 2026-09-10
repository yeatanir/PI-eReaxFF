#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/projects/jamil_replication"
BASE_CFG="$SCRIPT_DIR/config_jamil_112c_113el_p4_pool_verify.yaml"
SWEEP_ROOT="${JAMIL_P_SWEEP_ROOT:-$ROOT_DIR/build/jamil_replication/p_sweep_$(date +%Y%m%d_%H%M%S)}"
SUMMARY="$SCRIPT_DIR/jamil_p_sweep_summary.out"
WORKERS="${JAMIL_P_SWEEP_WORKERS:-4}"
OMP_THREADS="${JAMIL_P_SWEEP_OMP_THREADS:-$WORKERS}"
STEPS="${JAMIL_P_SWEEP_STEPS:-1}"
DT_FS="${JAMIL_P_SWEEP_DT_FS:-0.1}"

if [[ ! -f "$BASE_CFG" ]]; then
  echo "run_jamil_p_sweep: missing base config: $BASE_CFG" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/build/ereaxff_mixed_baoab_main" ]]; then
  echo "run_jamil_p_sweep: missing executable: $ROOT_DIR/build/ereaxff_mixed_baoab_main" >&2
  exit 1
fi

mkdir -p "$SWEEP_ROOT"

if [[ -n "${JAMIL_P_SWEEP_CASES:-}" ]]; then
  read -r -a case_specs <<<"$JAMIL_P_SWEEP_CASES"
else
  case_specs=(4 8 16)
fi

extract_last_column() {
  local file="$1"
  local col="$2"
  awk -v col="$col" 'NF && $1 !~ /^#/ {v=$col} END {print v}' "$file"
}

summarize_case() {
  local case_name="$1"
  local p="$2"
  local case_dir="$3"
  local result_json="$case_dir/result.json"
  local out="$case_dir/out"
  local timing="$out/adapter_work/persistent_oracle_timing.log"

  python3 - "$case_name" "$p" "$WORKERS" "$OMP_THREADS" "$result_json" "$out" "$timing" "$SUMMARY" <<'PY'
from __future__ import annotations

import json
import math
import pathlib
import statistics
import sys

case_name = sys.argv[1]
p = int(sys.argv[2])
workers = int(sys.argv[3])
omp_threads = int(sys.argv[4])
result_path = pathlib.Path(sys.argv[5])
out = pathlib.Path(sys.argv[6])
timing_path = pathlib.Path(sys.argv[7])
summary_path = pathlib.Path(sys.argv[8])


def last_data(path: pathlib.Path, col: int) -> float:
    value = math.nan
    if not path.exists():
        return value
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split()
        if not parts or parts[0].startswith("#"):
            continue
        if len(parts) >= col:
            value = float(parts[col - 1])
    return value


def timing_stats(path: pathlib.Path) -> tuple[int, int, float, float, float, float]:
    calls = 0
    batches = 0
    batch_elapsed: list[float] = []
    if not path.exists():
        return calls, batches, math.nan, math.nan, math.nan, math.nan
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "REQ":
            calls += 1
        elif parts[0] == "BATCH":
            batches += 1
            batch_elapsed.append(float(parts[8]))
    if not batch_elapsed:
        return calls, batches, math.nan, math.nan, math.nan, math.nan
    return (
        calls,
        batches,
        min(batch_elapsed),
        statistics.fmean(batch_elapsed),
        max(batch_elapsed),
        sum(batch_elapsed),
    )


result = json.loads(result_path.read_text(encoding="utf-8")) if result_path.exists() else {}
status = "PASS" if result.get("status") == 0 else f"FAIL:{result.get('status', 'missing')}"
calls, batches, min_batch, mean_batch, max_batch, sum_batch = timing_stats(timing_path)

row = {
    "case": case_name,
    "P": p,
    "workers": workers,
    "omp_threads": omp_threads,
    "wall_s": float(result.get("wall_s", math.nan)),
    "oracle_calls": calls,
    "batches": batches,
    "batch_sum_s": sum_batch,
    "batch_mean_s": mean_batch,
    "batch_min_s": min_batch,
    "batch_max_s": max_batch,
    "electron_T_K": last_data(out / "electron_temperature.out", 2),
    "electron_extended_T_K": last_data(out / "electron_temperature.out", 3),
    "electron_expected_extended_T_K": last_data(out / "electron_temperature.out", 4),
    "electron_ratio": last_data(out / "electron_temperature.out", 5),
    "carbon_T_K": last_data(out / "carbon_temperature.out", 2),
    "thermostat_heat_eV": last_data(out / "thermostat_heat.out", 2),
    "quantum_heat_eV": last_data(out / "thermostat_heat.out", 4),
    "bead_rms_spread_A": last_data(out / "bead_spread.out", 5),
    "status": status,
}

summary_path.parent.mkdir(parents=True, exist_ok=True)
exists = summary_path.exists()
with summary_path.open("a", encoding="utf-8") as handle:
    if not exists:
        handle.write("# Jamil P Sweep Summary\n\n")
        handle.write(
            "| case | P | workers | OMP | wall_s | oracle_calls | batches | batch_sum_s | "
            "batch_mean_s | electron_T_K | extended_T_K | expected_extended_T_K | ratio | "
            "carbon_T_K | thermostat_heat_eV | quantum_heat_eV | bead_rms_spread_A | status |\n"
        )
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |\n")
    handle.write(
        f"| {row['case']} | {row['P']} | {row['workers']} | {row['omp_threads']} | "
        f"{row['wall_s']:.3f} | {row['oracle_calls']} | {row['batches']} | {row['batch_sum_s']:.3f} | "
        f"{row['batch_mean_s']:.3f} | {row['electron_T_K']:.3f} | {row['electron_extended_T_K']:.3f} | "
        f"{row['electron_expected_extended_T_K']:.3f} | {row['electron_ratio']:.3f} | "
        f"{row['carbon_T_K']:.3f} | {row['thermostat_heat_eV']:.6f} | {row['quantum_heat_eV']:.6f} | "
        f"{row['bead_rms_spread_A']:.6f} | {row['status']} |\n"
    )
print(json.dumps(row, sort_keys=True))
PY
}

rm -f "$SUMMARY"

for p in "${case_specs[@]}"; do
  if [[ ! "$p" =~ ^[0-9]+$ ]]; then
    echo "run_jamil_p_sweep: P value must be an integer, got '$p'" >&2
    exit 2
  fi
  case_name="p${p}_w${WORKERS}_steps${STEPS}_dt${DT_FS//./p}"
  case_dir="$SWEEP_ROOT/$case_name"
  cfg="$case_dir/config.yaml"
  out="$case_dir/out"
  stdout_log="$case_dir/stdout.log"
  stderr_log="$case_dir/stderr.log"
  result_json="$case_dir/result.json"
  mkdir -p "$case_dir" "$out"

  python3 - "$BASE_CFG" "$cfg" "$p" "$WORKERS" "$STEPS" "$DT_FS" "$out" <<'PY'
from __future__ import annotations

import pathlib
import re
import sys

src = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
dst = pathlib.Path(sys.argv[2])
p = sys.argv[3]
workers = sys.argv[4]
steps = sys.argv[5]
dt_fs = sys.argv[6]
out = sys.argv[7]

def replace(pattern: str, value: str, text: str) -> str:
    return re.sub(pattern, value, text, flags=re.M)

text = src
text = replace(r"^  dt_fs:.*$", f"  dt_fs: {dt_fs}", text)
text = replace(r"^  steps:.*$", f"  steps: {steps}", text)
text = replace(r"^  output_dir:.*$", f'  output_dir: "{out}"', text)
text = replace(r"^  oracle_backend_mode:.*$", "  oracle_backend_mode: persistent_worker_pool", text)
text = replace(r"^  oracle_dispatch_mode:.*$", "  oracle_dispatch_mode: bounded_pool", text)
text = replace(r"^  oracle_max_workers:.*$", f"  oracle_max_workers: {workers}", text)
text = replace(r"^  oracle_timing_log:.*$", "  oracle_timing_log: true", text)
text = replace(r"^  oracle_timeout_s:.*$", "  oracle_timeout_s: 180", text)
text = replace(r"^  repeatability_check:.*$", "  repeatability_check: false", text)
text = replace(r"^    beads:.*$", f"    beads: {p}", text)
dst.write_text(text, encoding="utf-8")
PY

  python3 - "$ROOT_DIR" "$cfg" "$out" "$OMP_THREADS" "$stdout_log" "$stderr_log" "$result_json" <<'PY'
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys
import time

root = pathlib.Path(sys.argv[1])
cfg = sys.argv[2]
out = sys.argv[3]
omp_threads = sys.argv[4]
stdout_log = pathlib.Path(sys.argv[5])
stderr_log = pathlib.Path(sys.argv[6])
result_json = pathlib.Path(sys.argv[7])

env = os.environ.copy()
env["OMP_NUM_THREADS"] = omp_threads
cmd = [str(root / "build/ereaxff_mixed_baoab_main"), cfg, out]
start = time.time()
status = 0
note = "ok"
try:
    with stdout_log.open("w", encoding="utf-8") as stdout, stderr_log.open("w", encoding="utf-8") as stderr:
        proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env)
    status = proc.returncode
    if status != 0:
        note = f"exit_{status}"
except Exception as exc:
    status = 125
    note = f"error:{exc.__class__.__name__}"
end = time.time()
result_json.write_text(
    json.dumps(
        {
            "command": cmd,
            "omp_threads": int(omp_threads),
            "status": status,
            "note": note,
            "wall_s": round(end - start, 6),
            "config": cfg,
            "out_dir": out,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
print(note)
PY

  summarize_case "$case_name" "$p" "$case_dir"
done

echo "sweep_root=$SWEEP_ROOT"
echo "summary=$SUMMARY"
