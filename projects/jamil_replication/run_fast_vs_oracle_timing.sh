#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/projects/jamil_replication"
BASE_CFG="$SCRIPT_DIR/config_jamil_112c_113el_p4_pool_verify.yaml"
RUN_ROOT="${JAMIL_FAST_TIMING_ROOT:-$ROOT_DIR/build/jamil_replication/fast_vs_oracle_timing_$(date +%Y%m%d_%H%M%S)}"
SUMMARY="$SCRIPT_DIR/fast_vs_oracle_timing_summary.out"
STEPS="${JAMIL_FAST_TIMING_STEPS:-1}"
DT_FS="${JAMIL_FAST_TIMING_DT_FS:-0.1}"
OMP_THREADS="${JAMIL_FAST_TIMING_OMP_THREADS:-4}"
ORACLE_WORKERS="${JAMIL_FAST_TIMING_ORACLE_WORKERS:-4}"
if [[ ${JAMIL_FAST_TIMING_FAST_CASES+x} ]]; then
  FAST_CASES="$JAMIL_FAST_TIMING_FAST_CASES"
else
  FAST_CASES="4 8"
fi
if [[ ${JAMIL_FAST_TIMING_ORACLE_CASES+x} ]]; then
  ORACLE_CASES="$JAMIL_FAST_TIMING_ORACLE_CASES"
else
  ORACLE_CASES="4"
fi

if [[ ! -f "$BASE_CFG" ]]; then
  echo "run_fast_vs_oracle_timing: missing base config: $BASE_CFG" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/build/ereaxff_mixed_baoab_main" ]]; then
  echo "run_fast_vs_oracle_timing: missing executable; run make ereaxff_mixed_baoab_main" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/.venv-fast/bin/python" ]]; then
  echo "run_fast_vs_oracle_timing: missing .venv-fast/bin/python" >&2
  exit 77
fi

mkdir -p "$RUN_ROOT"
rm -f "$SUMMARY"

write_config() {
  local backend="$1"
  local p="$2"
  local out="$3"
  local cfg="$4"

  python3 - "$BASE_CFG" "$cfg" "$backend" "$p" "$STEPS" "$DT_FS" "$out" "$ORACLE_WORKERS" <<'PY'
from __future__ import annotations

import pathlib
import re
import sys

src = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
dst = pathlib.Path(sys.argv[2])
backend = sys.argv[3]
p = sys.argv[4]
steps = sys.argv[5]
dt_fs = sys.argv[6]
out = sys.argv[7]
workers = sys.argv[8]


def replace(pattern: str, value: str, text: str) -> str:
    return re.sub(pattern, value, text, flags=re.M)


text = src
text = replace(r"^  dt_fs:.*$", f"  dt_fs: {dt_fs}", text)
text = replace(r"^  steps:.*$", f"  steps: {steps}", text)
text = replace(r"^  output_stride:.*$", "  output_stride: 1", text)
if "rescale_initial_pimd_velocities:" in text:
    text = replace(r"^  rescale_initial_pimd_velocities:.*$", "  rescale_initial_pimd_velocities: true", text)
else:
    text = replace(r"^  restart_stride:.*$", "  restart_stride: 1\n  rescale_initial_pimd_velocities: true", text)
text = replace(r"^  output_dir:.*$", f'  output_dir: "{out}"', text)
text = replace(r"^    beads:.*$", f"    beads: {p}", text)

if backend == "fast":
    text = replace(r"^  type:.*$", "  type: fast_jax_oracle", text)
    insertion = (
        "  repeatability_check: false\n"
        '  fast_backend_command: ".venv-fast/bin/python projects/jamil_replication/fast_backend/jax_live_worker.py"\n'
        '  fast_backend_model: "reaxff_corrected_bond_nonbonded_valence"\n'
        "  fast_backend_allow_unvalidated: true"
    )
    text = replace(r"^  repeatability_check:.*$", insertion, text)
elif backend == "oracle":
    text = replace(r"^  type:.*$", "  type: ereaxff_oracle", text)
    text = replace(r"^  oracle_backend_mode:.*$", "  oracle_backend_mode: persistent_worker_pool", text)
    text = replace(r"^  oracle_dispatch_mode:.*$", "  oracle_dispatch_mode: bounded_pool", text)
    text = replace(r"^  oracle_max_workers:.*$", f"  oracle_max_workers: {workers}", text)
    text = replace(r"^  oracle_timing_log:.*$", "  oracle_timing_log: true", text)
    text = replace(r"^  repeatability_check:.*$", "  repeatability_check: false", text)
else:
    raise SystemExit(f"unknown backend: {backend}")

dst.write_text(text, encoding="utf-8")
PY
}

summarize_case() {
  local backend="$1"
  local p="$2"
  local case_dir="$3"
  local out="$case_dir/out"
  local result_json="$case_dir/result.json"
  local timing="$out/adapter_work/persistent_oracle_timing.log"

  python3 - "$backend" "$p" "$STEPS" "$OMP_THREADS" "$ORACLE_WORKERS" "$result_json" "$out" "$timing" "$SUMMARY" <<'PY'
from __future__ import annotations

import json
import math
import pathlib
import statistics
import sys

backend = sys.argv[1]
p = int(sys.argv[2])
steps = int(sys.argv[3])
omp = int(sys.argv[4])
workers = int(sys.argv[5])
result_path = pathlib.Path(sys.argv[6])
out = pathlib.Path(sys.argv[7])
timing = pathlib.Path(sys.argv[8])
summary = pathlib.Path(sys.argv[9])


def last_data(path: pathlib.Path, col: int) -> float:
    value = math.nan
    if not path.exists():
        return value
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split()
        if parts and not parts[0].startswith("#") and len(parts) >= col:
            value = float(parts[col - 1])
    return value


def oracle_timing(path: pathlib.Path) -> tuple[int, int, float, float]:
    calls = 0
    batches = 0
    elapsed = []
    if not path.exists():
        return 0, 0, math.nan, math.nan
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "REQ":
            calls += 1
        elif parts[0] == "BATCH":
            batches += 1
            elapsed.append(float(parts[8]))
    return calls, batches, sum(elapsed) if elapsed else math.nan, statistics.fmean(elapsed) if elapsed else math.nan


result = json.loads(result_path.read_text(encoding="utf-8"))
calls, batches, batch_sum, batch_mean = oracle_timing(timing)
wall_s = float(result["wall_s"])
force_evals = 1 + 4 * steps
bead_force_cases = force_evals * p
row = {
    "backend": backend,
    "P": p,
    "steps": steps,
    "wall_s": wall_s,
    "force_evals": force_evals,
    "bead_force_cases": bead_force_cases,
    "s_per_force_eval": wall_s / force_evals,
    "s_per_bead_force_case": wall_s / bead_force_cases,
    "omp_threads": omp,
    "oracle_workers": workers if backend == "oracle" else 0,
    "oracle_calls": calls,
    "oracle_batches": batches,
    "oracle_batch_sum_s": batch_sum,
    "oracle_batch_mean_s": batch_mean,
    "electron_T_K": last_data(out / "electron_temperature.out", 2),
    "electron_extended_T_K": last_data(out / "electron_temperature.out", 3),
    "carbon_T_K": last_data(out / "carbon_temperature.out", 2),
    "energy_eV": last_data(out / "pimd_energy.out", 3),
    "status": "PASS" if result["status"] == 0 else f"FAIL:{result['status']}",
}

exists = summary.exists()
with summary.open("a", encoding="utf-8") as handle:
    if not exists:
        handle.write("# Fast-vs-Oracle Timing Summary\n\n")
        handle.write(
            "| backend | P | steps | wall_s | force_evals | bead_force_cases | s_per_force_eval | "
            "s_per_bead_force_case | OMP | oracle_workers | oracle_calls | oracle_batches | "
            "oracle_batch_sum_s | oracle_batch_mean_s | electron_T_K | electron_extended_T_K | carbon_T_K | "
            "energy_eV | status |\n"
        )
        handle.write(
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
            "---: | ---: | ---: | ---: | ---: | ---: | --- |\n"
        )
    handle.write(
        f"| {row['backend']} | {row['P']} | {row['steps']} | {row['wall_s']:.3f} | "
        f"{row['force_evals']} | {row['bead_force_cases']} | {row['s_per_force_eval']:.3f} | "
        f"{row['s_per_bead_force_case']:.3f} | {row['omp_threads']} | {row['oracle_workers']} | "
        f"{row['oracle_calls']} | {row['oracle_batches']} | {row['oracle_batch_sum_s']:.3f} | "
        f"{row['oracle_batch_mean_s']:.3f} | {row['electron_T_K']:.3f} | "
        f"{row['electron_extended_T_K']:.3f} | {row['carbon_T_K']:.3f} | "
        f"{row['energy_eV']:.6f} | {row['status']} |\n"
    )
print(json.dumps(row, sort_keys=True))
PY
}

run_case() {
  local backend="$1"
  local p="$2"
  local case_dir="$RUN_ROOT/${backend}_p${p}_steps${STEPS}"
  local cfg="$case_dir/config.yaml"
  local out="$case_dir/out"
  mkdir -p "$case_dir" "$out"
  write_config "$backend" "$p" "$out" "$cfg"

  python3 - "$ROOT_DIR" "$cfg" "$out" "$OMP_THREADS" "$case_dir/stdout.log" "$case_dir/stderr.log" "$case_dir/result.json" <<'PY'
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
omp = sys.argv[4]
stdout_log = pathlib.Path(sys.argv[5])
stderr_log = pathlib.Path(sys.argv[6])
result_json = pathlib.Path(sys.argv[7])
env = os.environ.copy()
env["OMP_NUM_THREADS"] = omp
cmd = [str(root / "build/ereaxff_mixed_baoab_main"), cfg, out]
start = time.time()
with stdout_log.open("w", encoding="utf-8") as stdout, stderr_log.open("w", encoding="utf-8") as stderr:
    proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env)
wall_s = time.time() - start
result_json.write_text(
    json.dumps(
        {
            "command": cmd,
            "omp_threads": int(omp),
            "status": proc.returncode,
            "wall_s": round(wall_s, 6),
            "config": cfg,
            "out_dir": out,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
raise SystemExit(proc.returncode)
PY
  summarize_case "$backend" "$p" "$case_dir"
}

for p in $FAST_CASES; do
  run_case fast "$p"
done

for p in $ORACLE_CASES; do
  run_case oracle "$p"
done

echo "summary: $SUMMARY"
echo "run_root: $RUN_ROOT"
