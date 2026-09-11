#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/projects/jamil_replication"
BASE_CFG="$SCRIPT_DIR/config_jamil_112c_113el_p4_pool_verify.yaml"
RUN_ROOT="${JAMIL_FAST_STABILITY_ROOT:-$ROOT_DIR/build/jamil_replication/fast_stability_$(date +%Y%m%d_%H%M%S)}"
SUMMARY="$SCRIPT_DIR/fast_stability_summary.out"
P="${JAMIL_FAST_STABILITY_P:-16}"
MAX_CARBON_T_K="${JAMIL_FAST_STABILITY_MAX_CARBON_T_K:-1000.0}"
MAX_BEAD_SPREAD_A="${JAMIL_FAST_STABILITY_MAX_BEAD_SPREAD_A:-2.0}"
DT_CASES="${JAMIL_FAST_STABILITY_DT_CASES:-0.05 0.1}"
C_THERMOSTAT="${JAMIL_FAST_STABILITY_C_THERMOSTAT:-}"
C_TAU_FS="${JAMIL_FAST_STABILITY_C_TAU_FS:-}"
C_GAMMA_FS_INV="${JAMIL_FAST_STABILITY_C_GAMMA_FS_INV:-}"
if [[ "${JAMIL_FAST_STABILITY_RUN_LONG:-0}" == "1" ]]; then
  STEP_CASES="${JAMIL_FAST_STABILITY_STEP_CASES:-10 100}"
else
  STEP_CASES="${JAMIL_FAST_STABILITY_STEP_CASES:-10}"
fi

if [[ ! -f "$BASE_CFG" ]]; then
  echo "run_fast_stability_check: missing base config: $BASE_CFG" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/build/ereaxff_mixed_baoab_main" ]]; then
  echo "run_fast_stability_check: missing executable; run make ereaxff_mixed_baoab_main" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/.venv-fast/bin/python" ]]; then
  echo "run_fast_stability_check: missing .venv-fast/bin/python" >&2
  exit 77
fi

mkdir -p "$RUN_ROOT"
rm -f "$SUMMARY"

write_config() {
  local steps="$1"
  local dt_fs="$2"
  local out="$3"
  local cfg="$4"

  python3 - "$BASE_CFG" "$cfg" "$P" "$steps" "$dt_fs" "$out" "$C_THERMOSTAT" "$C_TAU_FS" "$C_GAMMA_FS_INV" <<'PY'
from __future__ import annotations

import pathlib
import re
import sys

src = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
dst = pathlib.Path(sys.argv[2])
p = sys.argv[3]
steps = sys.argv[4]
dt_fs = sys.argv[5]
out = sys.argv[6]
c_thermostat = sys.argv[7]
c_tau_fs = sys.argv[8]
c_gamma_fs_inv = sys.argv[9]


def replace(pattern: str, value: str, text: str) -> str:
    return re.sub(pattern, value, text, flags=re.M)


text = src
text = replace(r"^  type:.*$", "  type: fast_jax_oracle", text)
text = replace(r"^  dt_fs:.*$", f"  dt_fs: {dt_fs}", text)
text = replace(r"^  steps:.*$", f"  steps: {steps}", text)
text = replace(r"^  output_stride:.*$", "  output_stride: 1", text)
if "rescale_initial_pimd_velocities:" in text:
    text = replace(r"^  rescale_initial_pimd_velocities:.*$", "  rescale_initial_pimd_velocities: true", text)
else:
    text = replace(r"^  restart_stride:.*$", "  restart_stride: 1\n  rescale_initial_pimd_velocities: true", text)
text = replace(r"^  output_dir:.*$", f'  output_dir: "{out}"', text)
text = replace(r"^    beads:.*$", f"    beads: {p}", text)
if c_thermostat:
    text = re.sub(
        r"(  - name: graphene_carbons\n"
        r"    selector: \"element == C\"\n"
        r"    dynamics: classical\n"
        r"    mobile: true\n)"
        r"    thermostat: .*\n"
        r"    temperature_K: .*\n"
        r"(?:    tau_fs: .*\n|    gamma_fs_inv: .*\n)?",
        lambda match: (
            match.group(1)
            + f"    thermostat: {c_thermostat}\n"
            + "    temperature_K: 300.0\n"
            + (
                f"    tau_fs: {c_tau_fs or '100.0'}\n"
                if c_thermostat == "berendsen"
                else f"    gamma_fs_inv: {c_gamma_fs_inv or '0.001'}\n"
            )
        ),
        text,
        count=1,
    )
insertion = (
    "  repeatability_check: false\n"
    '  fast_backend_command: ".venv-fast/bin/python projects/jamil_replication/fast_backend/jax_live_worker.py"\n'
    '  fast_backend_model: "reaxff_corrected_bond_nonbonded_valence"\n'
    "  fast_backend_allow_unvalidated: true"
)
text = replace(r"^  repeatability_check:.*$", insertion, text)
dst.write_text(text, encoding="utf-8")
PY
}

summarize_case() {
  local steps="$1"
  local dt_fs="$2"
  local case_dir="$3"
  local out="$case_dir/out"
  local result_json="$case_dir/result.json"

  python3 - "$P" "$steps" "$dt_fs" "$MAX_CARBON_T_K" "$MAX_BEAD_SPREAD_A" "$result_json" "$out" "$SUMMARY" <<'PY'
from __future__ import annotations

import csv
import json
import math
import pathlib
import sys

p = int(sys.argv[1])
steps = int(sys.argv[2])
dt_fs = float(sys.argv[3])
max_carbon_t_k = float(sys.argv[4])
max_bead_spread_a = float(sys.argv[5])
result = json.loads(pathlib.Path(sys.argv[6]).read_text(encoding="utf-8"))
out = pathlib.Path(sys.argv[7])
summary = pathlib.Path(sys.argv[8])


def rows(path: pathlib.Path) -> list[list[float]]:
    data = []
    if not path.exists():
        return data
    for raw in path.read_text(encoding="utf-8").splitlines():
        parts = raw.split()
        if parts and not parts[0].startswith("#"):
            data.append([float(x) for x in parts])
    return data


def last_col(path: pathlib.Path, col: int) -> float:
    data = rows(path)
    return data[-1][col - 1] if data else math.nan


def minmax_col(path: pathlib.Path, col: int) -> tuple[float, float]:
    data = rows(path)
    if not data:
        return math.nan, math.nan
    values = [row[col - 1] for row in data]
    return min(values), max(values)


def all_finite(paths: list[pathlib.Path]) -> bool:
    for path in paths:
        if not path.exists():
            return False
        for row in rows(path):
            if any(not math.isfinite(value) for value in row):
                return False
    return True


timing_path = out / "fast_jax_work" / "worker_timing.csv"
timing_rows = []
if timing_path.exists():
    with timing_path.open(newline="", encoding="utf-8") as fh:
        timing_rows = list(csv.DictReader(fh))
force_total_s = sum(float(row["force_s"]) for row in timing_rows)
request_total_s = sum(float(row["total_s"]) for row in timing_rows)
startup_s = max(0.0, float(result["wall_s"]) - request_total_s)
force_evals = 1 + 4 * steps
bead_force_cases = force_evals * p
spread_min, spread_max = minmax_col(out / "bead_spread.out", 5)
conserved_min, conserved_max = minmax_col(out / "pimd_energy.out", 5)
max_force = max(abs(row[6]) for row in rows(out / "pimd_energy.out")) if rows(out / "pimd_energy.out") else math.nan
finite = all_finite(
    [
        out / "pimd_energy.out",
        out / "bead_spread.out",
        out / "electron_temperature.out",
        out / "carbon_temperature.out",
        out / "thermostat_heat.out",
    ]
)
xyz_ok = (out / "trajectory_centroid.xyz").stat().st_size > 0 and (out / "trajectory_beads.xyz").stat().st_size > 0
worker_done = (out / "fast_jax_work" / "worker.done").exists()
carbon_last = last_col(out / "carbon_temperature.out", 2)
spread_last = last_col(out / "bead_spread.out", 5)
thermal_ok = carbon_last <= max_carbon_t_k and spread_last <= max_bead_spread_a
status = "PASS" if result["status"] == 0 and finite and xyz_ok and worker_done and thermal_ok else "FAIL"
row = {
    "P": p,
    "steps": steps,
    "dt_fs": dt_fs,
    "wall_s": float(result["wall_s"]),
    "force_evals": force_evals,
    "bead_force_cases": bead_force_cases,
    "s_per_step": float(result["wall_s"]) / max(steps, 1),
    "s_per_force_eval": float(result["wall_s"]) / force_evals,
    "worker_force_sum_s": force_total_s,
    "worker_request_sum_s": request_total_s,
    "estimated_startup_s": startup_s,
    "electron_T_last_K": last_col(out / "electron_temperature.out", 2),
    "electron_extended_T_last_K": last_col(out / "electron_temperature.out", 3),
    "carbon_T_last_K": carbon_last,
    "bead_spread_last_A": spread_last,
    "bead_spread_min_A": spread_min,
    "bead_spread_max_A": spread_max,
    "conserved_span_eV": conserved_max - conserved_min,
    "max_force_eV_per_A": max_force,
    "finite": finite,
    "thermal_ok": thermal_ok,
    "xyz_ok": xyz_ok,
    "worker_done": worker_done,
    "status": status,
}

exists = summary.exists()
with summary.open("a", encoding="utf-8") as handle:
    if not exists:
        handle.write("# Fast JAX Stability Summary\n\n")
        handle.write(
            "| P | steps | dt_fs | wall_s | s/step | s/force_eval | worker_force_sum_s | "
            "worker_request_sum_s | startup_s | electron_T_last_K | extended_T_last_K | carbon_T_last_K | "
            "bead_spread_last_A | bead_spread_min_A | bead_spread_max_A | conserved_span_eV | "
            "max_force_eV_per_A | finite | thermal_ok | xyz | worker_done | status |\n"
        )
        handle.write(
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
            "---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- | --- |\n"
        )
    handle.write(
        f"| {row['P']} | {row['steps']} | {row['dt_fs']:.3f} | {row['wall_s']:.3f} | "
        f"{row['s_per_step']:.3f} | {row['s_per_force_eval']:.3f} | {row['worker_force_sum_s']:.3f} | "
        f"{row['worker_request_sum_s']:.3f} | {row['estimated_startup_s']:.3f} | "
        f"{row['electron_T_last_K']:.3f} | {row['electron_extended_T_last_K']:.3f} | "
        f"{row['carbon_T_last_K']:.3f} | {row['bead_spread_last_A']:.6f} | "
        f"{row['bead_spread_min_A']:.6f} | {row['bead_spread_max_A']:.6f} | "
        f"{row['conserved_span_eV']:.6f} | {row['max_force_eV_per_A']:.6e} | "
        f"{row['finite']} | {row['thermal_ok']} | {row['xyz_ok']} | {row['worker_done']} | {row['status']} |\n"
    )
print(json.dumps(row, sort_keys=True))
if status != "PASS":
    raise SystemExit(2)
PY
}

run_case() {
  local steps="$1"
  local dt_fs="$2"
  local dt_token="${dt_fs//./p}"
  local case_dir="$RUN_ROOT/fast_p${P}_steps${steps}_dt${dt_token}"
  local cfg="$case_dir/config.yaml"
  local out="$case_dir/out"
  rm -rf "$case_dir"
  mkdir -p "$case_dir" "$out"
  write_config "$steps" "$dt_fs" "$out" "$cfg"

  python3 - "$ROOT_DIR" "$cfg" "$out" "$case_dir/stdout.log" "$case_dir/stderr.log" "$case_dir/result.json" <<'PY'
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
stdout_log = pathlib.Path(sys.argv[4])
stderr_log = pathlib.Path(sys.argv[5])
result_json = pathlib.Path(sys.argv[6])
env = os.environ.copy()
env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "4")
cmd = [str(root / "build/ereaxff_mixed_baoab_main"), cfg, out]
start = time.time()
with stdout_log.open("w", encoding="utf-8") as stdout, stderr_log.open("w", encoding="utf-8") as stderr:
    proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env)
wall_s = time.time() - start
result_json.write_text(
    json.dumps(
        {
            "command": cmd,
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
  summarize_case "$steps" "$dt_fs" "$case_dir"
}

for dt_fs in $DT_CASES; do
  for steps in $STEP_CASES; do
    run_case "$steps" "$dt_fs"
  done
done

echo "summary: $SUMMARY"
echo "run_root: $RUN_ROOT"
