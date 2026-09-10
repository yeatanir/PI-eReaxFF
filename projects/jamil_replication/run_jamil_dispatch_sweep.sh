#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/projects/jamil_replication"
BASE_CFG="$SCRIPT_DIR/config_jamil_112c_113el_p16_300K_prod.yaml"
SUMMARY_SCRIPT="$SCRIPT_DIR/summarize_jamil_dispatch_sweep.py"
SWEEP_ROOT="${JAMIL_DISPATCH_SWEEP_ROOT:-$ROOT_DIR/build/jamil_replication/dispatch_sweep_$(date +%Y%m%d_%H%M%S)}"
FINAL_SUMMARY="$SCRIPT_DIR/jamil_dispatch_sweep_summary.out"
TIMEOUT_S="${JAMIL_DISPATCH_CASE_TIMEOUT_S:-0}"

if [[ ! -f "$BASE_CFG" ]]; then
  echo "run_jamil_dispatch_sweep: missing base config: $BASE_CFG" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/build/ereaxff_mixed_baoab_main" ]]; then
  echo "run_jamil_dispatch_sweep: missing executable: $ROOT_DIR/build/ereaxff_mixed_baoab_main" >&2
  exit 1
fi
if [[ ! -x "$SUMMARY_SCRIPT" ]]; then
  chmod +x "$SUMMARY_SCRIPT"
fi

mkdir -p "$SWEEP_ROOT"

declare -a cases=(
  "pool_omp1_w1|persistent_worker_pool|1|bounded_pool|1"
  "pool_omp1_w2|persistent_worker_pool|1|bounded_pool|2"
  "pool_omp1_w4|persistent_worker_pool|1|bounded_pool|4"
  "pool_omp1_w6|persistent_worker_pool|1|bounded_pool|6"
  "pool_omp1_w8|persistent_worker_pool|1|bounded_pool|8"
)

if [[ "${JAMIL_DISPATCH_SKIP_W1:-0}" == "1" ]]; then
  cases=(
    "pool_omp1_w2|persistent_worker_pool|1|bounded_pool|2"
    "pool_omp1_w4|persistent_worker_pool|1|bounded_pool|4"
    "pool_omp1_w6|persistent_worker_pool|1|bounded_pool|6"
    "pool_omp1_w8|persistent_worker_pool|1|bounded_pool|8"
  )
fi

if [[ "${JAMIL_INCLUDE_OMP2:-0}" == "1" ]]; then
  cases+=(
    "pool_omp2_w2|persistent_worker_pool|2|bounded_pool|2"
    "pool_omp2_w4|persistent_worker_pool|2|bounded_pool|4"
  )
fi

if [[ "${JAMIL_INCLUDE_OPENMP_ALL:-0}" == "1" ]]; then
  cases+=(
    "external_omp1_openmp_all|external_process|1|openmp_all|1"
  )
fi

run_case() {
  local case_name="$1"
  local backend_mode="$2"
  local omp_threads="$3"
  local dispatch_mode="$4"
  local max_workers="$5"
  local case_dir="$SWEEP_ROOT/$case_name"
  local out_dir="$case_dir/out"
  local cfg="$case_dir/config.yaml"
  local stdout_log="$case_dir/stdout.log"
  local stderr_log="$case_dir/stderr.log"
  local result_json="$case_dir/result.json"
  local start_s end_s wall_s status note

  mkdir -p "$case_dir" "$out_dir"

  python3 - "$BASE_CFG" "$cfg" "$backend_mode" "$dispatch_mode" "$max_workers" <<'PY'
from __future__ import annotations

import pathlib
import re
import sys

src = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
dst = pathlib.Path(sys.argv[2])
backend_mode = sys.argv[3]
dispatch_mode = sys.argv[4]
max_workers = sys.argv[5]

def replace(pattern: str, value: str, text: str) -> str:
    return re.sub(pattern, value, text, flags=re.M)

text = src
text = replace(r"^  oracle_backend_mode:.*$", f"  oracle_backend_mode: {backend_mode}", text)
text = replace(r"^  oracle_dispatch_mode:.*$", f"  oracle_dispatch_mode: {dispatch_mode}", text)
text = replace(r"^  oracle_max_workers:.*$", f"  oracle_max_workers: {max_workers}", text)
text = replace(r"^  oracle_timing_log:.*$", "  oracle_timing_log: true", text)
text = replace(r"^  oracle_stagger_launch_ms:.*$", "  oracle_stagger_launch_ms: 0", text)
text = replace(r"^  oracle_fail_fast:.*$", "  oracle_fail_fast: true", text)
text = replace(r"^  oracle_timeout_s:.*$", "  oracle_timeout_s: 90", text)
dst.write_text(text, encoding="utf-8")
PY

  start_s="$(python3 - <<'PY'
import time
print(f"{time.time():.6f}")
PY
)"
  set +e
  if [[ "$TIMEOUT_S" -gt 0 ]]; then
    python3 - "$ROOT_DIR" "$cfg" "$out_dir" "$omp_threads" "$TIMEOUT_S" "$stdout_log" "$stderr_log" "$result_json" <<'PY'
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys
import time

root = pathlib.Path(sys.argv[1])
cfg = sys.argv[2]
out_dir = sys.argv[3]
omp_threads = sys.argv[4]
timeout_s = float(sys.argv[5])
stdout_log = pathlib.Path(sys.argv[6])
stderr_log = pathlib.Path(sys.argv[7])
result_json = pathlib.Path(sys.argv[8])

env = os.environ.copy()
env["OMP_NUM_THREADS"] = omp_threads
cmd = [str(root / "build/ereaxff_mixed_baoab_main"), cfg, out_dir]
start = time.time()
status = 0
note = "ok"
try:
    with stdout_log.open("w", encoding="utf-8") as out, stderr_log.open("w", encoding="utf-8") as err:
        proc = subprocess.run(cmd, stdout=out, stderr=err, env=env, timeout=timeout_s)
    status = proc.returncode
    if status != 0:
        note = f"exit_{status}"
except subprocess.TimeoutExpired:
    status = 124
    note = "timeout"
end = time.time()
payload = {
    "command": cmd,
    "omp_threads": int(omp_threads),
    "status": status,
    "note": note,
    "wall_s": round(end - start, 6),
    "stdout_log": str(stdout_log),
    "stderr_log": str(stderr_log),
    "config": cfg,
    "out_dir": out_dir,
}
result_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(note)
sys.exit(0)
PY
    status=$?
  else
    python3 - "$ROOT_DIR" "$cfg" "$out_dir" "$omp_threads" "$stdout_log" "$stderr_log" "$result_json" <<'PY'
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import sys
import time

root = pathlib.Path(sys.argv[1])
cfg = sys.argv[2]
out_dir = sys.argv[3]
omp_threads = sys.argv[4]
stdout_log = pathlib.Path(sys.argv[5])
stderr_log = pathlib.Path(sys.argv[6])
result_json = pathlib.Path(sys.argv[7])

env = os.environ.copy()
env["OMP_NUM_THREADS"] = omp_threads
cmd = [str(root / "build/ereaxff_mixed_baoab_main"), cfg, out_dir]
start = time.time()
status = 0
note = "ok"
try:
    with stdout_log.open("w", encoding="utf-8") as out, stderr_log.open("w", encoding="utf-8") as err:
        proc = subprocess.run(cmd, stdout=out, stderr=err, env=env)
    status = proc.returncode
    if status != 0:
        note = f"exit_{status}"
except Exception as exc:
    status = 125
    note = f"error:{exc.__class__.__name__}"
end = time.time()
payload = {
    "command": cmd,
    "omp_threads": int(omp_threads),
    "status": status,
    "note": note,
    "wall_s": round(end - start, 6),
    "stdout_log": str(stdout_log),
    "stderr_log": str(stderr_log),
    "config": cfg,
    "out_dir": out_dir,
}
result_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(note)
sys.exit(0)
PY
    status=$?
  fi
  set -e
  end_s="$(python3 - <<'PY'
import time
print(f"{time.time():.6f}")
PY
)"
  wall_s="$(python3 - "$start_s" "$end_s" <<'PY'
from __future__ import annotations

import sys
start = float(sys.argv[1])
end = float(sys.argv[2])
print(f"{end - start:.6f}")
PY
)"
  note="$(python3 - "$result_json" <<'PY'
from __future__ import annotations

import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    print("missing_result")
else:
    print(json.loads(path.read_text(encoding="utf-8"))["note"])
PY
)"
  echo "$case_name $backend_mode $omp_threads $dispatch_mode $max_workers wall_s=$wall_s status=$status note=$note"
}

for spec in "${cases[@]}"; do
  IFS='|' read -r case_name backend_mode omp_threads dispatch_mode max_workers <<<"$spec"
  run_case "$case_name" "$backend_mode" "$omp_threads" "$dispatch_mode" "$max_workers" || true
done

"$SUMMARY_SCRIPT" "$SWEEP_ROOT" "$FINAL_SUMMARY"

echo "sweep_root=$SWEEP_ROOT"
echo "summary=$FINAL_SUMMARY"
