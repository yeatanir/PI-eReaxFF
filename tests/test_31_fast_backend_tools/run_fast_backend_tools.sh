#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

WORK_DIR="build/test_31_fast_backend_tools/adapter_work"
CALL_DIR="$WORK_DIR/persistent_worker_001/adapter_call_0001"
OUT_JSONL="build/test_31_fast_backend_tools/labeled_cases.jsonl"

mkdir -p "$CALL_DIR"
cat > "$WORK_DIR/persistent_oracle_timing.log" <<'EOF'
# kind request step bead stage worker_state elapsed_s status natom
CALL 0 1 1 1 1 bounded_pool bead_local initial_force pool 0.0 0.1 0.1 0 build/test_31_fast_backend_tools/adapter_work/adapter_call_0001
EOF
cat > "$CALL_DIR/geo" <<'EOF'
BIOGRF 200
DESCRP tiny
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
CRYSTX      9.00000      9.00000      9.00000     90.00000     90.00000     90.00000
HETATM     1 C           C       0.00000   0.00000   0.00000        0 0  0.00000
HETATM     2 El          El      1.00000   0.00000   0.00000        0 0  0.00000
MOLCHARGE    1    1   0.50
MOLCHARGE    2    2  -0.50
END
EOF
cat > "$CALL_DIR/charge_memory.in" <<'EOF'
1 0.0 0.0
2 0.0 0.0
EOF
cat > "$CALL_DIR/phase4_single_point_oracle.out" <<'EOF'
# energy_kcal max_force charge_sum status natom
 -1.000000000000000E+01  2.000000000000000E+00  0.000000000000000E+00 0 2
EOF
cat > "$CALL_DIR/control" <<'EOF'
      8 icharg
EOF
cat > "$CALL_DIR/single_point_forces.out" <<'EOF'
# atom Fx Fy Fz charge solver_guess
1  1.0  0.0  0.0  0.5  0.0
2 -1.0  0.0  0.0 -0.5  0.0
EOF

python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py \
  --adapter-work "$WORK_DIR" \
  --output "$OUT_JSONL"
python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py --check "$OUT_JSONL"

if [[ -x ".venv-fast/bin/python" ]]; then
  .venv-fast/bin/python projects/jamil_replication/fast_backend/jax_backend.py \
    --parity-check \
    --dataset "$OUT_JSONL"
  .venv-fast/bin/python - <<'PY'
import json
import subprocess
from pathlib import Path

case = json.loads(Path("build/test_31_fast_backend_tools/labeled_cases.jsonl").read_text())
request = json.dumps({"command": "charge_batch", "cases": [case, case]})
completed = subprocess.run(
    [".venv-fast/bin/python", "projects/jamil_replication/fast_backend/jax_backend.py"],
    input=request,
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["charge_models"] == ["fixed_molcharge"]
assert response["charge_out"] == [[0.5, -0.5], [0.5, -0.5]]
assert response["solver_guess_out"] == [[0.0, 0.0], [0.0, 0.0]]

failed = subprocess.run(
    [".venv-fast/bin/python", "projects/jamil_replication/fast_backend/jax_backend.py"],
    input=json.dumps({"command": "force_batch", "cases": [case]}),
    text=True,
    capture_output=True,
)
assert failed.returncode == 2
assert "physics force kernel is not implemented" in failed.stdout

completed = subprocess.run(
    [".venv-fast/bin/python", "projects/jamil_replication/fast_backend/jax_backend.py"],
    input=json.dumps({"command": "force_batch", "allow_oracle_replay": True, "cases": [case, case]}),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["oracle_replay"]
assert response["energy_eV"] == [-0.43364, -0.43364]
assert response["forces_eV_per_A"][0] == [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
PY
else
  python3 projects/jamil_replication/fast_backend/jax_backend.py \
    --self-test \
    --dataset "$OUT_JSONL"
fi

echo "test_31_fast_backend_tools: PASS"
