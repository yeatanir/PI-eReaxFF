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
cat > "$CALL_DIR/ffield" <<'EOF'
test fast backend ffield
 39       ! Number of general parameters
   50.0000
    9.5469
    1.6725
    1.7224
    6.8702
   54.6742
    1.0588
    4.6000
   12.1176
   13.3056
  -34.5448
    0.0000 !Lower Taper-radius
   10.0000 !Upper Taper-radius
    2.8793
   33.8667
    6.0891
    1.0563
    2.0384
    6.1431
    6.9290
    0.3989
    3.9954
   -2.4837
    8.5385
    6.7491
    0.1414
    1.0000
    1.1348
    1.5591
    0.1000
    1.7602
    0.6991
   50.0000
    1.8512
  548.6451
    0.0001
    1.5000
    0.0000
    0.7903
 2    ! Nr of atoms
 header 1
 header 2
 header 3
 C    1.3417   4.0000  12.0000   2.0464   0.1064   0.3500   1.1681   4.0000
      9.0000   1.5000   4.0000  40.0000   0.8000   5.3422   4.5000   3.1838
      1.1437   0.0000 181.0000   5.4236  19.2788  13.5366   3.1838  22.0000
     -4.0998   4.8750   0.0027   4.0000   2.9663   0.0000   0.0000   0.0000
 El  -0.1000   0.0000   0.5000   1.3000   0.0050   0.2500  -0.1000   0.0000
      2.0000   4.0000   4.0000   0.0000   0.0000   1.0002   1.0009  10.0000
     -0.1000   0.0000  -2.3700   6.4918   8.5961   0.2368   2.0000   3.5000
     -5.0000   3.1873   0.1600   6.2998   2.5791   0.0000   0.0000   0.0000
 1    ! Nr of bonds
 header
  1  1  78.4266 115.3834  68.1631   0.5777  -0.2901   1.0000  34.9989   0.5101
         3.8560  -0.1640   8.2326   1.0000  -0.0585   6.7997   1.0000   0.1481
 0    ! Nr of off-diagonal terms
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

blocked = subprocess.run(
    [".venv-fast/bin/python", "projects/jamil_replication/fast_backend/jax_backend.py"],
    input=json.dumps({"command": "force_batch", "force_kernel": "fixed_charge_coulomb", "cases": [case]}),
    text=True,
    capture_output=True,
)
assert blocked.returncode == 2
assert "experimental force kernels require --allow-unvalidated" in blocked.stdout

completed = subprocess.run(
    [
        ".venv-fast/bin/python",
        "projects/jamil_replication/fast_backend/jax_backend.py",
        "--allow-unvalidated",
    ],
    input=json.dumps({"command": "force_batch", "force_kernel": "fixed_charge_coulomb", "cases": [case]}),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["fixed_charge_coulomb"]
assert abs(response["energy_eV"][0] + 3.599911369606375) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][0][0] - 3.599911369606375) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][1][0] + 3.599911369606375) < 1.0e-12

completed = subprocess.run(
    [
        ".venv-fast/bin/python",
        "projects/jamil_replication/fast_backend/jax_backend.py",
        "--allow-unvalidated",
    ],
    input=json.dumps({"command": "force_batch", "force_kernel": "reaxff_shielded_coulomb", "cases": [case]}),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["reaxff_shielded_coulomb"]
assert abs(response["energy_eV"][0] + 1.0529536968375994) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][0][0] - 0.03734168384345882) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][1][0] + 0.03734168384345882) < 1.0e-12

import sys
sys.path.insert(0, "projects/jamil_replication/fast_backend")
from jax_force_kernel import evaluate_reaxff_shielded_coulomb_case

eps = 1.0e-5
case_plus = json.loads(json.dumps(case))
case_minus = json.loads(json.dumps(case))
case_plus["positions_A"][0][0] += eps
case_minus["positions_A"][0][0] -= eps
e_plus = float(evaluate_reaxff_shielded_coulomb_case(case_plus)["energy_eV"])
e_minus = float(evaluate_reaxff_shielded_coulomb_case(case_minus)["energy_eV"])
fd_force = -(e_plus - e_minus) / (2.0 * eps)
analytic_force = float(evaluate_reaxff_shielded_coulomb_case(case)["forces_eV_per_A"][0][0])
assert abs(fd_force - analytic_force) < 1.0e-8

completed = subprocess.run(
    [
        ".venv-fast/bin/python",
        "projects/jamil_replication/fast_backend/jax_backend.py",
        "--allow-unvalidated",
    ],
    input=json.dumps(
        {
            "command": "force_batch",
            "force_kernel": "reaxff_nonbonded_no_inner_shielded",
            "cases": [case],
        }
    ),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["reaxff_nonbonded_no_inner_shielded"]
assert abs(response["energy_eV"][0] + 1.0453112428163907) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][0][0] - 0.02335580645549877) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][1][0] + 0.02335580645549877) < 1.0e-12

from jax_force_kernel import evaluate_reaxff_no_inner_shielded_nonbonded_case

e_plus = float(evaluate_reaxff_no_inner_shielded_nonbonded_case(case_plus)["energy_eV"])
e_minus = float(evaluate_reaxff_no_inner_shielded_nonbonded_case(case_minus)["energy_eV"])
fd_force = -(e_plus - e_minus) / (2.0 * eps)
analytic_force = float(evaluate_reaxff_no_inner_shielded_nonbonded_case(case)["forces_eV_per_A"][0][0])
assert abs(fd_force - analytic_force) < 1.0e-8

bond_case = json.loads(json.dumps(case))
bond_case["atom_types"] = ["C", "C"]
bond_case["positions_A"] = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
bond_case["charge_constraints"]["fixed_charges"] = [0.0, 0.0]
bond_case["charge_in"] = [0.0, 0.0]
bond_case["solver_guess_in"] = [0.0, 0.0]
completed = subprocess.run(
    [
        ".venv-fast/bin/python",
        "projects/jamil_replication/fast_backend/jax_backend.py",
        "--allow-unvalidated",
    ],
    input=json.dumps(
        {
            "command": "force_batch",
            "force_kernel": "reaxff_uncorrected_bond",
            "cases": [bond_case],
        }
    ),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["reaxff_uncorrected_bond"]
assert abs(response["energy_eV"][0] + 6.070757166945563) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][0][0] - 9.444906995307184) < 1.0e-12
assert abs(response["forces_eV_per_A"][0][1][0] + 9.444906995307184) < 1.0e-12

from jax_force_kernel import evaluate_reaxff_uncorrected_bond_case

bond_plus = json.loads(json.dumps(bond_case))
bond_minus = json.loads(json.dumps(bond_case))
bond_plus["positions_A"][0][0] += eps
bond_minus["positions_A"][0][0] -= eps
e_plus = float(evaluate_reaxff_uncorrected_bond_case(bond_plus)["energy_eV"])
e_minus = float(evaluate_reaxff_uncorrected_bond_case(bond_minus)["energy_eV"])
fd_force = -(e_plus - e_minus) / (2.0 * eps)
analytic_force = float(evaluate_reaxff_uncorrected_bond_case(bond_case)["forces_eV_per_A"][0][0])
assert abs(fd_force - analytic_force) < 1.0e-8

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
