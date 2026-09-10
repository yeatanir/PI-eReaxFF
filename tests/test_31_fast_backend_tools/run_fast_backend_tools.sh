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
cat > "$CALL_DIR/single_point_energy_terms.out" <<'EOF'
# term energy_kcal
bond -4.0
overcoordination 0.5
lone_pair 0.25
molecule 0.0
valence_penalty 1.0
angle_conjugation -0.25
hydrogen_bond 0.0
torsion 0.75
torsion_conjugation -0.125
vdw 2.0
coulomb -9.625
charge 0.0
external_field 0.0
electron -0.5
total -10.0
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
import math
import subprocess
from pathlib import Path

case = json.loads(Path("build/test_31_fast_backend_tools/labeled_cases.jsonl").read_text())
terms = case["labels"]["energy_terms_kcal"]
assert terms["bond"] == -4.0
assert terms["electron"] == -0.5
assert abs(sum(case["labels"]["energy_terms_eV"].values()) - case["labels"]["energy_eV"]) < 1.0e-12

import sys
sys.path.insert(0, "projects/jamil_replication/fast_backend")
from build_oracle_dataset import parse_energy_terms

fort73 = Path("build/test_31_fast_backend_tools/fort73_sample.out")
fort73.write_text("7 -4.0 0.5 0.25 0.0 1.0 -0.25 0.0 0.75 -0.125 2.0 -9.625 0.0 0.0 -0.5\n")
_, fort73_terms_eV = parse_energy_terms(fort73)
assert abs(sum(fort73_terms_eV.values()) - case["labels"]["energy_eV"]) < 1.0e-12

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

from jax_force_kernel import evaluate_reaxff_shielded_coulomb_case
from jax_force_kernel import evaluate_reaxff_corrected_bond_nonbonded_valence_case
from jax_force_kernel import parse_ffield_coulomb_parameters
from jax_force_kernel import _bond_energy_inputs
from jax_force_kernel import _orthorhombic_cell_lengths
from jax_force_kernel import _sparse_topology_from_bond_orders

parameters = parse_ffield_coulomb_parameters(Path(case["source_call_dir"]) / "ffield")
assert len(parameters["angle_params"]) == 0
assert len(parameters["torsion_params"]) == 0
assert len(parameters["hydrogen_bond_params"]) == 0
real_dataset = Path("build/jamil_replication/fast_backend/oracle_labeled_p4_latest.jsonl")
if real_dataset.exists():
    real_case = json.loads(real_dataset.read_text().splitlines()[0])
    real_parameters = parse_ffield_coulomb_parameters(Path(real_case["source_call_dir"]) / "ffield")
    assert len(real_parameters["angle_params"]) == 101
    assert len(real_parameters["torsion_params"]) == 56
    assert len(real_parameters["hydrogen_bond_params"]) == 0

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
assert set(response["energy_terms_eV"]) == {"coulomb", "vdw"}
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

corrected_case = json.loads(json.dumps(case))
r = 1.4
diag = r / math.sqrt(3.0)
corrected_case["atom_types"] = ["C"] * 9
corrected_case["positions_A"] = [
    [0.0, 0.0, 0.0],
    [r, 0.0, 0.0],
    [-r, 0.0, 0.0],
    [0.0, r, 0.0],
    [0.0, -r, 0.0],
    [0.0, 0.0, r],
    [0.0, 0.0, -r],
    [diag, diag, diag],
    [-diag, -diag, -diag],
]
corrected_case["charge_constraints"]["fixed_charges"] = [0.0] * 9
corrected_case["charge_in"] = [0.0] * 9
corrected_case["solver_guess_in"] = [0.0] * 9
completed = subprocess.run(
    [
        ".venv-fast/bin/python",
        "projects/jamil_replication/fast_backend/jax_backend.py",
        "--allow-unvalidated",
    ],
    input=json.dumps(
        {
            "command": "force_batch",
            "force_kernel": "reaxff_corrected_bond",
            "cases": [corrected_case],
        }
    ),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["reaxff_corrected_bond"]
assert set(response["energy_terms_eV"]) == {"bond"}
assert abs(response["energy_eV"][0] + 77.82789596632907) < 1.0e-10
assert abs(response["forces_eV_per_A"][0][1][0] - 6.588497625294385) < 1.0e-10
assert abs(response["forces_eV_per_A"][0][7][0] + 2.9441751399842717) < 1.0e-10

from jax_force_kernel import evaluate_reaxff_corrected_bond_case

corrected_plus = json.loads(json.dumps(corrected_case))
corrected_minus = json.loads(json.dumps(corrected_case))
corrected_plus["positions_A"][1][0] += eps
corrected_minus["positions_A"][1][0] -= eps
e_plus = float(evaluate_reaxff_corrected_bond_case(corrected_plus)["energy_eV"])
e_minus = float(evaluate_reaxff_corrected_bond_case(corrected_minus)["energy_eV"])
fd_force = -(e_plus - e_minus) / (2.0 * eps)
analytic_force = float(evaluate_reaxff_corrected_bond_case(corrected_case)["forces_eV_per_A"][1][0])
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
            "force_kernel": "reaxff_corrected_bond_nonbonded",
            "cases": [corrected_case],
        }
    ),
    text=True,
    capture_output=True,
    check=True,
)
response = json.loads(completed.stdout)
assert response["status"] == 0
assert response["kernel"] == ["reaxff_corrected_bond_nonbonded"]
assert set(response["energy_terms_eV"]) == {"bond", "coulomb", "vdw"}
term_sum = sum(response["energy_terms_eV"][key][0] for key in response["energy_terms_eV"])
assert abs(term_sum - response["energy_eV"][0]) < 1.0e-12
assert abs(response["energy_eV"][0] + 62.15205963429547) < 1.0e-10
assert abs(response["forces_eV_per_A"][0][1][0] - 11.077932259015485) < 1.0e-10
assert abs(response["forces_eV_per_A"][0][7][0] - 0.12878014525015002) < 1.0e-10

from jax_force_kernel import evaluate_reaxff_corrected_bond_nonbonded_case

combined_plus = json.loads(json.dumps(corrected_case))
combined_minus = json.loads(json.dumps(corrected_case))
combined_plus["positions_A"][1][0] += eps
combined_minus["positions_A"][1][0] -= eps
e_plus = float(evaluate_reaxff_corrected_bond_nonbonded_case(combined_plus)["energy_eV"])
e_minus = float(evaluate_reaxff_corrected_bond_nonbonded_case(combined_minus)["energy_eV"])
fd_force = -(e_plus - e_minus) / (2.0 * eps)
analytic_force = float(evaluate_reaxff_corrected_bond_nonbonded_case(corrected_case)["forces_eV_per_A"][1][0])
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

if real_dataset.exists():
    real_case = json.loads(real_dataset.read_text().splitlines()[0])
    angle_case = json.loads(json.dumps(real_case))
    angle = 120.0 * math.pi / 180.0
    angle_case["atom_types"] = ["C", "C", "C"]
    angle_case["positions_A"] = [[0.0, 0.0, 0.0], [r, 0.0, 0.0], [r + r * math.cos(angle), r * math.sin(angle), 0.0]]
    angle_case["charge_constraints"]["fixed_charges"] = [0.0] * 3
    angle_case["charge_in"] = [0.0] * 3
    angle_case["solver_guess_in"] = [0.0] * 3
    angle_case["crystx"] = None
    angle_topology = _sparse_topology_from_bond_orders(
        angle_case,
        real_parameters,
        _bond_energy_inputs(angle_case, real_parameters),
        _orthorhombic_cell_lengths(angle_case),
    )
    assert angle_topology["angle_indices"].shape == (3, 3)
    assert angle_topology["torsion_indices"].shape == (0, 4)
    completed = subprocess.run(
        [
            ".venv-fast/bin/python",
            "projects/jamil_replication/fast_backend/jax_backend.py",
            "--allow-unvalidated",
        ],
        input=json.dumps(
            {
                "command": "force_batch",
                "force_kernel": "reaxff_corrected_bond_nonbonded_valence",
                "cases": [angle_case],
            }
        ),
        text=True,
        capture_output=True,
        check=True,
    )
    response = json.loads(completed.stdout)
    assert response["kernel"] == ["reaxff_corrected_bond_nonbonded_valence"]
    assert set(response["energy_terms_eV"]) == {
        "angle_conjugation",
        "bond",
        "coulomb",
        "electron_core",
        "electron_electron",
        "penalty",
        "torsion",
        "torsion_conjugation",
        "valence",
        "vdw",
    }
    assert abs(response["energy_terms_eV"]["valence"][0] - 4.930837527913495) < 1.0e-10
    assert abs(response["energy_terms_eV"]["penalty"][0] - 0.04994792351223957) < 1.0e-10
    assert abs(response["energy_terms_eV"]["electron_core"][0]) < 1.0e-12
    assert abs(response["energy_terms_eV"]["electron_electron"][0]) < 1.0e-12
    assert abs(response["energy_eV"][0] + 10.791805565593878) < 1.0e-10
    assert abs(response["forces_eV_per_A"][0][2][1] + 12.263151967963509) < 1.0e-10

    completed = subprocess.run(
        [
            ".venv-fast/bin/python",
            "projects/jamil_replication/fast_backend/jax_backend.py",
            "--allow-unvalidated",
        ],
        input=json.dumps(
            {
                "command": "force_batch",
                "force_kernel": "reaxff_corrected_bond_nonbonded_valence",
                "cases": [angle_case, angle_case],
            }
        ),
        text=True,
        capture_output=True,
        check=True,
    )
    batch_response = json.loads(completed.stdout)
    assert batch_response["kernel"] == ["reaxff_corrected_bond_nonbonded_valence"]
    assert batch_response["energy_eV"] == response["energy_eV"] * 2
    for key in response["energy_terms_eV"]:
        assert batch_response["energy_terms_eV"][key] == response["energy_terms_eV"][key] * 2

    angle_eps = 1.0e-6
    angle_plus = json.loads(json.dumps(angle_case))
    angle_minus = json.loads(json.dumps(angle_case))
    angle_plus["positions_A"][2][1] += angle_eps
    angle_minus["positions_A"][2][1] -= angle_eps
    e_plus = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(angle_plus)["energy_eV"])
    e_minus = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(angle_minus)["energy_eV"])
    fd_force = -(e_plus - e_minus) / (2.0 * angle_eps)
    analytic_force = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(angle_case)["forces_eV_per_A"][2][1])
    assert abs(fd_force - analytic_force) < 1.0e-8

    torsion_case = json.loads(json.dumps(real_case))
    torsion_case["atom_types"] = ["C", "C", "C", "C"]
    torsion_case["positions_A"] = [
        [0.0, 0.0, 0.0],
        [1.4, 0.0, 0.0],
        [2.1, 1.2124355653, 0.0],
        [3.4, 1.4, 0.7],
    ]
    torsion_case["charge_constraints"]["fixed_charges"] = [0.0] * 4
    torsion_case["charge_in"] = [0.0] * 4
    torsion_case["solver_guess_in"] = [0.0] * 4
    torsion_case["crystx"] = None
    torsion_topology = _sparse_topology_from_bond_orders(
        torsion_case,
        real_parameters,
        _bond_energy_inputs(torsion_case, real_parameters),
        _orthorhombic_cell_lengths(torsion_case),
    )
    assert torsion_topology["angle_indices"].shape == (8, 3)
    assert torsion_topology["torsion_indices"].shape == (3, 4)
    torsion_output = evaluate_reaxff_corrected_bond_nonbonded_valence_case(torsion_case)
    torsion_terms = {key: float(value) for key, value in torsion_output["energy_terms_eV"].items()}
    assert abs(torsion_terms["torsion"] - 0.1283796427711645) < 1.0e-12
    assert abs(torsion_terms["torsion_conjugation"] + 0.09034874251558697) < 1.0e-12
    assert abs(sum(torsion_terms.values()) - float(torsion_output["energy_eV"])) < 1.0e-12
    assert abs(float(torsion_output["energy_eV"]) + 13.114098449985388) < 1.0e-12

    torsion_eps = 1.0e-6
    torsion_plus = json.loads(json.dumps(torsion_case))
    torsion_minus = json.loads(json.dumps(torsion_case))
    torsion_plus["positions_A"][3][2] += torsion_eps
    torsion_minus["positions_A"][3][2] -= torsion_eps
    e_plus = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(torsion_plus)["energy_eV"])
    e_minus = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(torsion_minus)["energy_eV"])
    fd_force = -(e_plus - e_minus) / (2.0 * torsion_eps)
    analytic_force = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(torsion_case)["forces_eV_per_A"][3][2])
    assert abs(fd_force - analytic_force) < 1.0e-8

electron_case = json.loads(json.dumps(case))
electron_case["atom_types"] = ["C", "El"]
electron_case["positions_A"] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
electron_case["charge_constraints"]["fixed_charges"] = [0.0, 0.0]
electron_case["charge_in"] = [0.0, 0.0]
electron_case["solver_guess_in"] = [0.0, 0.0]
electron_case["crystx"] = None
electron_output = evaluate_reaxff_corrected_bond_nonbonded_valence_case(electron_case)
electron_terms = {key: float(value) for key, value in electron_output["energy_terms_eV"].items()}
assert abs(electron_terms["electron_core"] + 0.2326373860870282) < 1.0e-12
assert abs(electron_terms["electron_electron"]) < 1.0e-12
assert abs(sum(electron_terms.values()) - float(electron_output["energy_eV"])) < 1.0e-12
assert abs(float(electron_output["energy_eV"]) + 0.2249949320658194) < 1.0e-12

periodic_electron_case = json.loads(json.dumps(electron_case))
periodic_electron_case["atom_types"] = ["C", "El", "El"]
periodic_electron_case["positions_A"] = [
    [0.0, 0.0, 0.0],
    [0.2, 0.0, 0.0],
    [9.8, 0.0, 0.0],
]
periodic_electron_case["charge_constraints"]["fixed_charges"] = [0.0, 0.0, 0.0]
periodic_electron_case["charge_in"] = [0.0, 0.0, 0.0]
periodic_electron_case["solver_guess_in"] = [0.0, 0.0, 0.0]
periodic_electron_case["crystx"] = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
periodic_electron_output = evaluate_reaxff_corrected_bond_nonbonded_valence_case(periodic_electron_case)
periodic_electron_terms = {
    key: float(value)
    for key, value in periodic_electron_output["energy_terms_eV"].items()
}
assert abs(periodic_electron_terms["electron_core"] + 1.0956432714254907) < 1.0e-12
assert abs(periodic_electron_terms["electron_electron"] - 1.9507120300148208e-05) < 1.0e-12
assert periodic_electron_terms["electron_electron"] < 1.0e-3

hole_case = json.loads(json.dumps(electron_case))
hole_case["atom_types"] = ["C", "Ho"]
try:
    evaluate_reaxff_corrected_bond_nonbonded_valence_case(hole_case)
except ValueError as exc:
    assert "electron-only" in str(exc)
else:
    raise AssertionError("Ho/hole input should fail fast in the JAX backend")

electron_eps = 1.0e-6
electron_plus = json.loads(json.dumps(electron_case))
electron_minus = json.loads(json.dumps(electron_case))
electron_plus["positions_A"][1][0] += electron_eps
electron_minus["positions_A"][1][0] -= electron_eps
e_plus = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(electron_plus)["energy_eV"])
e_minus = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(electron_minus)["energy_eV"])
fd_force = -(e_plus - e_minus) / (2.0 * electron_eps)
analytic_force = float(evaluate_reaxff_corrected_bond_nonbonded_valence_case(electron_case)["forces_eV_per_A"][1][0])
assert abs(fd_force - analytic_force) < 1.0e-8
PY
else
  python3 projects/jamil_replication/fast_backend/jax_backend.py \
    --self-test \
    --dataset "$OUT_JSONL"
fi

echo "test_31_fast_backend_tools: PASS"
