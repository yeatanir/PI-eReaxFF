#!/usr/bin/env python3
"""Build JSONL oracle cases from mixed BAOAB restart output.

This is a data-extraction tool, not a fast force implementation.  It creates
the stable reference cases that a JAX/native backend must match before it is
allowed to drive dynamics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

KCAL_TO_EV = 0.043364
SCHEMA = "ereaxff_fast_backend_oracle_case_v1"
ENERGY_TERM_NAMES = [
    "bond",
    "overcoordination",
    "lone_pair",
    "molecule",
    "valence_penalty",
    "angle_conjugation",
    "hydrogen_bond",
    "torsion",
    "torsion_conjugation",
    "vdw",
    "coulomb",
    "charge",
    "external_field",
    "electron",
]


def looks_float(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def parse_restart(path: Path) -> dict[str, Any]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or not lines[0].startswith("# PI-eReaxFF mixed BAOAB restart"):
        raise ValueError(f"{path} is not a mixed BAOAB restart")

    idx = 1
    header: dict[str, Any] = {}
    for key in ("step", "ndim", "nparticles", "nquantum", "beads"):
        got, value = lines[idx].split(maxsplit=1)
        if got != key:
            raise ValueError(f"expected {key}, got {got}")
        header[key] = int(value)
        idx += 1
    got, value = lines[idx].split(maxsplit=1)
    if got != "charge_policy":
        raise ValueError(f"expected charge_policy, got {got}")
    header["charge_policy"] = value
    idx += 1
    for key in ("classical_heat_eV", "rng_state_i8"):
        got, value = lines[idx].split(maxsplit=1)
        if got != key:
            raise ValueError(f"expected {key}, got {got}")
        header[key] = float(value) if key.endswith("_eV") else int(value)
        idx += 1

    ndim = header["ndim"]
    nparticles = header["nparticles"]
    beads = header["beads"]

    quantum: list[dict[str, Any]] = []
    if lines[idx] != "quantum_indices":
        raise ValueError("expected quantum_indices")
    idx += 1
    for _ in range(header["nquantum"]):
        q, particle, qbeads, mass, heat = lines[idx].split()
        quantum.append(
            {
                "q": int(q),
                "particle_index": int(particle),
                "beads": int(qbeads),
                "mass_amu": float(mass),
                "heat_eV": float(heat),
            }
        )
        idx += 1

    positions: list[list[list[float]]] = []
    if lines[idx] != "positions_A":
        raise ValueError("expected positions_A")
    idx += 1
    for bead in range(1, beads + 1):
        key, value = lines[idx].split()
        if key != "bead" or int(value) != bead:
            raise ValueError("position bead mismatch")
        idx += 1
        bead_positions = []
        for particle in range(1, nparticles + 1):
            values = lines[idx].split()
            if int(values[0]) != particle:
                raise ValueError("position particle mismatch")
            bead_positions.append([float(x) for x in values[1 : 1 + ndim]])
            idx += 1
        positions.append(bead_positions)

    if lines[idx] != "classical_velocities":
        raise ValueError("expected classical_velocities")
    idx += 1 + nparticles

    for block in ("q_nm", "v_nm"):
        if lines[idx] != block:
            raise ValueError(f"expected {block}")
        idx += 1
        for qinfo in quantum:
            if lines[idx] != f"quantum {qinfo['q']}":
                raise ValueError(f"{block} quantum mismatch")
            idx += 1
            idx += 1
            for coord in range(1, ndim + 1):
                if lines[idx] != f"coord {coord}":
                    raise ValueError(f"{block} coord mismatch")
                idx += 1 + beads

    charges_by_bead: list[list[float]] = []
    solver_by_bead: list[list[float]] = []
    if lines[idx] != "charges_by_bead":
        raise ValueError("expected charges_by_bead")
    idx += 1
    for bead in range(1, beads + 1):
        key, value = lines[idx].split()
        if key != "bead" or int(value) != bead:
            raise ValueError("charge bead mismatch")
        idx += 1
        bead_charges = []
        bead_solver = []
        for particle in range(1, nparticles + 1):
            values = lines[idx].split()
            if int(values[0]) != particle:
                raise ValueError("charge particle mismatch")
            bead_charges.append(float(values[1]))
            bead_solver.append(float(values[2]))
            idx += 1
        charges_by_bead.append(bead_charges)
        solver_by_bead.append(bead_solver)

    return {
        "header": header,
        "quantum": quantum,
        "positions_A_by_bead": positions,
        "charges_by_bead": charges_by_bead,
        "solver_guess_by_bead": solver_by_bead,
    }


def read_energy_tail(path: Path) -> dict[str, float] | None:
    if not path.exists():
        return None
    rows = [
        line.split()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not rows:
        return None
    last = rows[-1]
    return {
        "step": float(last[0]),
        "kinetic_eV": float(last[1]),
        "physical_potential_eV": float(last[2]),
        "spring_eV": float(last[3]),
        "conserved_like_eV": float(last[4]),
        "charge_diff": float(last[5]),
        "max_force_eV_per_A": float(last[6]),
    }


def build_cases(run_dir: Path) -> list[dict[str, Any]]:
    restart = parse_restart(run_dir / "restart.dat")
    energy_tail = read_energy_tail(run_dir / "pimd_energy.out")
    header = restart["header"]
    cases = []
    for bead_index, positions in enumerate(restart["positions_A_by_bead"], start=1):
        cases.append(
            {
                "schema": "ereaxff_fast_backend_oracle_case_v1",
                "source_run_dir": str(run_dir),
                "step": header["step"],
                "bead": bead_index,
                "stage": "restart",
                "call_slot": None,
                "ndim": header["ndim"],
                "nparticles": header["nparticles"],
                "nbeads": header["beads"],
                "charge_policy": header["charge_policy"],
                "quantum": restart["quantum"],
                "positions_A": positions,
                "charge_in": restart["charges_by_bead"][bead_index - 1],
                "solver_guess_in": restart["solver_guess_by_bead"][bead_index - 1],
                "run_energy_tail": energy_tail,
                "labels": None,
            }
        )
    return cases


def parse_timing_calls(adapter_work: Path) -> list[dict[str, Any]]:
    timing_path = adapter_work / "persistent_oracle_timing.log"
    if not timing_path.exists():
        timing_path = adapter_work / "oracle_timing.log"
    if not timing_path.exists():
        return []

    calls: list[dict[str, Any]] = []
    for line in timing_path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if fields[0] != "CALL":
            continue
        calls.append(
            {
                "step": int(fields[1]),
                "batch": int(fields[2]),
                "call_slot": int(fields[3]),
                "bead": int(fields[4]),
                "worker": int(fields[5]),
                "dispatch_mode": fields[6],
                "stage_kind": fields[7],
                "stage": fields[8],
                "worker_state": fields[9],
                "elapsed_s": float(fields[12]),
                "status": int(fields[13]),
            }
        )
    return calls


def find_call_dir(adapter_work: Path, call_slot: int) -> Path:
    name = f"adapter_call_{call_slot:04d}"
    direct = adapter_work / name
    if direct.is_dir():
        return direct
    matches = sorted(adapter_work.glob(f"**/{name}"))
    matches = [path for path in matches if path.is_dir()]
    if not matches:
        raise FileNotFoundError(f"could not find {name} under {adapter_work}")
    if len(matches) > 1:
        matches = sorted(matches, key=lambda path: len(path.parts))
    return matches[0]


def parse_control(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[1] == "icharg":
            values["icharg"] = int(fields[0])
    return values


def parse_geo(path: Path) -> tuple[list[str], list[list[float]], list[float] | None, list[dict[str, Any]]]:
    atom_types: list[str] = []
    positions: list[list[float]] = []
    crystx: list[float] | None = None
    molcharges: list[dict[str, Any]] = []
    for line in path.read_text().splitlines():
        if line.startswith("CRYSTX"):
            crystx = [float(value) for value in line.split()[1:7]]
        elif line.startswith("HETATM"):
            fields = line.split()
            if len(fields) < 7:
                raise ValueError(f"{path}: malformed HETATM line")
            atom_types.append(fields[3])
            positions.append([float(fields[4]), float(fields[5]), float(fields[6])])
        elif line.startswith("MOLCHARGE"):
            fields = line.split()
            if len(fields) != 4:
                raise ValueError(f"{path}: malformed MOLCHARGE line")
            molcharges.append(
                {
                    "start_atom": int(fields[1]),
                    "end_atom": int(fields[2]),
                    "charge": float(fields[3]),
                }
            )
    if not positions:
        raise ValueError(f"{path} contains no HETATM positions")
    return atom_types, positions, crystx, molcharges


def fixed_charges_from_molcharges(nparticles: int, molcharges: list[dict[str, Any]]) -> list[float] | None:
    if not molcharges:
        return None
    charges = [0.0] * nparticles
    touched = [False] * nparticles
    for group in molcharges:
        start = group["start_atom"]
        end = group["end_atom"]
        if start < 1 or end > nparticles or start > end:
            raise ValueError(f"invalid MOLCHARGE range {start}:{end} for {nparticles} atoms")
        per_atom = group["charge"] / float(end - start + 1)
        for index in range(start - 1, end):
            charges[index] = per_atom
            touched[index] = True
    if not all(touched):
        return None
    return charges


def parse_charge_memory(path: Path) -> tuple[list[float], list[float]]:
    charge_in: list[float] = []
    solver_in: list[float] = []
    for expected, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split()
        atom_index = int(fields[0])
        if atom_index != expected:
            raise ValueError(f"{path}: expected atom {expected}, got {atom_index}")
        charge_in.append(float(fields[1]))
        solver_in.append(float(fields[2]))
    return charge_in, solver_in


def parse_summary(path: Path) -> dict[str, Any]:
    rows = [
        line.split()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one summary row")
    row = rows[0]
    energy_kcal = float(row[0])
    return {
        "energy_kcal": energy_kcal,
        "energy_eV": energy_kcal * KCAL_TO_EV,
        "max_force": float(row[1]),
        "charge_sum": float(row[2]),
        "status": int(row[3]),
        "nparticles": int(row[4]),
    }


def parse_energy_terms(path: Path) -> tuple[dict[str, float], dict[str, float]] | None:
    if not path.exists():
        return None
    rows = [
        line.split()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not rows:
        return None

    first_tokens = [row[0] for row in rows if row]
    if first_tokens and not any(looks_float(token) for token in first_tokens):
        terms_kcal = {}
        for row in rows:
            if len(row) < 2:
                raise ValueError(f"{path}: malformed named energy-term row")
            name = row[0]
            if name == "total":
                continue
            if name not in ENERGY_TERM_NAMES:
                raise ValueError(f"{path}: unexpected energy term {name!r}")
            terms_kcal[name] = float(row[1])
    else:
        row = rows[-1]
        values = [float(value) for value in row[1:]]
        if len(values) < len(ENERGY_TERM_NAMES):
            raise ValueError(f"{path}: expected at least {len(ENERGY_TERM_NAMES)} term values, got {len(values)}")
        terms_kcal = dict(zip(ENERGY_TERM_NAMES, values[: len(ENERGY_TERM_NAMES)]))

    terms_eV = {name: value * KCAL_TO_EV for name, value in terms_kcal.items()}
    return terms_kcal, terms_eV


def parse_optional_energy_terms(call_dir: Path) -> tuple[dict[str, float], dict[str, float]] | None:
    named = parse_energy_terms(call_dir / "single_point_energy_terms.out")
    if named is not None:
        return named
    return parse_energy_terms(call_dir / "fort.73")


def parse_bond_orders(path: Path, nparticles: int) -> dict[str, Any] | None:
    if not path.exists():
        return None
    abo = [None] * nparticles
    bonds: list[dict[str, Any]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        kind = fields[0]
        if kind == "ATOM":
            if len(fields) != 3:
                raise ValueError(f"{path}: malformed ATOM bond-order row")
            atom_index = int(fields[1])
            if atom_index < 1 or atom_index > nparticles:
                raise ValueError(f"{path}: ATOM index out of range")
            abo[atom_index - 1] = float(fields[2])
        elif kind == "BOND":
            if len(fields) != 7:
                raise ValueError(f"{path}: malformed BOND bond-order row")
            bonds.append(
                {
                    "index": int(fields[1]),
                    "atom_i": int(fields[2]),
                    "atom_j": int(fields[3]),
                    "bo": float(fields[4]),
                    "bopi": float(fields[5]),
                    "bopi2": float(fields[6]),
                }
            )
        else:
            raise ValueError(f"{path}: unexpected bond-order row kind {kind!r}")
    if any(value is None for value in abo):
        raise ValueError(f"{path}: missing ATOM bond-order rows")
    for expected, row in enumerate(bonds, start=1):
        if row["index"] != expected:
            raise ValueError(f"{path}: expected BOND index {expected}, got {row['index']}")
        if row["atom_i"] < 1 or row["atom_i"] > nparticles or row["atom_j"] < 1 or row["atom_j"] > nparticles:
            raise ValueError(f"{path}: BOND atom index out of range")
    return {
        "abo": [float(value) for value in abo],
        "bonds": bonds,
    }


def parse_forces(path: Path, nparticles: int) -> tuple[list[list[float]], list[float], list[float]]:
    forces: list[list[float]] = []
    charges: list[float] = []
    solver: list[float] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        atom_index = int(fields[0])
        expected = len(forces) + 1
        if atom_index != expected:
            raise ValueError(f"{path}: expected atom {expected}, got {atom_index}")
        forces.append([float(fields[1]), float(fields[2]), float(fields[3])])
        charges.append(float(fields[4]))
        solver.append(float(fields[5]))
    if len(forces) != nparticles:
        raise ValueError(f"{path}: expected {nparticles} atoms, got {len(forces)}")
    return forces, charges, solver


def build_labeled_cases(adapter_work: Path) -> list[dict[str, Any]]:
    calls = parse_timing_calls(adapter_work)
    if not calls:
        slots = []
        for path in sorted(adapter_work.glob("**/adapter_call_*")):
            if path.is_dir():
                slots.append(int(path.name.rsplit("_", 1)[1]))
        calls = [
            {
                "step": -1,
                "batch": -1,
                "call_slot": slot,
                "bead": -1,
                "worker": -1,
                "dispatch_mode": "unknown",
                "stage_kind": "unknown",
                "stage": "unknown",
                "worker_state": "unknown",
                "elapsed_s": 0.0,
                "status": 0,
            }
            for slot in sorted(slots)
        ]

    cases: list[dict[str, Any]] = []
    for call in sorted(calls, key=lambda item: item["call_slot"]):
        call_dir = find_call_dir(adapter_work, call["call_slot"])
        atom_types, positions, crystx, molcharges = parse_geo(call_dir / "geo")
        control = parse_control(call_dir / "control")
        charge_in, solver_in = parse_charge_memory(call_dir / "charge_memory.in")
        summary = parse_summary(call_dir / "phase4_single_point_oracle.out")
        energy_terms = parse_optional_energy_terms(call_dir)
        forces, charge_out, solver_out = parse_forces(call_dir / "single_point_forces.out", summary["nparticles"])
        nparticles = summary["nparticles"]
        bond_orders = parse_bond_orders(call_dir / "single_point_bond_orders.out", nparticles)
        fixed_charges = fixed_charges_from_molcharges(nparticles, molcharges)
        charge_model = "fixed_molcharge" if control.get("icharg") == 8 and fixed_charges is not None else "acks2"
        if len(positions) != nparticles or len(charge_in) != nparticles or len(solver_in) != nparticles:
            raise ValueError(f"{call_dir}: inconsistent particle counts")
        labels = {
            "energy_kcal": summary["energy_kcal"],
            "energy_eV": summary["energy_eV"],
            "forces_eV_per_A": forces,
            "charge_out": charge_out,
            "solver_guess_out": solver_out,
            "max_force": summary["max_force"],
            "charge_sum": summary["charge_sum"],
            "status": summary["status"],
        }
        if energy_terms is not None:
            terms_kcal, terms_eV = energy_terms
            labels["energy_terms_kcal"] = terms_kcal
            labels["energy_terms_eV"] = terms_eV
        if bond_orders is not None:
            labels["bond_orders"] = bond_orders

        cases.append(
            {
                "schema": SCHEMA,
                "source_adapter_work": str(adapter_work),
                "source_call_dir": str(call_dir),
                "step": call["step"],
                "bead": call["bead"],
                "stage": call["stage"],
                "stage_kind": call["stage_kind"],
                "call_slot": call["call_slot"],
                "worker": call["worker"],
                "dispatch_mode": call["dispatch_mode"],
                "elapsed_s": call["elapsed_s"],
                "ndim": 3,
                "nparticles": nparticles,
                "nbeads": max(1, max(item["bead"] for item in calls)),
                "atom_types": atom_types,
                "crystx": crystx,
                "control": control,
                "charge_model": charge_model,
                "charge_constraints": {
                    "molcharges": molcharges,
                    "fixed_charges": fixed_charges,
                },
                "charge_policy": "bead_local_acks2",
                "quantum": [],
                "positions_A": positions,
                "charge_in": charge_in,
                "solver_guess_in": solver_in,
                "run_energy_tail": None,
                "labels": labels,
            }
        )
    return cases


def check_jsonl(path: Path) -> int:
    count = 0
    labeled = 0
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        if case.get("schema") != SCHEMA:
            raise ValueError(f"{path}:{line_number}: unexpected schema")
        if len(case["positions_A"]) != case["nparticles"]:
            raise ValueError(f"{path}:{line_number}: position particle count mismatch")
        if len(case["charge_in"]) != case["nparticles"]:
            raise ValueError(f"{path}:{line_number}: charge particle count mismatch")
        if len(case["solver_guess_in"]) != case["nparticles"]:
            raise ValueError(f"{path}:{line_number}: solver particle count mismatch")
        labels = case.get("labels")
        if labels is not None:
            if len(labels["forces_eV_per_A"]) != case["nparticles"]:
                raise ValueError(f"{path}:{line_number}: force label particle count mismatch")
            if len(labels["charge_out"]) != case["nparticles"]:
                raise ValueError(f"{path}:{line_number}: charge label particle count mismatch")
            if len(labels["solver_guess_out"]) != case["nparticles"]:
                raise ValueError(f"{path}:{line_number}: solver label particle count mismatch")
            if "energy_terms_eV" in labels:
                energy_terms = labels["energy_terms_eV"]
                if not isinstance(energy_terms, dict) or not energy_terms:
                    raise ValueError(f"{path}:{line_number}: invalid energy_terms_eV label")
                for name, value in energy_terms.items():
                    if name not in ENERGY_TERM_NAMES:
                        raise ValueError(f"{path}:{line_number}: unexpected energy term {name!r}")
                    float(value)
            if "bond_orders" in labels:
                bond_orders = labels["bond_orders"]
                if len(bond_orders["abo"]) != case["nparticles"]:
                    raise ValueError(f"{path}:{line_number}: bond-order abo particle count mismatch")
                for value in bond_orders["abo"]:
                    float(value)
                previous_index = 0
                for row in bond_orders["bonds"]:
                    previous_index += 1
                    if row["index"] != previous_index:
                        raise ValueError(f"{path}:{line_number}: bond-order index mismatch")
                    if row["atom_i"] < 1 or row["atom_i"] > case["nparticles"]:
                        raise ValueError(f"{path}:{line_number}: bond-order atom_i out of range")
                    if row["atom_j"] < 1 or row["atom_j"] > case["nparticles"]:
                        raise ValueError(f"{path}:{line_number}: bond-order atom_j out of range")
                    float(row["bo"])
                    float(row["bopi"])
                    float(row["bopi2"])
            labeled += 1
        count += 1
    print(f"checked {count} cases in {path}; labeled={labeled}; labels_ready={labeled == count and count > 0}")
    return count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, help="Directory containing restart.dat")
    parser.add_argument("--adapter-work", type=Path, help="Adapter work directory containing oracle call outputs")
    parser.add_argument("--output", type=Path, help="JSONL output path")
    parser.add_argument("--check", type=Path, help="Validate an existing JSONL dataset")
    args = parser.parse_args()

    if args.check is not None:
        check_jsonl(args.check)
        return 0
    if args.output is None:
        parser.error("--output is required unless --check is used")
    if args.adapter_work is None and args.run_dir is None:
        parser.error("--adapter-work or --run-dir is required unless --check is used")

    if args.adapter_work is not None:
        cases = build_labeled_cases(args.adapter_work)
    else:
        cases = build_cases(args.run_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for case in cases:
            handle.write(json.dumps(case, sort_keys=True, separators=(",", ":")) + "\n")
    print(f"wrote {len(cases)} cases to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
