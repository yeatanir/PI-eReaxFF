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


def check_jsonl(path: Path) -> int:
    count = 0
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        if case.get("schema") != "ereaxff_fast_backend_oracle_case_v1":
            raise ValueError(f"{path}:{line_number}: unexpected schema")
        if len(case["positions_A"]) != case["nparticles"]:
            raise ValueError(f"{path}:{line_number}: position particle count mismatch")
        if len(case["charge_in"]) != case["nparticles"]:
            raise ValueError(f"{path}:{line_number}: charge particle count mismatch")
        if len(case["solver_guess_in"]) != case["nparticles"]:
            raise ValueError(f"{path}:{line_number}: solver particle count mismatch")
        count += 1
    print(f"checked {count} cases in {path}")
    return count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, help="Directory containing restart.dat")
    parser.add_argument("--output", type=Path, help="JSONL output path")
    parser.add_argument("--check", type=Path, help="Validate an existing JSONL dataset")
    args = parser.parse_args()

    if args.check is not None:
        check_jsonl(args.check)
        return 0
    if args.run_dir is None or args.output is None:
        parser.error("--run-dir and --output are required unless --check is used")

    cases = build_cases(args.run_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for case in cases:
            handle.write(json.dumps(case, sort_keys=True, separators=(",", ":")) + "\n")
    print(f"wrote {len(cases)} cases to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
