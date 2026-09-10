#!/usr/bin/env python3
"""Persistent file-protocol worker for the experimental JAX eReaxFF force path."""

from __future__ import annotations

import argparse
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

from jax_force_kernel import evaluate_force_batch, parse_ffield_coulomb_parameters


def parse_request(path: Path) -> dict[str, Any]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    cursor = 0

    def take() -> list[str]:
        nonlocal cursor
        if cursor >= len(lines):
            raise ValueError("unexpected end of request")
        parts = lines[cursor].split()
        cursor += 1
        return parts

    key, value = take()
    if key != "request_id":
        raise ValueError("request must start with request_id")
    request_id = int(value)
    key, value = take()
    if key != "natoms":
        raise ValueError("request missing natoms")
    natoms = int(value)
    key, value = take()
    if key != "nbeads":
        raise ValueError("request missing nbeads")
    nbeads = int(value)
    parts = take()
    if parts[0] != "crystx" or len(parts) != 7:
        raise ValueError("request missing crystx")
    crystx = [float(x) for x in parts[1:]]
    if take()[0] != "atoms":
        raise ValueError("request missing atoms block")

    atom_types: list[str] = []
    fixed_charges: list[float] = []
    for expected in range(1, natoms + 1):
        parts = take()
        if int(parts[0]) != expected:
            raise ValueError("atom index mismatch")
        atom_types.append(parts[1])
        fixed_charges.append(float(parts[2]))

    positions_by_bead: list[list[list[float]]] = []
    charges_by_bead: list[list[float]] = []
    solver_by_bead: list[list[float]] = []
    for bead in range(1, nbeads + 1):
        parts = take()
        if parts[0] != "bead" or int(parts[1]) != bead:
            raise ValueError("bead block mismatch")
        positions: list[list[float]] = []
        charges: list[float] = []
        solver: list[float] = []
        for expected in range(1, natoms + 1):
            parts = take()
            if int(parts[0]) != expected:
                raise ValueError("bead atom index mismatch")
            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
            charges.append(float(parts[4]))
            solver.append(float(parts[5]))
        positions_by_bead.append(positions)
        charges_by_bead.append(charges)
        solver_by_bead.append(solver)

    return {
        "request_id": request_id,
        "natoms": natoms,
        "nbeads": nbeads,
        "crystx": crystx,
        "atom_types": atom_types,
        "fixed_charges": fixed_charges,
        "positions_by_bead": positions_by_bead,
        "charges_by_bead": charges_by_bead,
        "solver_by_bead": solver_by_bead,
    }


def write_response(path: Path, response: dict[str, Any]) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w") as fh:
        fh.write(f"request_id {response['request_id']}\n")
        fh.write(f"status {response['status']}\n")
        fh.write(f"diagnostic {response.get('diagnostic', 'ok')}\n")
        fh.write(f"nbeads {response.get('nbeads', 0)}\n")
        fh.write(f"natoms {response.get('natoms', 0)}\n")
        if response["status"] == 0:
            charges_by_bead = response["charges_by_bead"]
            solver_by_bead = response["solver_by_bead"]
            for bead_index, energy in enumerate(response["energy_eV"], start=1):
                fh.write(f"BEAD {bead_index} 0 {float(energy):.17e}\n")
                fh.write("CHARGES\n")
                for atom_index, charge in enumerate(charges_by_bead[bead_index - 1], start=1):
                    fh.write(f"{atom_index} {float(charge):.17e}\n")
                fh.write("SOLVER\n")
                for atom_index, solver in enumerate(solver_by_bead[bead_index - 1], start=1):
                    fh.write(f"{atom_index} {float(solver):.17e}\n")
                fh.write("FORCES\n")
                forces = response["forces_eV_per_A"][bead_index - 1]
                for atom_index, force in enumerate(forces, start=1):
                    fh.write(
                        f"{atom_index} {float(force[0]):.17e} {float(force[1]):.17e} {float(force[2]):.17e}\n"
                    )
                fh.write("END_BEAD\n")
        fh.write("END\n")
    tmp_path.replace(path)


def evaluate_request(request: dict[str, Any], parameters: dict[str, Any], kernel: str) -> dict[str, Any]:
    cases = []
    for positions, charges, solver in zip(
        request["positions_by_bead"], request["charges_by_bead"], request["solver_by_bead"], strict=True
    ):
        cases.append(
            {
                "atom_types": request["atom_types"],
                "positions_A": positions,
                "crystx": request["crystx"],
                "charge_constraints": {"fixed_charges": request["fixed_charges"]},
                "charge_in": charges,
                "solver_guess_in": solver,
                "coulomb_parameters": parameters,
            }
        )
    batch = evaluate_force_batch(cases, kernel=kernel)
    energy = np.asarray(batch["energy_eV"], dtype=float)
    forces = np.asarray(batch["forces_eV_per_A"], dtype=float)
    fixed = np.asarray(request["fixed_charges"], dtype=float)
    solver = np.asarray(request["solver_by_bead"], dtype=float)
    return {
        "request_id": request["request_id"],
        "status": 0,
        "diagnostic": "ok",
        "nbeads": request["nbeads"],
        "natoms": request["natoms"],
        "energy_eV": energy.tolist(),
        "forces_eV_per_A": forces.tolist(),
        "charges_by_bead": np.tile(fixed, (request["nbeads"], 1)).tolist(),
        "solver_by_bead": solver.tolist(),
    }


def run_worker(worker_root: Path, ffield: Path, kernel: str, poll_s: float) -> int:
    worker_root.mkdir(parents=True, exist_ok=True)
    request_path = worker_root / "request.in"
    request_ready = worker_root / "request.ready"
    response_path = worker_root / "response.out"
    response_ready = worker_root / "response.ready"
    shutdown_flag = worker_root / "shutdown.flag"
    done_path = worker_root / "worker.done"
    ready_path = worker_root / "worker.ready"
    parameters = parse_ffield_coulomb_parameters(ffield)
    ready_path.write_text("ready\n")
    try:
        while not shutdown_flag.exists():
            if not request_ready.exists():
                time.sleep(poll_s)
                continue
            request_ready.unlink(missing_ok=True)
            try:
                request = parse_request(request_path)
                response = evaluate_request(request, parameters, kernel)
            except Exception as exc:  # pragma: no cover - diagnostic path
                response = {
                    "request_id": -1,
                    "status": 2,
                    "diagnostic": f"{exc.__class__.__name__}:{exc}",
                    "nbeads": 0,
                    "natoms": 0,
                }
                (worker_root / "last_error.txt").write_text(traceback.format_exc())
            write_response(response_path, response)
            response_ready.write_text("ready\n")
    finally:
        done_path.write_text("done\n")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-root", required=True, type=Path)
    parser.add_argument("--ffield", required=True, type=Path)
    parser.add_argument("--kernel", default="reaxff_corrected_bond_nonbonded_valence")
    parser.add_argument("--poll-s", type=float, default=0.025)
    args = parser.parse_args()
    return run_worker(args.worker_root, args.ffield, args.kernel, args.poll_s)


if __name__ == "__main__":
    raise SystemExit(main())
