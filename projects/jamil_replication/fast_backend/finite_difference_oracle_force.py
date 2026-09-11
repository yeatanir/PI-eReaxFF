#!/usr/bin/env python3
"""Finite-difference exact-oracle energies for one atom coordinate."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from build_oracle_dataset import parse_forces, parse_geo
from jax_force_kernel import evaluate_force_batch


KCAL_TO_EV = 0.043364


def parse_summary_energy_eV(path: Path) -> float:
    rows = [line.split() for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if not rows:
        raise ValueError(f"{path}: missing summary row")
    return float(rows[-1][0]) * KCAL_TO_EV


def update_geo_coordinate(source: Path, target: Path, atom_index_1based: int, component: int, displacement_A: float) -> None:
    atom_cursor = 0
    out_lines = []
    for raw in source.read_text().splitlines():
        if raw.startswith("HETATM"):
            atom_cursor += 1
            if atom_cursor == atom_index_1based:
                values = [float(raw[30:40]), float(raw[40:50]), float(raw[50:60])]
                values[component] += displacement_A
                raw = f"{raw[:30]}{values[0]:10.5f}{values[1]:10.5f}{values[2]:10.5f}{raw[60:]}"
        out_lines.append(raw)
    if atom_cursor < atom_index_1based:
        raise ValueError(f"{source}: atom index {atom_index_1based} out of range")
    target.write_text("\n".join(out_lines) + "\n")


def run_oracle(binary: Path, source_call_dir: Path, work_dir: Path, atom_index_1based: int, component: int, displacement_A: float) -> dict[str, Any]:
    if work_dir.exists():
        shutil.rmtree(work_dir)
    shutil.copytree(source_call_dir, work_dir)
    update_geo_coordinate(source_call_dir / "geo", work_dir / "geo", atom_index_1based, component, displacement_A)
    for name in (
        "phase4_single_point_oracle.out",
        "single_point_forces.out",
        "single_point_energy_terms.out",
        "single_point_bond_orders.out",
        "single_point_nonbonded_buckets.out",
        "run_single_point.log",
    ):
        path = work_dir / name
        if path.exists():
            path.unlink()
    with (work_dir / "run_single_point.log").open("w") as log:
        proc = subprocess.run([str(binary)], cwd=work_dir, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"oracle failed in {work_dir}; see run_single_point.log")
    _, positions, _crystx, _molcharges = parse_geo(work_dir / "geo")
    forces, charges, solver = parse_forces(work_dir / "single_point_forces.out", len(positions))
    return {
        "work_dir": str(work_dir),
        "energy_eV": parse_summary_energy_eV(work_dir / "phase4_single_point_oracle.out"),
        "force_eV_per_A": forces[atom_index_1based - 1][component],
        "charge": charges[atom_index_1based - 1],
        "solver": solver[atom_index_1based - 1],
    }


def load_matching_case(dataset: Path, source_call_dir: Path) -> dict[str, Any] | None:
    wanted = source_call_dir.resolve()
    for line in dataset.read_text().splitlines():
        if not line.strip():
            continue
        case = json.loads(line)
        candidate = case.get("source_call_dir")
        if candidate is not None and Path(candidate).resolve() == wanted:
            return case
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--call-dir", type=Path, required=True)
    parser.add_argument("--atom", type=int, required=True, help="1-based atom index")
    parser.add_argument("--component", choices=("x", "y", "z"), required=True)
    parser.add_argument("--delta", type=float, default=1.0e-4)
    parser.add_argument("--work-root", type=Path, default=Path("build/jamil_replication/finite_difference_oracle_force"))
    parser.add_argument("--oracle-binary", type=Path, default=Path("build/oracle_single_point_force_export"))
    parser.add_argument("--dataset", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    component = {"x": 0, "y": 1, "z": 2}[args.component]
    call_dir = args.call_dir.resolve()
    work_root = (args.work_root / f"atom{args.atom}_{args.component}_d{args.delta:g}").resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    plus = run_oracle(args.oracle_binary.resolve(), call_dir, work_root / "plus", args.atom, component, args.delta)
    minus = run_oracle(args.oracle_binary.resolve(), call_dir, work_root / "minus", args.atom, component, -args.delta)
    central = run_oracle(args.oracle_binary.resolve(), call_dir, work_root / "central", args.atom, component, 0.0)

    derivative_eV_per_A = (plus["energy_eV"] - minus["energy_eV"]) / (2.0 * args.delta)
    finite_difference_force = -derivative_eV_per_A

    result: dict[str, Any] = {
        "call_dir": str(call_dir),
        "atom_index_1based": args.atom,
        "component": args.component,
        "delta_A": args.delta,
        "energy_plus_eV": plus["energy_eV"],
        "energy_minus_eV": minus["energy_eV"],
        "energy_central_eV": central["energy_eV"],
        "finite_difference_force_eV_per_A": finite_difference_force,
        "oracle_exported_force_eV_per_A": central["force_eV_per_A"],
        "fd_minus_oracle_exported_eV_per_A": finite_difference_force - central["force_eV_per_A"],
    }

    if args.dataset is not None:
        case = load_matching_case(args.dataset, call_dir)
        if case is None:
            raise ValueError(f"{args.dataset}: no case has source_call_dir={call_dir}")
        jax_output = evaluate_force_batch([case], kernel="reaxff_corrected_bond_nonbonded_valence")
        jax_force = float(np.asarray(jax_output["forces_eV_per_A"])[0, args.atom - 1, component])
        label_force = float(case["labels"]["forces_eV_per_A"][args.atom - 1][component])
        result.update(
            {
                "jax_force_eV_per_A": jax_force,
                "dataset_label_force_eV_per_A": label_force,
                "jax_minus_fd_eV_per_A": jax_force - finite_difference_force,
                "jax_minus_oracle_exported_eV_per_A": jax_force - central["force_eV_per_A"],
            }
        )

    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
