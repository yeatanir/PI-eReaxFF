#!/usr/bin/env python3
"""Diagnose JAX force contributions against labeled legacy oracle cases."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import jax
import jax.numpy as jnp

from jax_force_kernel import (
    _angle_parameter_tensors,
    _atom_type_index_array,
    _bond_energy_inputs,
    _coulomb_parameters_from_case,
    _electron_energy_inputs,
    _fixed_charges,
    _gamma_array,
    _orthorhombic_cell_lengths,
    _parameter_array,
    _sparse_topology_from_bond_orders,
    _taper_coefficients,
    _torsion_parameter_tensors,
    evaluate_force_batch,
    reaxff_corrected_bond_nonbonded_valence_terms,
)


def load_cases(path: Path, limit: int | None) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        labels = case.get("labels")
        if labels is None or "forces_eV_per_A" not in labels:
            raise ValueError(f"{path}:{line_number}: missing force labels")
        cases.append(case)
        if limit is not None and len(cases) >= limit:
            break
    return cases


def _case_term_forces(case: dict[str, Any]) -> dict[str, Any]:
    jax.config.update("jax_enable_x64", True)
    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    charges = _fixed_charges(case)
    cell_lengths = _orthorhombic_cell_lengths(case)
    parameters = _coulomb_parameters_from_case(case)
    gamma = _gamma_array(case, parameters)
    rvdw = _parameter_array(case, parameters["rvdw_by_symbol"], "Rvdw")
    eps = _parameter_array(case, parameters["eps_by_symbol"], "Evdw")
    alfa = _parameter_array(case, parameters["alfa_by_symbol"], "alfa")
    gammavdw = _parameter_array(case, parameters["gammavdw_by_symbol"], "gammavdW")
    coefficients = _taper_coefficients(parameters["lower_taper_A"], parameters["upper_taper_A"])
    upper_taper = float(parameters["upper_taper_A"])
    vdw_power = float(parameters["vdw_shielding_power"])
    bond_inputs = _bond_energy_inputs(case, parameters)
    atom_type_indices = _atom_type_index_array(case, parameters)
    angle_inputs = _angle_parameter_tensors(parameters)
    torsion_inputs = _torsion_parameter_tensors(parameters)
    electron_inputs = _electron_energy_inputs(case, parameters)
    topology = _sparse_topology_from_bond_orders(case, parameters, bond_inputs, cell_lengths)

    def terms_fn(x):
        return reaxff_corrected_bond_nonbonded_valence_terms(
            x,
            charges,
            gamma,
            rvdw,
            eps,
            alfa,
            gammavdw,
            vdw_power,
            coefficients,
            upper_taper,
            bond_inputs,
            atom_type_indices,
            angle_inputs,
            torsion_inputs,
            topology,
            electron_inputs,
            cell_lengths=cell_lengths,
        )

    terms = terms_fn(positions)
    term_forces = {}
    for name in sorted(terms):
        force = -jax.grad(lambda x, key=name: terms_fn(x)[key])(positions)
        term_forces[name] = force
    total_force = jnp.sum(jnp.stack(list(term_forces.values()), axis=0), axis=0)
    return {
        "terms": {key: float(value) for key, value in terms.items()},
        "term_forces": term_forces,
        "total_force": total_force,
    }


def _worst_component(error: jnp.ndarray) -> tuple[int, int, float]:
    abs_error = jnp.abs(error)
    flat_index = int(jnp.argmax(abs_error))
    atom, component = divmod(flat_index, 3)
    return atom, component, float(error[atom, component])


def summarize(cases: list[dict[str, Any]]) -> dict[str, Any]:
    batch = evaluate_force_batch(cases, kernel="reaxff_corrected_bond_nonbonded_valence")
    term_component_sums: dict[str, float] = defaultdict(float)
    term_component_max: dict[str, float] = defaultdict(float)
    leave_out_max: dict[str, float] = defaultdict(float)
    worst_cases = []
    global_max = 0.0

    for case_index, case in enumerate(cases):
        oracle = jnp.asarray(case["labels"]["forces_eV_per_A"], dtype=jnp.float64)
        jax_total = jnp.asarray(batch["forces_eV_per_A"][case_index], dtype=jnp.float64)
        error = jax_total - oracle
        atom, component, signed_error = _worst_component(error)
        max_abs_error = abs(signed_error)
        global_max = max(global_max, max_abs_error)
        term_report = _case_term_forces(case)
        component_terms = {}
        for name, forces in term_report["term_forces"].items():
            value = float(forces[atom, component])
            component_terms[name] = value
            term_component_sums[name] += abs(value)
            term_component_max[name] = max(term_component_max[name], abs(value))
            residual_without = term_report["total_force"] - forces - oracle
            leave_out_max[name] = max(leave_out_max[name], float(jnp.max(jnp.abs(residual_without))))
        worst_cases.append(
            {
                "case_index": case_index,
                "source_call_dir": case.get("source_call_dir"),
                "stage": case.get("stage"),
                "bead": case.get("bead"),
                "atom_index_0based": atom,
                "component": "xyz"[component],
                "jax_force_eV_per_A": float(jax_total[atom, component]),
                "oracle_force_eV_per_A": float(oracle[atom, component]),
                "signed_error_eV_per_A": signed_error,
                "abs_error_eV_per_A": max_abs_error,
                "jax_term_forces_at_component": dict(
                    sorted(component_terms.items(), key=lambda item: abs(item[1]), reverse=True)
                ),
            }
        )

    count = len(cases)
    return {
        "cases": count,
        "max_abs_error_eV_per_A": global_max,
        "largest_component_terms_mean_abs_eV_per_A": {
            key: value / count
            for key, value in sorted(term_component_sums.items(), key=lambda item: item[1], reverse=True)
        },
        "largest_component_terms_max_abs_eV_per_A": dict(
            sorted(term_component_max.items(), key=lambda item: item[1], reverse=True)
        ),
        "leave_one_term_out_max_abs_error_eV_per_A": dict(
            sorted(leave_out_max.items(), key=lambda item: item[1])
        ),
        "worst_cases": sorted(worst_cases, key=lambda row: row["abs_error_eV_per_A"], reverse=True)[:5],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = summarize(load_cases(args.dataset, args.limit))
    text = json.dumps(report, indent=2)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
