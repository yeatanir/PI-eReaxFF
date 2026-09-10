#!/usr/bin/env python3
"""Compare current JAX-style nonbonded buckets against legacy oracle labels."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from jax_force_kernel import (
    KCAL_TO_EV,
    REAXFF_COULOMB_KCAL_A,
    _bond_energy_inputs,
    _orthorhombic_cell_lengths,
    _parameter_array,
    _taper_coefficients,
    parse_ffield_coulomb_parameters,
    taper_value,
)


def load_cases(path: Path, limit: int | None) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        labels = case.get("labels") or {}
        if "nonbonded_buckets" not in labels:
            raise ValueError(f"{path}:{line_number}: missing nonbonded_buckets labels")
        cases.append(case)
        if limit is not None and len(cases) >= limit:
            break
    return cases


def _atom_type_indices(case: dict[str, Any], parameters: dict[str, Any]) -> np.ndarray:
    index_by_symbol = {
        symbol: index
        for index, symbol in enumerate(parameters["amas_by_symbol"].keys(), start=1)
    }
    return np.asarray([index_by_symbol[symbol] for symbol in case["atom_types"]], dtype=int)


def _taper(r_A: float, coefficients: np.ndarray) -> float:
    return float(taper_value(np.asarray(r_A, dtype=float), coefficients))


def jax_style_buckets(case: dict[str, Any]) -> dict[tuple[int, int, bool], dict[str, float]]:
    parameters = parse_ffield_coulomb_parameters(Path(case["source_call_dir"]) / "ffield")
    bond_inputs = _bond_energy_inputs(case, parameters)
    positions = np.asarray(case["positions_A"], dtype=float)
    charges = np.asarray(case["labels"]["charge_out"], dtype=float)
    gamma_eem = np.asarray(_parameter_array(case, parameters["gamma_eem_by_symbol"], "gamma EEM"), dtype=float)
    gammavdw = np.asarray(_parameter_array(case, parameters["gammavdw_by_symbol"], "gamma vdW"), dtype=float)
    p1co = np.asarray(bond_inputs["p1co"], dtype=float)
    p2co = np.asarray(bond_inputs["p2co"], dtype=float)
    p3co = np.asarray(bond_inputs["p3co"], dtype=float)
    type_indices = _atom_type_indices(case, parameters)
    coefficients = np.asarray(_taper_coefficients(parameters["lower_taper_A"], parameters["upper_taper_A"]), dtype=float)
    upper = float(parameters["upper_taper_A"])
    power = float(parameters["vdw_shielding_power"])
    # Current Jamil oracle labels use central-cell nonbonded pairs.
    _orthorhombic_cell_lengths(case)
    offsets = np.zeros((1, 3))

    rows: dict[tuple[int, int, bool], dict[str, float]] = defaultdict(
        lambda: {"count": 0, "vdw_eV": 0.0, "coulomb_eV": 0.0}
    )
    n_atoms = positions.shape[0]
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            left = int(min(type_indices[i], type_indices[j]))
            right = int(max(type_indices[i], type_indices[j]))
            key = (left, right, False)
            for offset in offsets:
                disp = positions[i] - positions[j] + offset
                r = float(np.linalg.norm(disp))
                if r > upper or r < 1.0e-12:
                    continue
                rrw = (r**power + 1.0 / (np.sqrt(gammavdw[i] * gammavdw[j]) ** power)) ** (1.0 / power)
                h1 = np.exp(p3co[i, j] * (1.0 - rrw / p1co[i, j]))
                h2 = np.exp(0.5 * p3co[i, j] * (1.0 - rrw / p1co[i, j]))
                shielded = (r**3 + 1.0 / (np.sqrt(gamma_eem[i] * gamma_eem[j]) ** 3)) ** (1.0 / 3.0)
                sw = _taper(r, coefficients)
                rows[key]["count"] += 1
                rows[key]["vdw_eV"] += KCAL_TO_EV * p2co[i, j] * (h1 - 2.0 * h2) * sw
                rows[key]["coulomb_eV"] += KCAL_TO_EV * REAXFF_COULOMB_KCAL_A * charges[i] * charges[j] / shielded * sw
    return rows


def oracle_buckets(case: dict[str, Any]) -> dict[tuple[int, int, bool], dict[str, float]]:
    rows = {}
    for bucket in case["labels"]["nonbonded_buckets"]:
        key = (int(bucket["left_type"]), int(bucket["right_type"]), bool(bucket["self_image"]))
        rows[key] = {
            "count": int(bucket["count"]),
            "vdw_eV": float(bucket["vdw_eV"]),
            "coulomb_eV": float(bucket["coulomb_eV"]),
        }
    return rows


def summarize(cases: list[dict[str, Any]]) -> dict[str, Any]:
    stats: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "cases": 0,
            "oracle_vdw_eV": 0.0,
            "jax_vdw_eV": 0.0,
            "oracle_coulomb_eV": 0.0,
            "jax_coulomb_eV": 0.0,
            "max_vdw_abs_error_eV": 0.0,
            "max_coulomb_abs_error_eV": 0.0,
            "oracle_count": 0,
            "jax_count": 0,
        }
    )
    for case in cases:
        oracle = oracle_buckets(case)
        jax = jax_style_buckets(case)
        for key in sorted(set(oracle) | set(jax)):
            label = f"{key[0]}-{key[1]}:{'self' if key[2] else 'pair'}"
            o = oracle.get(key, {"count": 0, "vdw_eV": 0.0, "coulomb_eV": 0.0})
            j = jax.get(key, {"count": 0, "vdw_eV": 0.0, "coulomb_eV": 0.0})
            row = stats[label]
            row["cases"] += 1
            row["oracle_count"] += int(o["count"])
            row["jax_count"] += int(j["count"])
            row["oracle_vdw_eV"] += float(o["vdw_eV"])
            row["jax_vdw_eV"] += float(j["vdw_eV"])
            row["oracle_coulomb_eV"] += float(o["coulomb_eV"])
            row["jax_coulomb_eV"] += float(j["coulomb_eV"])
            row["max_vdw_abs_error_eV"] = max(row["max_vdw_abs_error_eV"], abs(float(j["vdw_eV"]) - float(o["vdw_eV"])))
            row["max_coulomb_abs_error_eV"] = max(row["max_coulomb_abs_error_eV"], abs(float(j["coulomb_eV"]) - float(o["coulomb_eV"])))

    by_bucket = {}
    for label, row in sorted(stats.items()):
        case_count = row["cases"]
        by_bucket[label] = {
            "mean_oracle_vdw_eV": row["oracle_vdw_eV"] / case_count,
            "mean_jax_vdw_eV": row["jax_vdw_eV"] / case_count,
            "mean_vdw_error_eV": (row["jax_vdw_eV"] - row["oracle_vdw_eV"]) / case_count,
            "max_vdw_abs_error_eV": row["max_vdw_abs_error_eV"],
            "mean_oracle_coulomb_eV": row["oracle_coulomb_eV"] / case_count,
            "mean_jax_coulomb_eV": row["jax_coulomb_eV"] / case_count,
            "mean_coulomb_error_eV": (row["jax_coulomb_eV"] - row["oracle_coulomb_eV"]) / case_count,
            "max_coulomb_abs_error_eV": row["max_coulomb_abs_error_eV"],
            "mean_oracle_count": row["oracle_count"] / case_count,
            "mean_jax_count": row["jax_count"] / case_count,
        }
    return {"cases": len(cases), "buckets": by_bucket}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = summarize(load_cases(args.dataset, args.limit))
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
