#!/usr/bin/env python3
"""Compare JAX corrected bond orders against legacy oracle bond-order labels."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from jax_force_kernel import evaluate_bond_order_state_case


def load_cases(path: Path, limit: int | None) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        labels = case.get("labels") or {}
        if "bond_orders" not in labels:
            raise ValueError(f"{path}:{line_number}: missing bond_orders labels")
        cases.append(case)
        if limit is not None and len(cases) >= limit:
            break
    return cases


def _update(row: dict[str, float], error: float) -> None:
    row["count"] += 1
    row["abs_error_sum"] += abs(error)
    row["signed_error_sum"] += error
    row["max_abs_error"] = max(row["max_abs_error"], abs(error))


def summarize(cases: list[dict[str, Any]], bond_threshold: float) -> dict[str, Any]:
    stats: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "count": 0,
            "abs_error_sum": 0.0,
            "signed_error_sum": 0.0,
            "max_abs_error": 0.0,
        }
    )
    missing_jax_bonds = 0
    extra_jax_bonds = 0
    max_case: dict[str, Any] | None = None

    for case_index, case in enumerate(cases):
        labels = case["labels"]["bond_orders"]
        output = evaluate_bond_order_state_case(case)
        abo = np.asarray(output["abo"], dtype=float)
        bo_matrix = np.asarray(output["bo_matrix"], dtype=float)
        bopi_matrix = np.asarray(output["bopi_matrix"], dtype=float)
        bopi2_matrix = np.asarray(output["bopi2_matrix"], dtype=float)
        active_matrix = np.asarray(output["active_matrix"], dtype=bool)

        for atom_index, oracle in enumerate(labels["abo"]):
            error = float(abo[atom_index]) - float(oracle)
            _update(stats["abo"], error)
            if max_case is None or abs(error) > max_case["abs_error"]:
                max_case = {
                    "case_index": case_index,
                    "call_slot": case.get("call_slot"),
                    "stage": case.get("stage"),
                    "kind": "abo",
                    "atom": atom_index + 1,
                    "oracle": float(oracle),
                    "jax": float(abo[atom_index]),
                    "error": error,
                    "abs_error": abs(error),
                }

        seen_pairs = set()
        for row in labels["bonds"]:
            i = int(row["atom_i"]) - 1
            j = int(row["atom_j"]) - 1
            meaningful = abs(row["bo"]) > bond_threshold or abs(row["bopi"]) > bond_threshold or abs(row["bopi2"]) > bond_threshold
            if meaningful:
                seen_pairs.add(tuple(sorted((i, j))))
            for key, matrix in (("bo", bo_matrix), ("bopi", bopi_matrix), ("bopi2", bopi2_matrix)):
                error = float(matrix[i, j]) - float(row[key])
                _update(stats[key], error)
                if max_case is None or abs(error) > max_case["abs_error"]:
                    max_case = {
                        "case_index": case_index,
                        "call_slot": case.get("call_slot"),
                        "stage": case.get("stage"),
                        "kind": key,
                        "atom_i": i + 1,
                        "atom_j": j + 1,
                        "oracle": float(row[key]),
                        "jax": float(matrix[i, j]),
                        "error": error,
                        "abs_error": abs(error),
                    }
            if meaningful and not active_matrix[i, j]:
                missing_jax_bonds += 1

        jax_pairs = {
            tuple(sorted((i, j)))
            for i, j in np.argwhere(np.triu(bo_matrix > bond_threshold, k=1))
        }
        extra_jax_bonds += len(jax_pairs - seen_pairs)

    rows = {}
    for name, row in sorted(stats.items()):
        count = int(row["count"])
        rows[name] = {
            "count": count,
            "mean_abs_error": row["abs_error_sum"] / count if count else 0.0,
            "mean_signed_error": row["signed_error_sum"] / count if count else 0.0,
            "max_abs_error": row["max_abs_error"],
        }
    return {
        "cases": len(cases),
        "bond_threshold": bond_threshold,
        "terms": rows,
        "missing_jax_bonds": missing_jax_bonds,
        "extra_jax_bonds": extra_jax_bonds,
        "max_error_case": max_case,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--bond-threshold", type=float, default=1.0e-8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = summarize(load_cases(args.dataset, args.limit), args.bond_threshold)
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
