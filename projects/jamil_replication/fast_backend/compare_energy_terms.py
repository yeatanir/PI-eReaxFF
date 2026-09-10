#!/usr/bin/env python3
"""Compare unvalidated JAX energy terms against labeled legacy oracle cases."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from jax_force_kernel import evaluate_force_batch


LEGACY_TO_JAX = {
    "bond": "bond",
    "overcoordination": (),
    "lone_pair": (),
    "molecule": (),
    "valence_penalty": ("valence", "penalty"),
    "angle_conjugation": "angle_conjugation",
    "hydrogen_bond": (),
    "torsion": "torsion",
    "torsion_conjugation": "torsion_conjugation",
    "vdw": "vdw",
    "coulomb": "coulomb",
    "charge": (),
    "external_field": (),
    "electron": ("electron_core", "electron_electron"),
}


def load_cases(path: Path, limit: int | None) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        labels = case.get("labels")
        if labels is None or "energy_terms_eV" not in labels:
            raise ValueError(f"{path}:{line_number}: missing energy_terms_eV labels")
        cases.append(case)
        if limit is not None and len(cases) >= limit:
            break
    return cases


def jax_term_value(output: dict[str, Any], mapping: str | tuple[str, ...]) -> float:
    terms = output.get("energy_terms_eV", {})
    keys = (mapping,) if isinstance(mapping, str) else mapping
    return sum(float(terms.get(key, 0.0)) for key in keys)


def summarize(cases: list[dict[str, Any]], outputs: list[dict[str, Any]]) -> dict[str, Any]:
    stats: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "count": 0,
            "oracle_sum_eV": 0.0,
            "jax_sum_eV": 0.0,
            "signed_error_sum_eV": 0.0,
            "abs_error_sum_eV": 0.0,
            "max_abs_error_eV": 0.0,
        }
    )
    total_abs_error_sum = 0.0
    total_max_abs_error = 0.0
    force_max_abs = 0.0

    for case, output in zip(cases, outputs, strict=True):
        labels = case["labels"]
        total_error = float(output["energy_eV"]) - float(labels["energy_eV"])
        total_abs_error_sum += abs(total_error)
        total_max_abs_error = max(total_max_abs_error, abs(total_error))
        oracle_forces = labels["forces_eV_per_A"]
        jax_forces = output["forces_eV_per_A"]
        for atom_oracle, atom_jax in zip(oracle_forces, jax_forces, strict=True):
            for oracle_value, jax_value in zip(atom_oracle, atom_jax, strict=True):
                force_max_abs = max(force_max_abs, abs(float(jax_value) - float(oracle_value)))
        for legacy_name, mapping in LEGACY_TO_JAX.items():
            oracle = float(labels["energy_terms_eV"].get(legacy_name, 0.0))
            jax = jax_term_value(output, mapping)
            error = jax - oracle
            row = stats[legacy_name]
            row["count"] += 1
            row["oracle_sum_eV"] += oracle
            row["jax_sum_eV"] += jax
            row["signed_error_sum_eV"] += error
            row["abs_error_sum_eV"] += abs(error)
            row["max_abs_error_eV"] = max(row["max_abs_error_eV"], abs(error))

    term_rows = {}
    for name, row in sorted(stats.items()):
        count = row["count"]
        term_rows[name] = {
            "mean_oracle_eV": row["oracle_sum_eV"] / count,
            "mean_jax_eV": row["jax_sum_eV"] / count,
            "mean_signed_error_eV": row["signed_error_sum_eV"] / count,
            "mean_abs_error_eV": row["abs_error_sum_eV"] / count,
            "max_abs_error_eV": row["max_abs_error_eV"],
        }

    return {
        "cases": len(cases),
        "kernel": outputs[0].get("kernel", "unknown") if outputs else "unknown",
        "total_energy": {
            "mean_abs_error_eV": total_abs_error_sum / len(cases),
            "max_abs_error_eV": total_max_abs_error,
        },
        "forces": {
            "max_abs_error_eV_per_A": force_max_abs,
        },
        "terms": term_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--kernel", default="reaxff_corrected_bond_nonbonded_valence")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    cases = load_cases(args.dataset, args.limit)
    if not cases:
        raise ValueError(f"{args.dataset} contains no labeled cases")
    batch = evaluate_force_batch(cases, kernel=args.kernel)
    kernels = batch.get("kernel", [])
    if len(kernels) == 1 and len(cases) > 1:
        kernels = kernels * len(cases)
    outputs = []
    for index, case in enumerate(cases):
        outputs.append(
            {
                "kernel": kernels[index],
                "energy_eV": float(batch["energy_eV"][index]),
                "energy_terms_eV": {
                    key: float(values[index])
                    for key, values in batch.get("energy_terms_eV", {}).items()
                },
                "forces_eV_per_A": batch["forces_eV_per_A"][index].tolist(),
                "status": batch["status"][index],
                "case": {
                    "source_call_dir": case.get("source_call_dir"),
                    "stage": case.get("stage"),
                    "bead": case.get("bead"),
                },
            }
        )
    report = summarize(cases, outputs)
    report["dataset"] = str(args.dataset)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
