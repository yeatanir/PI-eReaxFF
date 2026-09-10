#!/usr/bin/env python3
"""Fail-closed entrypoint for the future JAX eReaxFF backend."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from jax_charge_kernel import charge_parity_metrics, evaluate_charge_batch


SCHEMA = "ereaxff_fast_backend_oracle_case_v1"


def load_cases(path: Path) -> list[dict[str, Any]]:
    cases = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        case = json.loads(line)
        if case.get("schema") != SCHEMA:
            raise ValueError(f"{path}:{line_number}: unexpected schema")
        cases.append(case)
    return cases


def jax_status() -> dict[str, Any]:
    try:
        import jax  # type: ignore

        jax.config.update("jax_enable_x64", True)
        return {
            "available": True,
            "version": getattr(jax, "__version__", "unknown"),
            "backend": jax.default_backend(),
            "x64": bool(jax.config.read("jax_enable_x64")),
        }
    except Exception as exc:  # pragma: no cover - depends on local optional deps.
        return {"available": False, "error": str(exc)}


def dataset_summary(cases: list[dict[str, Any]]) -> dict[str, Any]:
    labeled = [case for case in cases if case.get("labels") is not None]
    stages = sorted({case.get("stage", "unknown") for case in cases})
    beads = sorted({case.get("bead", -1) for case in cases})
    first = cases[0]
    return {
        "cases": len(cases),
        "labeled_cases": len(labeled),
        "labels_ready": len(labeled) == len(cases) and len(cases) > 0,
        "nparticles": first["nparticles"],
        "nbeads": first["nbeads"],
        "stages": stages,
        "beads": beads,
    }


def parity_check(dataset: Path) -> int:
    status = jax_status()
    if not status.get("available"):
        print(json.dumps({"jax": status, "parity_ready": False}, sort_keys=True))
        return 2

    import jax.numpy as jnp  # type: ignore

    cases = load_cases(dataset)
    summary = dataset_summary(cases)
    if not summary["labels_ready"]:
        print(json.dumps({"dataset": summary, "parity_ready": False}, sort_keys=True))
        return 2

    max_force_abs = 0.0
    max_charge_abs = 0.0
    max_solver_abs = 0.0
    max_charge_error = 0.0
    max_solver_error = 0.0
    energy_min = None
    energy_max = None
    charge_models = set()
    for case in cases:
        labels = case["labels"]
        charge_models.add(case.get("charge_model", "acks2"))
        forces = jnp.asarray(labels["forces_eV_per_A"], dtype=jnp.float64)
        charges = jnp.asarray(labels["charge_out"], dtype=jnp.float64)
        solver = jnp.asarray(labels["solver_guess_out"], dtype=jnp.float64)
        charge_metrics = charge_parity_metrics(case)
        max_force_abs = max(max_force_abs, float(jnp.max(jnp.abs(forces))))
        max_charge_abs = max(max_charge_abs, float(jnp.max(jnp.abs(charges))))
        max_solver_abs = max(max_solver_abs, float(jnp.max(jnp.abs(solver))))
        max_charge_error = max(max_charge_error, charge_metrics["max_charge_abs_error"])
        max_solver_error = max(max_solver_error, charge_metrics["max_solver_abs_error"])
        energy = float(labels["energy_eV"])
        energy_min = energy if energy_min is None else min(energy_min, energy)
        energy_max = energy if energy_max is None else max(energy_max, energy)

    charge_passes = max_charge_error <= 1.0e-12 and max_solver_error <= 1.0e-12
    print(
        json.dumps(
            {
                "dataset": str(dataset),
                "jax": status,
                "summary": summary,
                "oracle_label_stats": {
                    "energy_min_eV": energy_min,
                    "energy_max_eV": energy_max,
                    "max_force_abs_eV_per_A": max_force_abs,
                    "max_charge_abs": max_charge_abs,
                    "max_solver_abs": max_solver_abs,
                },
                "charge_kernel": {
                    "models": sorted(charge_models),
                    "max_charge_abs_error": max_charge_error,
                    "max_solver_abs_error": max_solver_error,
                    "passes": charge_passes,
                },
                "parity_ready": True,
            },
            sort_keys=True,
        )
    )
    return 0 if charge_passes else 2


def self_test(dataset: Path | None) -> int:
    status = jax_status()
    print(json.dumps({"jax": status}, sort_keys=True))
    if dataset is not None:
        cases = load_cases(dataset)
        if not cases:
            raise ValueError(f"{dataset} contains no cases")
        summary = dataset_summary(cases)
        summary["dataset"] = str(dataset)
        print(json.dumps(summary, sort_keys=True))
    return 0


def evaluate_from_stdin(allow_unvalidated: bool) -> int:
    request = json.load(sys.stdin)
    if request.get("command") == "charge_batch":
        batch = evaluate_charge_batch(request.get("cases", []))
        response = {
            "status": 0,
            "command": "charge_batch",
            "charge_models": batch["charge_models"],
            "charge_out": batch["charge_out"].tolist(),
            "solver_guess_out": batch["solver_guess_out"].tolist(),
            "case_status": batch["status"],
        }
        print(json.dumps(response, sort_keys=True))
        return 0
    response = {
        "status": 2,
        "diagnostic": "fast_jax_oracle force evaluation is not implemented or validated yet",
        "request_command": request.get("command"),
    }
    if allow_unvalidated:
        response["diagnostic"] = "unvalidated execution was requested, but no force kernel exists yet"
    print(json.dumps(response, sort_keys=True))
    return 2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, help="Oracle JSONL dataset for self-test/parity work")
    parser.add_argument("--self-test", action="store_true", help="Report optional JAX and dataset status")
    parser.add_argument("--parity-check", action="store_true", help="Validate labeled oracle cases with JAX arrays")
    parser.add_argument(
        "--allow-unvalidated",
        action="store_true",
        help="Acknowledge unvalidated mode; force evaluation still fails until a kernel exists",
    )
    args = parser.parse_args()

    if args.self_test:
        return self_test(args.dataset)
    if args.parity_check:
        if args.dataset is None:
            parser.error("--parity-check requires --dataset")
        return parity_check(args.dataset)
    return evaluate_from_stdin(args.allow_unvalidated)


if __name__ == "__main__":
    raise SystemExit(main())
