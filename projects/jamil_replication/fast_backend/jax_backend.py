#!/usr/bin/env python3
"""Fail-closed entrypoint for the future JAX eReaxFF backend."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


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


def self_test(dataset: Path | None) -> int:
    status = jax_status()
    print(json.dumps({"jax": status}, sort_keys=True))
    if dataset is not None:
        cases = load_cases(dataset)
        if not cases:
            raise ValueError(f"{dataset} contains no cases")
        first = cases[0]
        print(
            json.dumps(
                {
                    "dataset": str(dataset),
                    "cases": len(cases),
                    "nparticles": first["nparticles"],
                    "nbeads": first["nbeads"],
                    "labels_ready": all(case.get("labels") is not None for case in cases),
                },
                sort_keys=True,
            )
        )
    return 0


def evaluate_from_stdin(allow_unvalidated: bool) -> int:
    request = json.load(sys.stdin)
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
    parser.add_argument(
        "--allow-unvalidated",
        action="store_true",
        help="Acknowledge unvalidated mode; force evaluation still fails until a kernel exists",
    )
    args = parser.parse_args()

    if args.self_test:
        return self_test(args.dataset)
    return evaluate_from_stdin(args.allow_unvalidated)


if __name__ == "__main__":
    raise SystemExit(main())

