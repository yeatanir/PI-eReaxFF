"""JAX force-kernel interfaces for the fast eReaxFF backend."""

from __future__ import annotations

import hashlib
import json
from typing import Any


def geometry_fingerprint(case: dict[str, Any]) -> str:
    payload = {
        "atom_types": case.get("atom_types"),
        "crystx": case.get("crystx"),
        "positions_A": case.get("positions_A"),
        "charge_in": case.get("charge_in"),
        "solver_guess_in": case.get("solver_guess_in"),
        "charge_model": case.get("charge_model"),
        "charge_constraints": case.get("charge_constraints"),
    }
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def evaluate_force_case(case: dict[str, Any], *, allow_replay: bool = False) -> dict[str, Any]:
    """Evaluate one force case.

    The current implementation supports only strict oracle-label replay.  That
    is useful for validating batching, protocol shape, and Fortran integration
    before a physics kernel is enabled.  Unknown geometries fail closed.
    """

    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    labels = case.get("labels")
    if labels is None:
        raise ValueError("force evaluation requires labels until a physics kernel exists")
    if not allow_replay:
        raise NotImplementedError("JAX eReaxFF physics force kernel is not implemented yet")
    forces = jnp.asarray(labels["forces_eV_per_A"], dtype=jnp.float64)
    return {
        "kernel": "oracle_replay",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": jnp.asarray(labels["energy_eV"], dtype=jnp.float64),
        "forces_eV_per_A": forces,
        "status": int(labels["status"]),
    }


def evaluate_force_batch(cases: list[dict[str, Any]], *, allow_replay: bool = False) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    if not cases:
        raise ValueError("force batch is empty")
    outputs = [evaluate_force_case(case, allow_replay=allow_replay) for case in cases]
    energies = jnp.stack([output["energy_eV"] for output in outputs], axis=0)
    forces = jnp.stack([output["forces_eV_per_A"] for output in outputs], axis=0)
    identity = jax.jit(lambda energy_array, force_array: (energy_array, force_array))
    energies, forces = identity(energies, forces)
    return {
        "kernel": sorted({output["kernel"] for output in outputs}),
        "fingerprints": [str(output["fingerprint"]) for output in outputs],
        "energy_eV": energies,
        "forces_eV_per_A": forces,
        "status": [int(output["status"]) for output in outputs],
    }


def force_parity_metrics(case: dict[str, Any], *, allow_replay: bool = False) -> dict[str, float]:
    import jax.numpy as jnp

    labels = case.get("labels")
    if labels is None:
        raise ValueError("force parity requires labels")
    predicted = evaluate_force_case(case, allow_replay=allow_replay)
    energy_ref = jnp.asarray(labels["energy_eV"], dtype=jnp.float64)
    force_ref = jnp.asarray(labels["forces_eV_per_A"], dtype=jnp.float64)
    energy_err = predicted["energy_eV"] - energy_ref
    force_err = predicted["forces_eV_per_A"] - force_ref
    return {
        "abs_energy_error_eV": float(jnp.abs(energy_err)),
        "max_force_abs_error_eV_per_A": float(jnp.max(jnp.abs(force_err))),
    }

