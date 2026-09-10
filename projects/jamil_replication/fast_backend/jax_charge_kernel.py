"""JAX charge kernels for the fast eReaxFF backend."""

from __future__ import annotations

from typing import Any


def evaluate_charge_case(case: dict[str, Any]) -> dict[str, Any]:
    """Evaluate the charge model for one oracle case.

    The current validated Jamil 112C/113El target uses legacy `icharg=8` with
    explicit `MOLCHARGE` constraints. Full dynamic ACKS2 matrix solving is the
    next kernel and fails closed here.
    """

    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    charge_model = case.get("charge_model", "acks2")
    if charge_model == "fixed_molcharge":
        fixed = case.get("charge_constraints", {}).get("fixed_charges")
        if fixed is None:
            raise ValueError("fixed_molcharge case is missing fixed_charges")
        charge_out = jnp.asarray(fixed, dtype=jnp.float64)
        solver_out = jnp.zeros_like(charge_out)
        return {
            "charge_model": charge_model,
            "charge_out": charge_out,
            "solver_guess_out": solver_out,
            "status": 0,
        }
    raise NotImplementedError(f"charge model {charge_model!r} is not implemented in JAX yet")


def evaluate_charge_batch(cases: list[dict[str, Any]]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    if not cases:
        raise ValueError("charge batch is empty")
    outputs = [evaluate_charge_case(case) for case in cases]
    models = sorted({output["charge_model"] for output in outputs})
    charges = jnp.stack([output["charge_out"] for output in outputs], axis=0)
    solver = jnp.stack([output["solver_guess_out"] for output in outputs], axis=0)
    identity = jax.jit(lambda charge_array, solver_array: (charge_array, solver_array))
    charges, solver = identity(charges, solver)
    return {
        "charge_models": models,
        "charge_out": charges,
        "solver_guess_out": solver,
        "status": [int(output["status"]) for output in outputs],
    }


def charge_parity_metrics(case: dict[str, Any]) -> dict[str, float]:
    import jax.numpy as jnp

    labels = case.get("labels")
    if labels is None:
        raise ValueError("charge parity requires labels")
    predicted = evaluate_charge_case(case)
    charge_ref = jnp.asarray(labels["charge_out"], dtype=jnp.float64)
    solver_ref = jnp.asarray(labels["solver_guess_out"], dtype=jnp.float64)
    charge_err = predicted["charge_out"] - charge_ref
    solver_err = predicted["solver_guess_out"] - solver_ref
    return {
        "max_charge_abs_error": float(jnp.max(jnp.abs(charge_err))),
        "max_solver_abs_error": float(jnp.max(jnp.abs(solver_err))),
    }
