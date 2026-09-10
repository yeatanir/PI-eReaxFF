"""JAX force-kernel interfaces for the fast eReaxFF backend."""

from __future__ import annotations

import hashlib
import json
from typing import Any

COULOMB_EV_A = 14.3996454784255


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


def _fixed_charges(case: dict[str, Any]):
    import jax.numpy as jnp

    constraints = case.get("charge_constraints") or {}
    charges = constraints.get("fixed_charges")
    if charges is None:
        charges = case.get("charge_in")
    if charges is None:
        raise ValueError("fixed-charge Coulomb kernel requires charge_in or fixed_charges")
    return jnp.asarray(charges, dtype=jnp.float64)


def _orthorhombic_cell_lengths(case: dict[str, Any]):
    import jax.numpy as jnp

    crystx = case.get("crystx")
    if crystx is None:
        return None
    if len(crystx) < 6:
        raise ValueError("CRYSTX must contain a, b, c, alpha, beta, gamma")
    angles = [float(value) for value in crystx[3:6]]
    if any(abs(angle - 90.0) > 1.0e-8 for angle in angles):
        raise NotImplementedError("fixed-charge Coulomb kernel currently supports only orthorhombic cells")
    lengths = jnp.asarray(crystx[:3], dtype=jnp.float64)
    if bool(jnp.any(lengths <= 0.0)):
        raise ValueError("periodic cell lengths must be positive")
    return lengths


def _minimum_image(displacements, cell_lengths):
    import jax.numpy as jnp

    if cell_lengths is None:
        return displacements
    return displacements - cell_lengths * jnp.round(displacements / cell_lengths)


def fixed_charge_coulomb_energy(positions_A, charges, cell_lengths=None):
    """Unshielded fixed-charge pair Coulomb energy in eV.

    This is an intentionally small physics slice for validating JAX batching,
    autodiff, and timing.  It is not the full eReaxFF Coulomb term.
    """

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    charges = jnp.asarray(charges, dtype=jnp.float64)
    displacements = positions[:, None, :] - positions[None, :, :]
    displacements = _minimum_image(displacements, cell_lengths)
    n_atoms = positions.shape[0]
    off_diagonal = ~jnp.eye(n_atoms, dtype=bool)
    r2 = jnp.sum(displacements * displacements, axis=-1)
    safe_r = jnp.sqrt(jnp.where(off_diagonal, r2, 1.0))
    pair_energy = COULOMB_EV_A * charges[:, None] * charges[None, :] / safe_r
    return 0.5 * jnp.sum(jnp.where(off_diagonal, pair_energy, 0.0))


def evaluate_fixed_charge_coulomb_case(case: dict[str, Any]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    charges = _fixed_charges(case)
    cell_lengths = _orthorhombic_cell_lengths(case)

    def energy_fn(x):
        return fixed_charge_coulomb_energy(x, charges, cell_lengths)

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    return {
        "kernel": "fixed_charge_coulomb",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": energy,
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def evaluate_force_case(
    case: dict[str, Any],
    *,
    allow_replay: bool = False,
    kernel: str = "validated",
) -> dict[str, Any]:
    """Evaluate one force case.

    The validated implementation supports only strict oracle-label replay.  The
    experimental fixed-charge Coulomb slice is useful for exercising real JAX
    autodiff and batching, but it is not a full eReaxFF force field.
    """

    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    if kernel == "fixed_charge_coulomb":
        return evaluate_fixed_charge_coulomb_case(case)
    if kernel != "validated":
        raise ValueError(f"unknown force kernel: {kernel}")
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


def evaluate_force_batch(
    cases: list[dict[str, Any]],
    *,
    allow_replay: bool = False,
    kernel: str = "validated",
) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    if not cases:
        raise ValueError("force batch is empty")
    if kernel == "fixed_charge_coulomb":
        positions = jnp.stack(
            [jnp.asarray(case["positions_A"], dtype=jnp.float64) for case in cases],
            axis=0,
        )
        charges = jnp.stack([_fixed_charges(case) for case in cases], axis=0)
        cell_lengths_by_case = [_orthorhombic_cell_lengths(case) for case in cases]
        has_cell = [cell is not None for cell in cell_lengths_by_case]
        if any(has_cell) and not all(has_cell):
            raise ValueError("mixed periodic and nonperiodic force batches are not supported")

        if all(has_cell):
            cells = jnp.stack(cell_lengths_by_case, axis=0)

            def energy_fn(x, q, cell):
                return fixed_charge_coulomb_energy(x, q, cell)

            energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
            energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions, charges, cells)
        else:

            def energy_fn(x, q):
                return fixed_charge_coulomb_energy(x, q, None)

            energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
            energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions, charges)

        return {
            "kernel": ["fixed_charge_coulomb"],
            "fingerprints": [geometry_fingerprint(case) for case in cases],
            "energy_eV": energies,
            "forces_eV_per_A": -grads,
            "status": [0 for _ in cases],
        }

    outputs = [evaluate_force_case(case, allow_replay=allow_replay, kernel=kernel) for case in cases]
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


def force_parity_metrics(
    case: dict[str, Any],
    *,
    allow_replay: bool = False,
    kernel: str = "validated",
) -> dict[str, float]:
    import jax.numpy as jnp

    labels = case.get("labels")
    if labels is None:
        raise ValueError("force parity requires labels")
    predicted = evaluate_force_case(case, allow_replay=allow_replay, kernel=kernel)
    energy_ref = jnp.asarray(labels["energy_eV"], dtype=jnp.float64)
    force_ref = jnp.asarray(labels["forces_eV_per_A"], dtype=jnp.float64)
    energy_err = predicted["energy_eV"] - energy_ref
    force_err = predicted["forces_eV_per_A"] - force_ref
    return {
        "abs_energy_error_eV": float(jnp.abs(energy_err)),
        "max_force_abs_error_eV_per_A": float(jnp.max(jnp.abs(force_err))),
    }
