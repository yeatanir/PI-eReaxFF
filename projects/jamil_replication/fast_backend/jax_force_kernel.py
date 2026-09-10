"""JAX force-kernel interfaces for the fast eReaxFF backend."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

COULOMB_EV_A = 14.3996454784255
KCAL_TO_EV = 0.043364
REAXFF_COULOMB_KCAL_A = 332.0638


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


def _leading_float(text: str) -> float:
    return float(text.split("!", 1)[0].split()[0])


def parse_ffield_coulomb_parameters(path: Path) -> dict[str, Any]:
    """Parse the ffield parameters needed by the eReaxFF Coulomb slice."""

    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: too short to be a ReaxFF ffield")
    idx = 1
    n_general = int(_leading_float(lines[idx]))
    idx += 1
    general = []
    for _ in range(n_general):
        general.append(_leading_float(lines[idx]))
        idx += 1
    if len(general) < 13:
        raise ValueError(f"{path}: expected at least 13 general parameters")
    n_atom_types = int(_leading_float(lines[idx]))
    idx += 4
    gamma_by_symbol: dict[str, float] = {}
    for _ in range(n_atom_types):
        first = lines[idx].split("!", 1)[0].split()
        if len(first) < 7:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 1}")
        gamma_by_symbol[first[0]] = float(first[6])
        idx += 4
    return {
        "lower_taper_A": general[11],
        "upper_taper_A": general[12],
        "gamma_eem_by_symbol": gamma_by_symbol,
    }


def _coulomb_parameters_from_case(case: dict[str, Any]) -> dict[str, Any]:
    supplied = case.get("coulomb_parameters")
    if supplied is not None:
        return supplied
    call_dir = case.get("source_call_dir")
    if call_dir is None:
        raise ValueError("reaxff_shielded_coulomb requires coulomb_parameters or source_call_dir")
    ffield = Path(call_dir) / "ffield"
    if not ffield.exists():
        raise FileNotFoundError(f"missing ffield for Coulomb parameters: {ffield}")
    return parse_ffield_coulomb_parameters(ffield)


def _taper_coefficients(lower_A: float, upper_A: float):
    import jax.numpy as jnp

    d1 = upper_A - lower_A
    if d1 <= 0.0:
        raise ValueError("upper taper radius must be greater than lower taper radius")
    d7 = d1**7.0
    a2 = lower_A * lower_A
    a3 = a2 * lower_A
    b2 = upper_A * upper_A
    b3 = b2 * upper_A
    return jnp.asarray(
        [
            (-35.0 * a3 * b2 * b2 + 21.0 * a2 * b3 * b2 + 7.0 * lower_A * b3 * b3 + b3 * b3 * upper_A) / d7,
            140.0 * a3 * b3 / d7,
            -210.0 * (a3 * b2 + a2 * b3) / d7,
            140.0 * (a3 * upper_A + 3.0 * a2 * b2 + lower_A * b3) / d7,
            -35.0 * (a3 + 9.0 * a2 * upper_A + 9.0 * lower_A * b2 + b3) / d7,
            84.0 * (a2 + 3.0 * lower_A * upper_A + b2) / d7,
            -70.0 * (lower_A + upper_A) / d7,
            20.0 / d7,
        ],
        dtype=jnp.float64,
    )


def taper_value(r_A, coefficients):
    """Evaluate the legacy 7th-order taper polynomial."""

    swc0, swc1, swc2, swc3, swc4, swc5, swc6, swc7 = coefficients
    r2 = r_A * r_A
    r3 = r2 * r_A
    return swc7 * r3 * r3 * r_A + swc6 * r3 * r3 + swc5 * r3 * r2 + swc4 * r2 * r2 + swc3 * r3 + swc2 * r2 + swc1 * r_A + swc0


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


def reaxff_shielded_coulomb_energy(positions_A, charges, gamma_eem, taper_coefficients, upper_taper_A, cell_lengths=None):
    """eReaxFF shielded/tapered Coulomb pair energy in eV.

    Formula source: Jamil eReaxFF `poten.f:nonbonnoinS`, using only the Coulomb
    `eph` term and `reac.f:tap7th` taper polynomial.
    """

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    charges = jnp.asarray(charges, dtype=jnp.float64)
    gamma_eem = jnp.asarray(gamma_eem, dtype=jnp.float64)
    displacements = positions[:, None, :] - positions[None, :, :]
    displacements = _minimum_image(displacements, cell_lengths)
    n_atoms = positions.shape[0]
    off_diagonal = ~jnp.eye(n_atoms, dtype=bool)
    r2 = jnp.sum(displacements * displacements, axis=-1)
    safe_r = jnp.sqrt(jnp.where(off_diagonal, r2, 1.0))
    gamma_pair = 1.0 / (jnp.sqrt(gamma_eem[:, None] * gamma_eem[None, :]) ** 3)
    shielded_r = (safe_r * r2 + gamma_pair) ** (1.0 / 3.0)
    taper = taper_value(safe_r, taper_coefficients)
    in_cutoff = safe_r <= upper_taper_A
    pair_energy_kcal = REAXFF_COULOMB_KCAL_A * charges[:, None] * charges[None, :] / shielded_r
    pair_energy_eV = KCAL_TO_EV * pair_energy_kcal * taper
    return 0.5 * jnp.sum(jnp.where(off_diagonal & in_cutoff, pair_energy_eV, 0.0))


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


def _gamma_array(case: dict[str, Any], coulomb_parameters: dict[str, Any]):
    import jax.numpy as jnp

    gamma_by_symbol = coulomb_parameters["gamma_eem_by_symbol"]
    try:
        gamma = [gamma_by_symbol[symbol] for symbol in case["atom_types"]]
    except KeyError as exc:
        raise ValueError(f"ffield is missing atom type {exc.args[0]!r}") from exc
    return jnp.asarray(gamma, dtype=jnp.float64)


def evaluate_reaxff_shielded_coulomb_case(case: dict[str, Any]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    charges = _fixed_charges(case)
    cell_lengths = _orthorhombic_cell_lengths(case)
    parameters = _coulomb_parameters_from_case(case)
    gamma = _gamma_array(case, parameters)
    coefficients = _taper_coefficients(parameters["lower_taper_A"], parameters["upper_taper_A"])
    upper_taper = float(parameters["upper_taper_A"])

    def energy_fn(x):
        return reaxff_shielded_coulomb_energy(x, charges, gamma, coefficients, upper_taper, cell_lengths)

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    return {
        "kernel": "reaxff_shielded_coulomb",
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
    if kernel == "reaxff_shielded_coulomb":
        return evaluate_reaxff_shielded_coulomb_case(case)
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
    if kernel in ("fixed_charge_coulomb", "reaxff_shielded_coulomb"):
        positions = jnp.stack(
            [jnp.asarray(case["positions_A"], dtype=jnp.float64) for case in cases],
            axis=0,
        )
        charges = jnp.stack([_fixed_charges(case) for case in cases], axis=0)
        cell_lengths_by_case = [_orthorhombic_cell_lengths(case) for case in cases]
        has_cell = [cell is not None for cell in cell_lengths_by_case]
        if any(has_cell) and not all(has_cell):
            raise ValueError("mixed periodic and nonperiodic force batches are not supported")

        if kernel == "reaxff_shielded_coulomb":
            parameters = [_coulomb_parameters_from_case(case) for case in cases]
            first_parameters = json.dumps(parameters[0], sort_keys=True)
            if any(json.dumps(item, sort_keys=True) != first_parameters for item in parameters):
                raise ValueError("mixed Coulomb parameter batches are not supported")
            gamma = jnp.stack([_gamma_array(case, parameters[0]) for case in cases], axis=0)
            coefficients = _taper_coefficients(parameters[0]["lower_taper_A"], parameters[0]["upper_taper_A"])
            upper_taper = float(parameters[0]["upper_taper_A"])

            if all(has_cell):
                cells = jnp.stack(cell_lengths_by_case, axis=0)

                def energy_fn(x, q, gam, cell):
                    return reaxff_shielded_coulomb_energy(x, q, gam, coefficients, upper_taper, cell)

                energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions, charges, gamma, cells)
            else:

                def energy_fn(x, q, gam):
                    return reaxff_shielded_coulomb_energy(x, q, gam, coefficients, upper_taper, None)

                energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions, charges, gamma)
        else:
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
            "kernel": [kernel],
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
