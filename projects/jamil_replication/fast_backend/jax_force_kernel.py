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
    """Parse the ffield parameters needed by the source-derived nonbonded slices."""

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
    atom_type_index_by_symbol: dict[str, int] = {}
    rat_by_symbol: dict[str, float] = {}
    rapt_by_symbol: dict[str, float] = {}
    vnq_by_symbol: dict[str, float] = {}
    rvdw_by_symbol: dict[str, float] = {}
    eps_by_symbol: dict[str, float] = {}
    alfa_by_symbol: dict[str, float] = {}
    gammavdw_by_symbol: dict[str, float] = {}
    symbols_by_index: dict[int, str] = {}
    for _ in range(n_atom_types):
        first = lines[idx].split("!", 1)[0].split()
        if len(first) < 7:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 1}")
        second = lines[idx + 1].split("!", 1)[0].split()
        if len(second) < 2:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 2}")
        symbol = first[0]
        type_index = len(symbols_by_index) + 1
        symbols_by_index[type_index] = symbol
        atom_type_index_by_symbol[symbol] = type_index
        rat_by_symbol[symbol] = float(first[1])
        rapt_by_symbol[symbol] = float(first[7])
        gamma_by_symbol[symbol] = float(first[6])
        rvdw_by_symbol[symbol] = float(first[4])
        eps_by_symbol[symbol] = float(first[5])
        alfa_by_symbol[symbol] = float(second[0])
        gammavdw_by_symbol[symbol] = float(second[1])
        third = lines[idx + 2].split("!", 1)[0].split()
        if not third:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 3}")
        vnq_by_symbol[symbol] = float(third[0])
        idx += 4
    rob1_by_pair = _default_pair_average(symbols_by_index, rat_by_symbol)
    rob2_by_pair = _default_pair_average(symbols_by_index, rapt_by_symbol)
    rob3_by_pair = _default_pair_average(symbols_by_index, vnq_by_symbol)
    bond_params_by_pair: dict[str, dict[str, float]] = {}

    n_bond_types = int(_leading_float(lines[idx]))
    idx += 2
    for _ in range(n_bond_types):
        first = lines[idx].split("!", 1)[0].split()
        second = lines[idx + 1].split("!", 1)[0].split()
        if len(first) < 10 or len(second) < 8:
            raise ValueError(f"{path}: malformed bond parameter block at line {idx + 1}")
        left_type = int(first[0])
        right_type = int(first[1])
        left_symbol = symbols_by_index[left_type]
        right_symbol = symbols_by_index[right_type]
        params = {
            "de1": float(first[2]),
            "de2": float(first[3]),
            "de3": float(first[4]),
            "psi": float(first[5]),
            "pdo": float(first[6]),
            "v13cor": float(first[7]),
            "popi": float(first[8]),
            "vover": float(first[9]),
            "psp": float(second[0]),
            "pdp": float(second[1]),
            "ptp": float(second[2]),
            "bom": float(second[3]),
            "bop1": float(second[4]),
            "bop2": float(second[5]),
            "ovc": float(second[6]),
            "pxel1": float(second[7]),
        }
        bond_params_by_pair[_pair_key(left_symbol, right_symbol)] = params
        bond_params_by_pair[_pair_key(right_symbol, left_symbol)] = params
        idx += 2

    n_offdiag = int(_leading_float(lines[idx]))
    idx += 1
    for _ in range(n_offdiag):
        fields = lines[idx].split("!", 1)[0].split()
        if len(fields) < 8:
            raise ValueError(f"{path}: malformed off-diagonal parameter line {idx + 1}")
        left_symbol = symbols_by_index[int(fields[0])]
        right_symbol = symbols_by_index[int(fields[1])]
        rodm = float(fields[3])
        rsig = float(fields[5])
        rpi = float(fields[6])
        rpi2 = float(fields[7])
        if rsig > 0.0:
            rob1_by_pair[_pair_key(left_symbol, right_symbol)] = rsig
            rob1_by_pair[_pair_key(right_symbol, left_symbol)] = rsig
        if rpi > 0.0:
            rob2_by_pair[_pair_key(left_symbol, right_symbol)] = rpi
            rob2_by_pair[_pair_key(right_symbol, left_symbol)] = rpi
        if rpi2 > 0.0:
            rob3_by_pair[_pair_key(left_symbol, right_symbol)] = rpi2
            rob3_by_pair[_pair_key(right_symbol, left_symbol)] = rpi2
        if rodm > 0.0:
            # Stored for the vdW path in the Fortran p1co table; the current
            # vdW parser keeps atom-combination values to preserve old behavior.
            pass
        idx += 1

    return {
        "bond_order_cutoff": 0.01 * general[29],
        "lower_taper_A": general[11],
        "upper_taper_A": general[12],
        "vdw_shielding_power": general[28],
        "atom_type_index_by_symbol": atom_type_index_by_symbol,
        "rat_by_symbol": rat_by_symbol,
        "rapt_by_symbol": rapt_by_symbol,
        "vnq_by_symbol": vnq_by_symbol,
        "gamma_eem_by_symbol": gamma_by_symbol,
        "rvdw_by_symbol": rvdw_by_symbol,
        "eps_by_symbol": eps_by_symbol,
        "alfa_by_symbol": alfa_by_symbol,
        "gammavdw_by_symbol": gammavdw_by_symbol,
        "rob1_by_pair": rob1_by_pair,
        "rob2_by_pair": rob2_by_pair,
        "rob3_by_pair": rob3_by_pair,
        "bond_params_by_pair": bond_params_by_pair,
    }


def _pair_key(left: str, right: str) -> str:
    return f"{left}|{right}"


def _default_pair_average(symbols_by_index: dict[int, str], values_by_symbol: dict[str, float]) -> dict[str, float]:
    pairs = {}
    for left in symbols_by_index.values():
        for right in symbols_by_index.values():
            pairs[_pair_key(left, right)] = 0.5 * (values_by_symbol[left] + values_by_symbol[right])
    return pairs


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


def reaxff_no_inner_shielded_vdw_energy(
    positions_A,
    rvdw,
    eps,
    alfa,
    gammavdw,
    vdw_shielding_power,
    taper_coefficients,
    upper_taper_A,
    cell_lengths=None,
):
    """eReaxFF shielded Morse vdW energy in eV for `nonbonnoinS`."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    rvdw = jnp.asarray(rvdw, dtype=jnp.float64)
    eps = jnp.asarray(eps, dtype=jnp.float64)
    alfa = jnp.asarray(alfa, dtype=jnp.float64)
    gammavdw = jnp.asarray(gammavdw, dtype=jnp.float64)
    displacements = positions[:, None, :] - positions[None, :, :]
    displacements = _minimum_image(displacements, cell_lengths)
    n_atoms = positions.shape[0]
    off_diagonal = ~jnp.eye(n_atoms, dtype=bool)
    r2 = jnp.sum(displacements * displacements, axis=-1)
    safe_r = jnp.sqrt(jnp.where(off_diagonal, r2, 1.0))
    p1 = jnp.sqrt(4.0 * rvdw[:, None] * rvdw[None, :])
    p2 = jnp.sqrt(eps[:, None] * eps[None, :])
    p3 = jnp.sqrt(alfa[:, None] * alfa[None, :])
    gamwco = 1.0 / (jnp.sqrt(gammavdw[:, None] * gammavdw[None, :]) ** vdw_shielding_power)
    rrw = (safe_r**vdw_shielding_power + gamwco) ** (1.0 / vdw_shielding_power)
    h1 = jnp.exp(p3 * (1.0 - rrw / p1))
    h2 = jnp.exp(0.5 * p3 * (1.0 - rrw / p1))
    taper = taper_value(safe_r, taper_coefficients)
    in_cutoff = safe_r <= upper_taper_A
    pair_energy_eV = KCAL_TO_EV * p2 * (h1 - 2.0 * h2) * taper
    return 0.5 * jnp.sum(jnp.where(off_diagonal & in_cutoff, pair_energy_eV, 0.0))


def reaxff_no_inner_shielded_nonbonded_energy(
    positions_A,
    charges,
    gamma_eem,
    rvdw,
    eps,
    alfa,
    gammavdw,
    vdw_shielding_power,
    taper_coefficients,
    upper_taper_A,
    cell_lengths=None,
):
    coulomb = reaxff_shielded_coulomb_energy(
        positions_A,
        charges,
        gamma_eem,
        taper_coefficients,
        upper_taper_A,
        cell_lengths,
    )
    vdw = reaxff_no_inner_shielded_vdw_energy(
        positions_A,
        rvdw,
        eps,
        alfa,
        gammavdw,
        vdw_shielding_power,
        taper_coefficients,
        upper_taper_A,
        cell_lengths,
    )
    return coulomb + vdw


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
    return _parameter_array(case, coulomb_parameters["gamma_eem_by_symbol"], "gammaEEM")


def _parameter_array(case: dict[str, Any], values_by_symbol: dict[str, float], name: str):
    import jax.numpy as jnp

    try:
        values = [values_by_symbol[symbol] for symbol in case["atom_types"]]
    except KeyError as exc:
        raise ValueError(f"ffield is missing {name} for atom type {exc.args[0]!r}") from exc
    return jnp.asarray(values, dtype=jnp.float64)


def _pair_parameter_matrix(case: dict[str, Any], values_by_pair: dict[str, float], name: str):
    import jax.numpy as jnp

    atom_types = case["atom_types"]
    rows = []
    for left in atom_types:
        row = []
        for right in atom_types:
            key = _pair_key(left, right)
            if key not in values_by_pair:
                raise ValueError(f"ffield is missing {name} for pair {key}")
            row.append(values_by_pair[key])
        rows.append(row)
    return jnp.asarray(rows, dtype=jnp.float64)


def _bond_parameter_matrix(case: dict[str, Any], parameters: dict[str, Any], name: str):
    import jax.numpy as jnp

    atom_types = case["atom_types"]
    bond_params = parameters["bond_params_by_pair"]
    rows = []
    mask_rows = []
    for left in atom_types:
        row = []
        mask_row = []
        for right in atom_types:
            item = bond_params.get(_pair_key(left, right))
            if item is None:
                row.append(0.0)
                mask_row.append(False)
            else:
                row.append(item[name])
                mask_row.append(True)
        rows.append(row)
        mask_rows.append(mask_row)
    return jnp.asarray(rows, dtype=jnp.float64), jnp.asarray(mask_rows, dtype=bool)


def raw_bond_orders(
    positions_A,
    rat,
    rapt,
    vnq,
    rob1,
    rob2,
    rob3,
    de1,
    pdo,
    popi,
    pdp,
    ptp,
    bop1,
    bop2,
    bond_type_mask,
    cutoff,
    cell_lengths=None,
):
    """Source-derived raw sigma/pi/pi-pi bond orders before `boncor`."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    displacements = positions[:, None, :] - positions[None, :, :]
    displacements = _minimum_image(displacements, cell_lengths)
    n_atoms = positions.shape[0]
    off_diagonal = ~jnp.eye(n_atoms, dtype=bool)
    r2 = jnp.sum(displacements * displacements, axis=-1)
    distance = jnp.sqrt(jnp.where(off_diagonal, r2, 1.0))
    sigma_enabled = (
        bond_type_mask
        & (de1 != 0.0)
        & (rat[:, None] > 0.0)
        & (rat[None, :] > 0.0)
        & (rob1 > 0.0)
    )
    pi_enabled = (
        bond_type_mask
        & (rapt[:, None] > 0.0)
        & (rapt[None, :] > 0.0)
        & (rob2 > 0.0)
    )
    pipi_enabled = (
        bond_type_mask
        & (vnq[:, None] > 0.0)
        & (vnq[None, :] > 0.0)
        & (rob3 > 0.0)
    )
    sigma_base = distance / jnp.where(sigma_enabled, rob1, 1.0)
    pi_base = distance / jnp.where(pi_enabled, rob2, 1.0)
    pipi_base = distance / jnp.where(pipi_enabled, rob3, 1.0)
    sigma = jnp.where(sigma_enabled, (1.0 + cutoff) * jnp.exp(bop1 * sigma_base**bop2), 0.0)
    pi = jnp.where(pi_enabled, jnp.exp(pdp * pi_base**ptp), 0.0)
    pipi = jnp.where(pipi_enabled, jnp.exp(pdo * pipi_base**popi), 0.0)
    total_pre_cutoff = sigma + pi + pipi
    active = off_diagonal & (total_pre_cutoff > cutoff)
    total = jnp.where(active, total_pre_cutoff - cutoff, 0.0)
    sigma_cut = jnp.where(active, sigma - cutoff, 0.0)
    pi = jnp.where(active, pi, 0.0)
    pipi = jnp.where(active, pipi, 0.0)
    return total, sigma_cut, pi, pipi, active


def reaxff_uncorrected_bond_energy(
    positions_A,
    rat,
    rapt,
    vnq,
    rob1,
    rob2,
    rob3,
    de1,
    de2,
    de3,
    psi,
    pdo,
    popi,
    pdp,
    ptp,
    psp,
    bop1,
    bop2,
    bond_type_mask,
    cutoff,
    cell_lengths=None,
):
    """Covalent bond energy using raw BO before overcoordination corrections."""

    import jax.numpy as jnp

    total, _, pi, pipi, active = raw_bond_orders(
        positions_A,
        rat,
        rapt,
        vnq,
        rob1,
        rob2,
        rob3,
        de1,
        pdo,
        popi,
        pdp,
        ptp,
        bop1,
        bop2,
        bond_type_mask,
        cutoff,
        cell_lengths,
    )
    sigma_order = jnp.maximum(total - pi - pipi, 0.0)
    sigma_power = sigma_order**psp
    pair_energy_kcal = -de1 * sigma_order * jnp.exp(psi * (1.0 - sigma_power)) - de2 * pi - de3 * pipi
    return 0.5 * jnp.sum(jnp.where(active, KCAL_TO_EV * pair_energy_kcal, 0.0))


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


def evaluate_reaxff_no_inner_shielded_nonbonded_case(case: dict[str, Any]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

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

    def energy_fn(x):
        return reaxff_no_inner_shielded_nonbonded_energy(
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
            cell_lengths,
        )

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    return {
        "kernel": "reaxff_nonbonded_no_inner_shielded",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": energy,
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def _bond_energy_inputs(case: dict[str, Any], parameters: dict[str, Any]):
    rat = _parameter_array(case, parameters["rat_by_symbol"], "covalent radius")
    rapt = _parameter_array(case, parameters["rapt_by_symbol"], "pi covalent radius")
    vnq = _parameter_array(case, parameters["vnq_by_symbol"], "double-pi covalent radius")
    rob1 = _pair_parameter_matrix(case, parameters["rob1_by_pair"], "sigma bond-order radius")
    rob2 = _pair_parameter_matrix(case, parameters["rob2_by_pair"], "pi bond-order radius")
    rob3 = _pair_parameter_matrix(case, parameters["rob3_by_pair"], "double-pi bond-order radius")
    de1, bond_mask = _bond_parameter_matrix(case, parameters, "de1")
    de2, _ = _bond_parameter_matrix(case, parameters, "de2")
    de3, _ = _bond_parameter_matrix(case, parameters, "de3")
    psi, _ = _bond_parameter_matrix(case, parameters, "psi")
    pdo, _ = _bond_parameter_matrix(case, parameters, "pdo")
    popi, _ = _bond_parameter_matrix(case, parameters, "popi")
    pdp, _ = _bond_parameter_matrix(case, parameters, "pdp")
    ptp, _ = _bond_parameter_matrix(case, parameters, "ptp")
    psp, _ = _bond_parameter_matrix(case, parameters, "psp")
    bop1, _ = _bond_parameter_matrix(case, parameters, "bop1")
    bop2, _ = _bond_parameter_matrix(case, parameters, "bop2")
    return {
        "rat": rat,
        "rapt": rapt,
        "vnq": vnq,
        "rob1": rob1,
        "rob2": rob2,
        "rob3": rob3,
        "de1": de1,
        "de2": de2,
        "de3": de3,
        "psi": psi,
        "pdo": pdo,
        "popi": popi,
        "pdp": pdp,
        "ptp": ptp,
        "psp": psp,
        "bop1": bop1,
        "bop2": bop2,
        "bond_mask": bond_mask,
        "cutoff": float(parameters["bond_order_cutoff"]),
    }


def evaluate_reaxff_uncorrected_bond_case(case: dict[str, Any]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    cell_lengths = _orthorhombic_cell_lengths(case)
    parameters = _coulomb_parameters_from_case(case)
    inputs = _bond_energy_inputs(case, parameters)

    def energy_fn(x):
        return reaxff_uncorrected_bond_energy(
            x,
            inputs["rat"],
            inputs["rapt"],
            inputs["vnq"],
            inputs["rob1"],
            inputs["rob2"],
            inputs["rob3"],
            inputs["de1"],
            inputs["de2"],
            inputs["de3"],
            inputs["psi"],
            inputs["pdo"],
            inputs["popi"],
            inputs["pdp"],
            inputs["ptp"],
            inputs["psp"],
            inputs["bop1"],
            inputs["bop2"],
            inputs["bond_mask"],
            inputs["cutoff"],
            cell_lengths,
        )

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    return {
        "kernel": "reaxff_uncorrected_bond",
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
    if kernel == "reaxff_nonbonded_no_inner_shielded":
        return evaluate_reaxff_no_inner_shielded_nonbonded_case(case)
    if kernel == "reaxff_uncorrected_bond":
        return evaluate_reaxff_uncorrected_bond_case(case)
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
    experimental_kernels = (
        "fixed_charge_coulomb",
        "reaxff_shielded_coulomb",
        "reaxff_nonbonded_no_inner_shielded",
        "reaxff_uncorrected_bond",
    )
    if kernel in experimental_kernels:
        positions = jnp.stack(
            [jnp.asarray(case["positions_A"], dtype=jnp.float64) for case in cases],
            axis=0,
        )
        charges = jnp.stack([_fixed_charges(case) for case in cases], axis=0)
        cell_lengths_by_case = [_orthorhombic_cell_lengths(case) for case in cases]
        has_cell = [cell is not None for cell in cell_lengths_by_case]
        if any(has_cell) and not all(has_cell):
            raise ValueError("mixed periodic and nonperiodic force batches are not supported")

        if kernel in ("reaxff_shielded_coulomb", "reaxff_nonbonded_no_inner_shielded"):
            parameters = [_coulomb_parameters_from_case(case) for case in cases]
            first_parameters = json.dumps(parameters[0], sort_keys=True)
            if any(json.dumps(item, sort_keys=True) != first_parameters for item in parameters):
                raise ValueError("mixed source-derived nonbonded parameter batches are not supported")
            gamma = jnp.stack([_gamma_array(case, parameters[0]) for case in cases], axis=0)
            rvdw = jnp.stack(
                [_parameter_array(case, parameters[0]["rvdw_by_symbol"], "Rvdw") for case in cases],
                axis=0,
            )
            eps = jnp.stack(
                [_parameter_array(case, parameters[0]["eps_by_symbol"], "Evdw") for case in cases],
                axis=0,
            )
            alfa = jnp.stack(
                [_parameter_array(case, parameters[0]["alfa_by_symbol"], "alfa") for case in cases],
                axis=0,
            )
            gammavdw = jnp.stack(
                [_parameter_array(case, parameters[0]["gammavdw_by_symbol"], "gammavdW") for case in cases],
                axis=0,
            )
            coefficients = _taper_coefficients(parameters[0]["lower_taper_A"], parameters[0]["upper_taper_A"])
            upper_taper = float(parameters[0]["upper_taper_A"])
            vdw_power = float(parameters[0]["vdw_shielding_power"])

            if all(has_cell):
                cells = jnp.stack(cell_lengths_by_case, axis=0)

                if kernel == "reaxff_shielded_coulomb":

                    def energy_fn(x, q, gam, _rvdw, _eps, _alfa, _gammavdw, cell):
                        return reaxff_shielded_coulomb_energy(x, q, gam, coefficients, upper_taper, cell)

                else:

                    def energy_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i, cell):
                        return reaxff_no_inner_shielded_nonbonded_energy(
                            x,
                            q,
                            gam,
                            rvdw_i,
                            eps_i,
                            alfa_i,
                            gammavdw_i,
                            vdw_power,
                            coefficients,
                            upper_taper,
                            cell,
                        )

                energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                energies, grads = jax.jit(jax.vmap(energy_and_grad))(
                    positions,
                    charges,
                    gamma,
                    rvdw,
                    eps,
                    alfa,
                    gammavdw,
                    cells,
                )
            else:

                if kernel == "reaxff_shielded_coulomb":

                    def energy_fn(x, q, gam, _rvdw, _eps, _alfa, _gammavdw):
                        return reaxff_shielded_coulomb_energy(x, q, gam, coefficients, upper_taper, None)

                else:

                    def energy_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i):
                        return reaxff_no_inner_shielded_nonbonded_energy(
                            x,
                            q,
                            gam,
                            rvdw_i,
                            eps_i,
                            alfa_i,
                            gammavdw_i,
                            vdw_power,
                            coefficients,
                            upper_taper,
                            None,
                        )

                energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                energies, grads = jax.jit(jax.vmap(energy_and_grad))(
                    positions,
                    charges,
                    gamma,
                    rvdw,
                    eps,
                    alfa,
                    gammavdw,
                )
        else:
            if kernel == "reaxff_uncorrected_bond":
                parameters = [_coulomb_parameters_from_case(case) for case in cases]
                first_parameters = json.dumps(parameters[0], sort_keys=True)
                if any(json.dumps(item, sort_keys=True) != first_parameters for item in parameters):
                    raise ValueError("mixed bond parameter batches are not supported")
                first_atom_types = cases[0]["atom_types"]
                if any(case["atom_types"] != first_atom_types for case in cases):
                    raise ValueError("mixed atom-type bond batches are not supported")
                bond_inputs = _bond_energy_inputs(cases[0], parameters[0])
                cutoff = float(bond_inputs["cutoff"])
                if all(has_cell):
                    cells = jnp.stack(cell_lengths_by_case, axis=0)

                    def energy_fn(x, cell):
                        return reaxff_uncorrected_bond_energy(
                            x,
                            bond_inputs["rat"],
                            bond_inputs["rapt"],
                            bond_inputs["vnq"],
                            bond_inputs["rob1"],
                            bond_inputs["rob2"],
                            bond_inputs["rob3"],
                            bond_inputs["de1"],
                            bond_inputs["de2"],
                            bond_inputs["de3"],
                            bond_inputs["psi"],
                            bond_inputs["pdo"],
                            bond_inputs["popi"],
                            bond_inputs["pdp"],
                            bond_inputs["ptp"],
                            bond_inputs["psp"],
                            bond_inputs["bop1"],
                            bond_inputs["bop2"],
                            bond_inputs["bond_mask"],
                            cutoff,
                            cell,
                        )

                    energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                    energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions, cells)
                else:

                    def energy_fn(x):
                        return reaxff_uncorrected_bond_energy(
                            x,
                            bond_inputs["rat"],
                            bond_inputs["rapt"],
                            bond_inputs["vnq"],
                            bond_inputs["rob1"],
                            bond_inputs["rob2"],
                            bond_inputs["rob3"],
                            bond_inputs["de1"],
                            bond_inputs["de2"],
                            bond_inputs["de3"],
                            bond_inputs["psi"],
                            bond_inputs["pdo"],
                            bond_inputs["popi"],
                            bond_inputs["pdp"],
                            bond_inputs["ptp"],
                            bond_inputs["psp"],
                            bond_inputs["bop1"],
                            bond_inputs["bop2"],
                            bond_inputs["bond_mask"],
                            cutoff,
                            None,
                        )

                    energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                    energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions)
            elif all(has_cell):
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
