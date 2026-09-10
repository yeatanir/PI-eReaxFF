"""JAX force-kernel interfaces for the fast eReaxFF backend."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

COULOMB_EV_A = 14.3996454784255
KCAL_TO_EV = 0.043364
REAXFF_COULOMB_KCAL_A = 332.0638


def _enable_jax_x64() -> None:
    import jax

    jax.config.update("jax_enable_x64", True)


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
    _enable_jax_x64()
    import jax.numpy as jnp

    constraints = case.get("charge_constraints") or {}
    charges = constraints.get("fixed_charges")
    if charges is None:
        charges = case.get("charge_in")
    if charges is None:
        raise ValueError("fixed-charge Coulomb kernel requires charge_in or fixed_charges")
    return jnp.asarray(charges, dtype=jnp.float64)


def _orthorhombic_cell_lengths(case: dict[str, Any]):
    _enable_jax_x64()
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
    if len(general) < 30:
        raise ValueError(f"{path}: expected at least 30 general parameters")
    n_atom_types = int(_leading_float(lines[idx]))
    idx += 4
    gamma_by_symbol: dict[str, float] = {}
    atom_type_index_by_symbol: dict[str, int] = {}
    aval_by_symbol: dict[str, float] = {}
    amas_by_symbol: dict[str, float] = {}
    rat_by_symbol: dict[str, float] = {}
    rapt_by_symbol: dict[str, float] = {}
    vnq_by_symbol: dict[str, float] = {}
    stlp_by_symbol: dict[str, float] = {}
    valf_by_symbol: dict[str, float] = {}
    valp1_by_symbol: dict[str, float] = {}
    pxel2_by_symbol: dict[str, float] = {}
    vlp1_by_symbol: dict[str, float] = {}
    vovun_by_symbol: dict[str, float] = {}
    vval1_by_symbol: dict[str, float] = {}
    betael_by_symbol: dict[str, float] = {}
    vval3_by_symbol: dict[str, float] = {}
    vval4_by_symbol: dict[str, float] = {}
    alphag_by_symbol: dict[str, float] = {}
    vincr_by_symbol: dict[str, float] = {}
    bo131_by_symbol: dict[str, float] = {}
    bo132_by_symbol: dict[str, float] = {}
    bo133_by_symbol: dict[str, float] = {}
    rvdw_by_symbol: dict[str, float] = {}
    eps_by_symbol: dict[str, float] = {}
    alfa_by_symbol: dict[str, float] = {}
    gammavdw_by_symbol: dict[str, float] = {}
    symbols_by_index: dict[int, str] = {}
    for _ in range(n_atom_types):
        first = lines[idx].split("!", 1)[0].split()
        if len(first) < 9:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 1}")
        second = lines[idx + 1].split("!", 1)[0].split()
        if len(second) < 5:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 2}")
        third = lines[idx + 2].split("!", 1)[0].split()
        if len(third) < 8:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 3}")
        fourth = lines[idx + 3].split("!", 1)[0].split()
        if len(fourth) < 5:
            raise ValueError(f"{path}: malformed atom parameter line {idx + 4}")
        symbol = first[0]
        type_index = len(symbols_by_index) + 1
        symbols_by_index[type_index] = symbol
        atom_type_index_by_symbol[symbol] = type_index
        aval_by_symbol[symbol] = float(first[2])
        amas_by_symbol[symbol] = float(first[3])
        rat_by_symbol[symbol] = float(first[1])
        rapt_by_symbol[symbol] = float(first[7])
        stlp_by_symbol[symbol] = float(first[8])
        gamma_by_symbol[symbol] = float(first[6])
        rvdw_by_symbol[symbol] = float(first[4])
        eps_by_symbol[symbol] = float(first[5])
        alfa_by_symbol[symbol] = float(second[0])
        gammavdw_by_symbol[symbol] = float(second[1])
        valf_by_symbol[symbol] = float(second[2])
        valp1_by_symbol[symbol] = float(second[3])
        pxel2_by_symbol[symbol] = float(second[4])
        vnq_by_symbol[symbol] = float(third[0])
        vlp1_by_symbol[symbol] = float(third[1])
        vincr_by_symbol[symbol] = float(third[2])
        bo131_by_symbol[symbol] = float(third[3])
        bo132_by_symbol[symbol] = float(third[4])
        bo133_by_symbol[symbol] = float(third[5])
        alphag_by_symbol[symbol] = float(third[7])
        vovun_by_symbol[symbol] = float(fourth[0])
        vval1_by_symbol[symbol] = float(fourth[1])
        betael_by_symbol[symbol] = float(fourth[2])
        vval3_by_symbol[symbol] = float(fourth[3])
        vval4_by_symbol[symbol] = float(fourth[4])
        idx += 4
    rob1_by_pair = _default_pair_average(symbols_by_index, rat_by_symbol)
    rob2_by_pair = _default_pair_average(symbols_by_index, rapt_by_symbol)
    rob3_by_pair = _default_pair_average(symbols_by_index, vnq_by_symbol)
    bond_params_by_pair: dict[str, dict[str, float]] = {}
    angle_params: list[dict[str, Any]] = []
    torsion_params: list[dict[str, Any]] = []
    hydrogen_bond_params: list[dict[str, Any]] = []

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

    if idx < len(lines):
        n_angle_types = int(_leading_float(lines[idx]))
        idx += 1
        for _ in range(n_angle_types):
            fields = lines[idx].split("!", 1)[0].split()
            if len(fields) < 10:
                raise ValueError(f"{path}: malformed valence-angle parameter line {idx + 1}")
            angle_params.append(
                {
                    "left": symbols_by_index[int(fields[0])],
                    "center": symbols_by_index[int(fields[1])],
                    "right": symbols_by_index[int(fields[2])],
                    "theta0_degrees": float(fields[3]),
                    "vka": float(fields[4]),
                    "vka3": float(fields[5]),
                    "vka8": float(fields[6]),
                    "vkac": float(fields[7]),
                    "vkap": float(fields[8]),
                    "vval2": float(fields[9]),
                }
            )
            idx += 1

    if idx < len(lines):
        n_torsion_types = int(_leading_float(lines[idx]))
        idx += 1
        for _ in range(n_torsion_types):
            fields = lines[idx].split("!", 1)[0].split()
            if len(fields) < 11:
                raise ValueError(f"{path}: malformed torsion parameter line {idx + 1}")
            torsion_params.append(
                {
                    "left": symbols_by_index.get(int(fields[0]), "*"),
                    "center_left": symbols_by_index[int(fields[1])],
                    "center_right": symbols_by_index[int(fields[2])],
                    "right": symbols_by_index.get(int(fields[3]), "*"),
                    "v1": float(fields[4]),
                    "v2": float(fields[5]),
                    "v3": float(fields[6]),
                    "v4": float(fields[7]),
                    "vconj": float(fields[8]),
                    "v2bo": float(fields[9]),
                    "v3bo": float(fields[10]),
                }
            )
            idx += 1

    if idx < len(lines):
        n_hydrogen_bond_types = int(_leading_float(lines[idx]))
        idx += 1
        for _ in range(n_hydrogen_bond_types):
            fields = lines[idx].split("!", 1)[0].split()
            if len(fields) < 7:
                raise ValueError(f"{path}: malformed hydrogen-bond parameter line {idx + 1}")
            hydrogen_bond_params.append(
                {
                    "donor": symbols_by_index[int(fields[0])],
                    "hydrogen": symbols_by_index[int(fields[1])],
                    "acceptor": symbols_by_index[int(fields[2])],
                    "rhb": float(fields[3]),
                    "dehb": float(fields[4]),
                    "vhb1": float(fields[5]),
                    "vhb2": float(fields[6]),
                }
            )
            idx += 1

    return {
        "general_parameters": general,
        "bond_order_cutoff": 0.01 * general[29],
        "overcoordination_parameter": general[0],
        "overcoordination_log_parameter": general[1],
        "lower_taper_A": general[11],
        "upper_taper_A": general[12],
        "vdw_shielding_power": general[28],
        "atom_type_index_by_symbol": atom_type_index_by_symbol,
        "aval_by_symbol": aval_by_symbol,
        "amas_by_symbol": amas_by_symbol,
        "rat_by_symbol": rat_by_symbol,
        "rapt_by_symbol": rapt_by_symbol,
        "vnq_by_symbol": vnq_by_symbol,
        "stlp_by_symbol": stlp_by_symbol,
        "valf_by_symbol": valf_by_symbol,
        "valp1_by_symbol": valp1_by_symbol,
        "pxel2_by_symbol": pxel2_by_symbol,
        "vlp1_by_symbol": vlp1_by_symbol,
        "vovun_by_symbol": vovun_by_symbol,
        "vval1_by_symbol": vval1_by_symbol,
        "betael_by_symbol": betael_by_symbol,
        "vval3_by_symbol": vval3_by_symbol,
        "vval4_by_symbol": vval4_by_symbol,
        "alphag_by_symbol": alphag_by_symbol,
        "vincr_by_symbol": vincr_by_symbol,
        "bo131_by_symbol": bo131_by_symbol,
        "bo132_by_symbol": bo132_by_symbol,
        "bo133_by_symbol": bo133_by_symbol,
        "gamma_eem_by_symbol": gamma_by_symbol,
        "rvdw_by_symbol": rvdw_by_symbol,
        "eps_by_symbol": eps_by_symbol,
        "alfa_by_symbol": alfa_by_symbol,
        "gammavdw_by_symbol": gammavdw_by_symbol,
        "rob1_by_pair": rob1_by_pair,
        "rob2_by_pair": rob2_by_pair,
        "rob3_by_pair": rob3_by_pair,
        "bond_params_by_pair": bond_params_by_pair,
        "angle_params": angle_params,
        "torsion_params": torsion_params,
        "hydrogen_bond_params": hydrogen_bond_params,
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


def _total_from_terms(terms: dict[str, Any]):
    total = None
    for value in terms.values():
        total = value if total is None else total + value
    if total is None:
        raise ValueError("energy term dictionary is empty")
    return total


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


def reaxff_electron_pair_energy(
    positions_A,
    is_electron,
    is_hole,
    core_charge,
    betael,
    alphag,
    taper_coefficients,
    upper_taper_A,
    cell_lengths=None,
):
    """Jamil eReaxFF `electron` pair term for core/El and El/El interactions."""

    return _total_from_terms(
        reaxff_electron_pair_terms(
            positions_A,
            is_electron,
            is_hole,
            core_charge,
            betael,
            alphag,
            taper_coefficients,
            upper_taper_A,
            cell_lengths,
        )
    )


def reaxff_electron_pair_terms(
    positions_A,
    is_electron,
    is_hole,
    core_charge,
    betael,
    alphag,
    taper_coefficients,
    upper_taper_A,
    cell_lengths=None,
):
    """Electron-only eReaxFF pair terms split by physical contribution."""

    import jax.numpy as jnp
    from jax.scipy.special import erf

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    # The Jamil legacy electron bucket matches central-cell distances here;
    # applying minimum image over-counts El-El repulsion for the 112C/113El cell.
    displacements = positions[:, None, :] - positions[None, :, :]
    n_atoms = positions.shape[0]
    off_diagonal = ~jnp.eye(n_atoms, dtype=bool)
    r2 = jnp.sum(displacements * displacements, axis=-1)
    safe_r = jnp.sqrt(jnp.where(off_diagonal, r2, 1.0))
    taper = taper_value(safe_r, taper_coefficients)
    in_cutoff = (safe_r <= upper_taper_A) & (safe_r >= 0.001)

    is_core = ~is_electron
    electron_core = (is_electron[:, None] & is_core[None, :]) | (is_core[:, None] & is_electron[None, :])
    core_beta = jnp.where(is_core[:, None], betael[:, None], betael[None, :])
    core_alpha = jnp.maximum(jnp.where(is_core[:, None], alphag[:, None], alphag[None, :]), 0.0)
    core_z = jnp.where(is_core[:, None], core_charge[:, None], core_charge[None, :])
    electron_core_kcal = (
        -REAXFF_COULOMB_KCAL_A
        * core_beta
        * core_z
        / safe_r
        * erf(jnp.sqrt(2.0 * core_alpha) * safe_r)
    )

    electron_electron = is_electron[:, None] & is_electron[None, :]
    electron_beta = jnp.sqrt(jnp.maximum(betael[:, None] * betael[None, :], 0.0))
    electron_alpha = jnp.sqrt(jnp.maximum(alphag[:, None] * alphag[None, :], 0.0))
    electron_electron_kcal = (
        REAXFF_COULOMB_KCAL_A
        * electron_beta
        / safe_r
        * erf(jnp.sqrt(electron_alpha) * safe_r)
    )

    active = off_diagonal & in_cutoff
    return {
        "electron_core": 0.5
        * KCAL_TO_EV
        * jnp.sum(jnp.where(active & electron_core, electron_core_kcal * taper, 0.0)),
        "electron_electron": 0.5
        * KCAL_TO_EV
        * jnp.sum(jnp.where(active & electron_electron, electron_electron_kcal * taper, 0.0)),
    }


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


def reaxff_no_inner_shielded_nonbonded_terms(
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
    return {
        "coulomb": reaxff_shielded_coulomb_energy(
            positions_A,
            charges,
            gamma_eem,
            taper_coefficients,
            upper_taper_A,
            cell_lengths,
        ),
        "vdw": reaxff_no_inner_shielded_vdw_energy(
            positions_A,
            rvdw,
            eps,
            alfa,
            gammavdw,
            vdw_shielding_power,
            taper_coefficients,
            upper_taper_A,
            cell_lengths,
        ),
    }


def reaxff_corrected_bond_nonbonded_energy(
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
    bond_inputs,
    cell_lengths=None,
):
    """Sum implemented source-derived slices: corrected bond + vdW + Coulomb."""

    return _total_from_terms(
        reaxff_corrected_bond_nonbonded_terms(
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
            bond_inputs,
            cell_lengths,
        )
    )


def reaxff_corrected_bond_nonbonded_terms(
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
    bond_inputs,
    cell_lengths=None,
):
    """Return implemented source-derived energy slices as named eV terms."""

    terms = reaxff_no_inner_shielded_nonbonded_terms(
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
        cell_lengths,
    )
    terms["bond"] = reaxff_corrected_bond_energy(
        positions_A,
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
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    return terms


def reaxff_corrected_bond_nonbonded_valence_terms(
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
    bond_inputs,
    atom_type_indices,
    angle_inputs,
    torsion_inputs=None,
    topology=None,
    electron_inputs=None,
    cell_lengths=None,
):
    terms = reaxff_corrected_bond_nonbonded_terms(
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
        bond_inputs,
        cell_lengths,
    )
    if electron_inputs is not None:
        terms.update(reaxff_lone_pair_overcoordination_terms(positions_A, bond_inputs, electron_inputs, cell_lengths))
    if topology is None:
        terms.update(reaxff_valence_angle_terms(positions_A, atom_type_indices, bond_inputs, angle_inputs, cell_lengths))
    else:
        terms.update(reaxff_sparse_valence_angle_terms(positions_A, bond_inputs, topology, cell_lengths))
    if torsion_inputs is not None and topology is None:
        terms.update(reaxff_torsion_terms(positions_A, atom_type_indices, bond_inputs, torsion_inputs, cell_lengths))
    elif topology is not None:
        terms.update(reaxff_sparse_torsion_terms(positions_A, bond_inputs, topology, cell_lengths))
    if electron_inputs is not None:
        terms.update(reaxff_electron_pair_terms(
            positions_A,
            electron_inputs["is_electron"],
            electron_inputs["is_hole"],
            electron_inputs["core_charge"],
            electron_inputs["betael"],
            electron_inputs["alphag"],
            taper_coefficients,
            upper_taper_A,
            cell_lengths,
        ))
    return terms


def reaxff_corrected_bond_nonbonded_valence_energy(
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
    bond_inputs,
    atom_type_indices,
    angle_inputs,
    torsion_inputs=None,
    topology=None,
    electron_inputs=None,
    cell_lengths=None,
):
    return _total_from_terms(
        reaxff_corrected_bond_nonbonded_valence_terms(
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
            bond_inputs,
            atom_type_indices,
            angle_inputs,
            torsion_inputs,
            topology,
            electron_inputs,
            cell_lengths,
        )
    )


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
        "energy_terms_eV": {"coulomb_unshielded": energy},
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def _gamma_array(case: dict[str, Any], coulomb_parameters: dict[str, Any]):
    return _parameter_array(case, coulomb_parameters["gamma_eem_by_symbol"], "gammaEEM")


def _atom_type_index_array(case: dict[str, Any], parameters: dict[str, Any]):
    _enable_jax_x64()
    import jax.numpy as jnp

    try:
        values = [parameters["atom_type_index_by_symbol"][symbol] - 1 for symbol in case["atom_types"]]
    except KeyError as exc:
        raise ValueError(f"ffield is missing atom type index for atom type {exc.args[0]!r}") from exc
    return jnp.asarray(values, dtype=jnp.int32)


def _parameter_array(case: dict[str, Any], values_by_symbol: dict[str, float], name: str):
    _enable_jax_x64()
    import jax.numpy as jnp

    try:
        values = [values_by_symbol[symbol] for symbol in case["atom_types"]]
    except KeyError as exc:
        raise ValueError(f"ffield is missing {name} for atom type {exc.args[0]!r}") from exc
    return jnp.asarray(values, dtype=jnp.float64)


def _pair_parameter_matrix(case: dict[str, Any], values_by_pair: dict[str, float], name: str):
    _enable_jax_x64()
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


def _angle_parameter_tensors(parameters: dict[str, Any]):
    _enable_jax_x64()
    import jax.numpy as jnp

    n_types = len(parameters["atom_type_index_by_symbol"])
    names = ("theta0_degrees", "vka", "vka3", "vka8", "vkac", "vkap", "vval2")
    arrays = {name: [[ [0.0 for _ in range(n_types)] for _ in range(n_types)] for _ in range(n_types)] for name in names}
    mask = [[[False for _ in range(n_types)] for _ in range(n_types)] for _ in range(n_types)]
    type_index = parameters["atom_type_index_by_symbol"]
    for item in parameters["angle_params"]:
        left = type_index[item["left"]] - 1
        center = type_index[item["center"]] - 1
        right = type_index[item["right"]] - 1
        if abs(item["vka"]) < 0.001:
            continue
        mask[left][center][right] = True
        for name in names:
            arrays[name][left][center][right] = item[name]
    return {name: jnp.asarray(value, dtype=jnp.float64) for name, value in arrays.items()} | {
        "mask": jnp.asarray(mask, dtype=bool)
    }


def _torsion_parameter_tensors(parameters: dict[str, Any]):
    _enable_jax_x64()
    import jax.numpy as jnp

    n_types = len(parameters["atom_type_index_by_symbol"])
    wildcard = n_types
    names = ("v1", "v2", "v3", "v4", "vconj", "v2bo", "v3bo")
    arrays = {
        name: [[[[0.0 for _ in range(n_types + 1)] for _ in range(n_types)] for _ in range(n_types)] for _ in range(n_types + 1)]
        for name in names
    }
    mask = [[[[False for _ in range(n_types + 1)] for _ in range(n_types)] for _ in range(n_types)] for _ in range(n_types + 1)]
    type_index = parameters["atom_type_index_by_symbol"]
    for item in parameters["torsion_params"]:
        left = wildcard if item["left"] == "*" else type_index[item["left"]] - 1
        center_left = type_index[item["center_left"]] - 1
        center_right = type_index[item["center_right"]] - 1
        right = wildcard if item["right"] == "*" else type_index[item["right"]] - 1
        mask[left][center_left][center_right][right] = True
        for name in names:
            arrays[name][left][center_left][center_right][right] = item[name]
    return {name: jnp.asarray(value, dtype=jnp.float64) for name, value in arrays.items()} | {
        "mask": jnp.asarray(mask, dtype=bool),
        "wildcard": jnp.asarray(wildcard, dtype=jnp.int32),
    }


def _lookup_angle_parameters(parameters: dict[str, Any], left_symbol: str, center_symbol: str, right_symbol: str):
    endpoints = tuple(sorted((left_symbol, right_symbol)))
    for item in parameters["angle_params"]:
        if tuple(sorted((item["left"], item["right"]))) == endpoints and item["center"] == center_symbol:
            if abs(item["vka"]) < 0.001:
                return None
            return item
    return None


def _lookup_torsion_parameters(parameters: dict[str, Any], left_symbol: str, center_left_symbol: str, center_right_symbol: str, right_symbol: str):
    if parameters["atom_type_index_by_symbol"][center_left_symbol] > parameters["atom_type_index_by_symbol"][center_right_symbol]:
        left_symbol, center_left_symbol, center_right_symbol, right_symbol = (
            right_symbol,
            center_right_symbol,
            center_left_symbol,
            left_symbol,
        )
    elif (
        center_left_symbol == center_right_symbol
        and parameters["atom_type_index_by_symbol"][right_symbol] < parameters["atom_type_index_by_symbol"][left_symbol]
    ):
        left_symbol, center_left_symbol, center_right_symbol, right_symbol = (
            right_symbol,
            center_right_symbol,
            center_left_symbol,
            left_symbol,
        )
    for item in parameters["torsion_params"]:
        if (
            item["left"] == left_symbol
            and item["center_left"] == center_left_symbol
            and item["center_right"] == center_right_symbol
            and item["right"] == right_symbol
        ):
            return item
    for item in parameters["torsion_params"]:
        if item["left"] == "*" and item["center_left"] == center_left_symbol and item["center_right"] == center_right_symbol and item["right"] == "*":
            return item
    return None


def _sparse_topology_from_bond_orders(case: dict[str, Any], parameters: dict[str, Any], bond_inputs: dict[str, Any], cell_lengths=None):
    _enable_jax_x64()
    import jax.numpy as jnp
    import numpy as np

    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    total, _, _, _, _ = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    total_np = np.asarray(total)
    atom_types = list(case["atom_types"])
    n_atoms = len(atom_types)
    neighbors = [[idx for idx in range(n_atoms) if idx != center and total_np[idx, center] > 0.001] for center in range(n_atoms)]

    angle_indices = []
    angle_values = {name: [] for name in ("theta0_degrees", "vka", "vka3", "vka8", "vkac", "vkap", "vval2")}
    for center in range(n_atoms):
        center_neighbors = neighbors[center]
        for left_pos, left in enumerate(center_neighbors):
            for right in center_neighbors[left_pos + 1 :]:
                if total_np[left, center] * total_np[center, right] < 0.00001:
                    continue
                item = _lookup_angle_parameters(parameters, atom_types[left], atom_types[center], atom_types[right])
                if item is None:
                    continue
                angle_indices.append((left, center, right))
                for name in angle_values:
                    angle_values[name].append(item[name])

    torsion_indices = []
    torsion_values = {name: [] for name in ("v1", "v2", "v3", "v4", "vconj", "v2bo", "v3bo")}
    for center_left in range(n_atoms):
        for center_right in neighbors[center_left]:
            if center_left >= center_right:
                continue
            for left in neighbors[center_left]:
                if left == center_right:
                    continue
                for right in neighbors[center_right]:
                    if right == center_left or right == left:
                        continue
                    if total_np[center_left, center_right] * total_np[left, center_left] * total_np[center_right, right] <= 0.001:
                        continue
                    item = _lookup_torsion_parameters(
                        parameters,
                        atom_types[left],
                        atom_types[center_left],
                        atom_types[center_right],
                        atom_types[right],
                    )
                    if item is None:
                        continue
                    torsion_indices.append((left, center_left, center_right, right))
                    for name in torsion_values:
                        torsion_values[name].append(item[name])

    topology = {
        "angle_indices": jnp.asarray(angle_indices, dtype=jnp.int32).reshape((-1, 3)),
        "torsion_indices": jnp.asarray(torsion_indices, dtype=jnp.int32).reshape((-1, 4)),
    }
    for name, values in angle_values.items():
        topology[f"angle_{name}"] = jnp.asarray(values, dtype=jnp.float64)
    for name, values in torsion_values.items():
        topology[f"torsion_{name}"] = jnp.asarray(values, dtype=jnp.float64)
    return topology


def _electron_energy_inputs(case: dict[str, Any], parameters: dict[str, Any]):
    _enable_jax_x64()
    import jax.numpy as jnp

    core_charges = {
        "H": 1.0,
        "Li": 3.0,
        "Lv": 3.0,
        "B": 5.0,
        "C": 6.0,
        "N": 7.0,
        "O": 8.0,
        "F": 9.0,
        "P": 15.0,
    }
    _assert_electron_only_case(case)
    atom_types = case["atom_types"]
    return {
        "is_electron": jnp.asarray([symbol == "El" for symbol in atom_types], dtype=bool),
        "is_hole": jnp.asarray([symbol == "Ho" for symbol in atom_types], dtype=bool),
        "core_charge": jnp.asarray([core_charges.get(symbol, 0.0) for symbol in atom_types], dtype=jnp.float64),
        "betael": _parameter_array(case, parameters["betael_by_symbol"], "electron beta"),
        "alphag": _parameter_array(case, parameters["alphag_by_symbol"], "Gaussian electron alpha"),
    }


def _electron_valence_sign_array(case: dict[str, Any]):
    _enable_jax_x64()
    import jax.numpy as jnp

    valence_down = {"C", "O", "H", "F", "N", "Os"}
    valence_up = {"Li", "Lv", "P", "B", "Av", "Al"}
    values = []
    for symbol in case["atom_types"]:
        if symbol in valence_down:
            values.append(-1.0)
        elif symbol in valence_up:
            values.append(1.0)
        else:
            values.append(0.0)
    return jnp.asarray(values, dtype=jnp.float64)


def _assert_electron_only_case(case: dict[str, Any]) -> None:
    if any(symbol == "Ho" for symbol in case["atom_types"]):
        raise ValueError("JAX eReaxFF fast backend supports electron-only inputs; Ho/hole particles are unsupported")


def _bond_parameter_matrix(case: dict[str, Any], parameters: dict[str, Any], name: str):
    _enable_jax_x64()
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
    # The current legacy oracle bond-order labels match central-cell pairs.
    # Full periodic-image covalent BO needs an explicit image-list path rather
    # than collapsing all images into one minimum-image matrix.
    del cell_lengths
    displacements = positions[:, None, :] - positions[None, :, :]
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


def corrected_bond_orders(
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
    aval,
    vval3,
    bo131,
    bo132,
    bo133,
    ovc,
    v13cor,
    overcoordination_parameter,
    overcoordination_log_parameter,
    cell_lengths=None,
):
    """Apply Jamil eReaxFF `boncor` corrections to source-derived raw BO."""

    import jax.numpy as jnp

    raw_total, raw_sigma, raw_pi, raw_pipi, raw_active = raw_bond_orders(
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
    abo = jnp.sum(raw_total, axis=1)
    ovi = abo[:, None] - aval[:, None]
    ovj = abo[None, :] - aval[None, :]

    overcoord_enabled = ovc > 0.001
    exphu1 = jnp.exp(-overcoordination_log_parameter * ovi)
    exphu2 = jnp.exp(-overcoordination_log_parameter * ovj)
    exp11 = jnp.exp(-overcoordination_parameter * ovi)
    exp21 = jnp.exp(-overcoordination_parameter * ovj)
    ovcor = -(1.0 / overcoordination_log_parameter) * jnp.log(0.5 * (exphu1 + exphu2))
    huli = aval[:, None] + exp11 + exp21
    hulj = aval[None, :] + exp11 + exp21
    corr1 = huli / (huli + ovcor)
    corr2 = hulj / (hulj + ovcor)
    corrtot = jnp.where(overcoord_enabled, 0.5 * (corr1 + corr2), 1.0)

    thirteen_enabled = v13cor > 0.001
    ovi2 = abo[:, None] - vval3[:, None]
    ovj2 = abo[None, :] - vval3[None, :]
    bo131_pair = jnp.sqrt(bo131[:, None] * bo131[None, :])
    bo132_pair = jnp.sqrt(bo132[:, None] * bo132[None, :])
    bo133_pair = jnp.sqrt(bo133[:, None] * bo133[None, :])
    cor1 = bo131_pair * raw_total * raw_total - ovi2
    cor2 = bo131_pair * raw_total * raw_total - ovj2
    bocor1 = jnp.where(thirteen_enabled, 1.0 / (1.0 + jnp.exp(-bo132_pair * cor1 + bo133_pair)), 1.0)
    bocor2 = jnp.where(thirteen_enabled, 1.0 / (1.0 + jnp.exp(-bo132_pair * cor2 + bo133_pair)), 1.0)

    corrected_total = raw_total * corrtot * bocor1 * bocor2
    corrected_pi = raw_pi * corrtot * corrtot * bocor1 * bocor2
    corrected_pipi = raw_pipi * corrtot * corrtot * bocor1 * bocor2
    corrected_total = jnp.where(corrected_total >= 1.0e-10, corrected_total, 0.0)
    corrected_pi = jnp.where(corrected_pi >= 1.0e-10, corrected_pi, 0.0)
    corrected_pipi = jnp.where(corrected_pipi >= 1.0e-10, corrected_pipi, 0.0)
    corrected_sigma = jnp.maximum(corrected_total - corrected_pi - corrected_pipi, 0.0)
    active = raw_active & ((corrected_total > 0.0) | (corrected_pi > 0.0) | (corrected_pipi > 0.0))
    return corrected_total, corrected_sigma, corrected_pi, corrected_pipi, active


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


def reaxff_corrected_bond_energy(
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
    aval,
    vval3,
    bo131,
    bo132,
    bo133,
    ovc,
    v13cor,
    overcoordination_parameter,
    overcoordination_log_parameter,
    cell_lengths=None,
):
    """Covalent bond energy using `boncor`-corrected BO values."""

    import jax.numpy as jnp

    _, sigma_order, pi, pipi, active = corrected_bond_orders(
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
        aval,
        vval3,
        bo131,
        bo132,
        bo133,
        ovc,
        v13cor,
        overcoordination_parameter,
        overcoordination_log_parameter,
        cell_lengths,
    )
    sigma_power = sigma_order**psp
    pair_energy_kcal = -de1 * sigma_order * jnp.exp(psi * (1.0 - sigma_power)) - de2 * pi - de3 * pipi
    return 0.5 * jnp.sum(jnp.where(active, KCAL_TO_EV * pair_energy_kcal, 0.0))


def electron_corrected_valencies(
    positions_A,
    aval,
    stlp,
    vincr,
    electron_valence_sign,
    is_electron,
    cell_lengths=None,
):
    """Legacy `elecpos` electron sharing for electron-only eReaxFF inputs."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    displacement = positions[:, None, :] - positions[None, :, :]
    displacement = _minimum_image(displacement, cell_lengths)
    distance2 = jnp.sum(displacement * displacement, axis=-1)
    is_core = ~is_electron
    share = jnp.exp(-vincr[:, None] * distance2)
    core_electron_share = jnp.where(is_core[:, None] & is_electron[None, :], share, 0.0)
    eshare2 = jnp.sum(core_electron_share, axis=1)
    avalc = aval + electron_valence_sign * eshare2
    stlpc = stlp + jnp.where(is_core, eshare2, 0.0)
    return avalc, stlpc, eshare2


def reaxff_lone_pair_overcoordination_terms(
    positions_A,
    bond_inputs,
    electron_inputs,
    cell_lengths=None,
):
    """Source-derived `lonpar` and `ovcor` energy terms."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    total, _, pi, pipi, active = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    avalc, stlpc, eshare2 = electron_corrected_valencies(
        positions,
        bond_inputs["aval"],
        bond_inputs["stlp"],
        bond_inputs["vincr"],
        bond_inputs["electron_valence_sign"],
        electron_inputs["is_electron"],
        cell_lengths,
    )
    abo = jnp.sum(total, axis=1)
    vund = abo - stlpc
    int_half_vund = jnp.trunc(vund / 2.0)
    vlpex = vund - 2.0 * int_half_vund
    vlp = jnp.exp(-bond_inputs["general_parameters"][15] * (2.0 + vlpex) * (2.0 + vlpex)) - int_half_vund
    voptlp = 0.5 * (stlpc - avalc)
    diffvlp = voptlp - vlp
    hulp1 = 1.0 / (1.0 + jnp.exp(-75.0 * diffvlp))
    lone_pair = bond_inputs["vlp1"] * diffvlp * hulp1

    first_row = bond_inputs["amas"] <= 31.0
    dfvl = jnp.where(first_row, 1.0, 0.0)
    vlptemp = jnp.where(first_row, vlp, voptlp)
    diffvlph = dfvl * (voptlp - vlptemp)

    safe_abo = jnp.maximum(abo, 1.0e-12)
    sumov33 = jnp.sum(jnp.where(active, bond_inputs["pxel1"] * total / safe_abo[:, None], 0.0), axis=1)
    neighbor_excess = (
        jnp.exp(-bond_inputs["pxel2"][None, :] * eshare2[None, :] * sumov33[None, :])
        * (abo[None, :] - avalc[None, :])
        - dfvl[:, None] * (voptlp[None, :] - vlptemp[None, :])
    )
    sumov = jnp.sum(jnp.where(active, (pi + pipi) * neighbor_excess, 0.0), axis=1)
    sumov2 = jnp.sum(jnp.where(active, bond_inputs["vover"] * bond_inputs["de1"] * total, 0.0), axis=1)
    sumov3 = jnp.sum(jnp.where(active, bond_inputs["pxel1"] * total / safe_abo[:, None], 0.0), axis=1)

    vho = 1.0 / (1.0 + bond_inputs["general_parameters"][32] * jnp.exp(bond_inputs["general_parameters"][31] * sumov))
    vov1 = (
        jnp.exp(-bond_inputs["pxel2"] * eshare2 * sumov3)
        * (abo - avalc)
        - diffvlph * vho
    )
    hulpo = 1.0 / (1.0 + jnp.exp(bond_inputs["vovun"] * vov1))
    hulpp = 1.0 / (vov1 + avalc + 1.0e-8)
    overcoordination = sumov2 * hulpp * hulpo * vov1
    core_mask = ~electron_inputs["is_electron"]
    return {
        "lone_pair": KCAL_TO_EV * jnp.sum(jnp.where(core_mask, lone_pair, 0.0)),
        "overcoordination": KCAL_TO_EV * jnp.sum(jnp.where(core_mask, overcoordination, 0.0)),
    }


def reaxff_valence_angle_terms(
    positions_A,
    atom_type_indices,
    bond_inputs,
    angle_inputs,
    cell_lengths=None,
):
    """Source-derived `valang` energy terms for dense masked angle triples."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    total, _, pi, pipi, _ = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    n_atoms = positions.shape[0]
    i_idx = jnp.arange(n_atoms)[:, None, None]
    j_idx = jnp.arange(n_atoms)[None, :, None]
    k_idx = jnp.arange(n_atoms)[None, None, :]
    endpoint_order = i_idx < k_idx
    distinct = (i_idx != j_idx) & (j_idx != k_idx) & (i_idx != k_idx)
    boa = total[:, :, None] - 0.001
    bob = total[None, :, :] - 0.001
    graph_mask = (boa >= 0.0) & (bob >= 0.0) & ((total[:, :, None] * total[None, :, :]) >= 0.00001)

    type_i = atom_type_indices[:, None, None]
    type_j = atom_type_indices[None, :, None]
    type_k = atom_type_indices[None, None, :]
    left_type = jnp.minimum(type_i, type_k)
    right_type = jnp.maximum(type_i, type_k)
    angle_mask = angle_inputs["mask"][left_type, type_j, right_type]
    mask = endpoint_order & distinct & graph_mask & angle_mask

    vec_a = positions[None, :, None, :] - positions[:, None, None, :]
    vec_b = positions[None, :, None, :] - positions[None, None, :, :]
    vec_a = _minimum_image(vec_a, cell_lengths)
    vec_b = _minimum_image(vec_b, cell_lengths)
    norm_a = jnp.sqrt(jnp.maximum(jnp.sum(vec_a * vec_a, axis=-1), 1.0e-24))
    norm_b = jnp.sqrt(jnp.maximum(jnp.sum(vec_b * vec_b, axis=-1), 1.0e-24))
    cosine = jnp.sum(vec_a * vec_b, axis=-1) / jnp.maximum(norm_a * norm_b, 1.0e-20)
    cosine = jnp.where(mask, jnp.clip(cosine, -1.0 + 1.0e-12, 1.0 - 1.0e-12), 0.0)
    angle = jnp.arccos(cosine)

    theta0 = angle_inputs["theta0_degrees"][left_type, type_j, right_type]
    vka = angle_inputs["vka"][left_type, type_j, right_type]
    vka3 = angle_inputs["vka3"][left_type, type_j, right_type]
    vka8 = angle_inputs["vka8"][left_type, type_j, right_type]
    vkac = angle_inputs["vkac"][left_type, type_j, right_type]
    vkap = angle_inputs["vkap"][left_type, type_j, right_type]
    vval2 = angle_inputs["vval2"][left_type, type_j, right_type]

    abo = jnp.sum(total, axis=1)
    center_abo = abo[None, :, None]
    center_valf = bond_inputs["valf"][None, :, None]
    center_vval1 = bond_inputs["vval1"][None, :, None]
    center_vval3 = bond_inputs["vval3"][None, :, None]
    center_vval4 = bond_inputs["vval4"][None, :, None]
    center_aval = bond_inputs["aval"][None, :, None]
    center_stlp = bond_inputs["stlp"][None, :, None]

    vmbo = jnp.exp(jnp.sum(-(total**8), axis=1))[None, :, None]
    sbo2_base = jnp.sum(pi + pipi, axis=1)[None, :, None]
    exbo = center_abo - center_valf
    exlp1 = center_abo - center_stlp
    exlp2 = 2.0 * jnp.trunc(exlp1 / 2.0)
    exlp = exlp1 - exlp2
    vlpnc = jnp.exp(-bond_inputs["general_parameters"][15] * (2.0 + exlp) * (2.0 + exlp)) - jnp.trunc(exlp1 / 2.0)
    vlpadj = jnp.where(exlp < 0.0, vlpnc, 0.0)
    sbo2 = sbo2_base + (1.0 - vmbo) * (-exbo - bond_inputs["general_parameters"][33] * vlpadj)
    powv = bond_inputs["general_parameters"][16]
    sbo2_pow_low = jnp.maximum(sbo2, 0.0) ** powv
    sbo2_pow_high = jnp.maximum(2.0 - sbo2, 0.0) ** powv
    sbo2h = jnp.where(sbo2 <= 0.0, 0.0, jnp.where(sbo2 <= 1.0, sbo2_pow_low, jnp.where(sbo2 < 2.0, 2.0 - sbo2_pow_high, 2.0)))
    expsbo = jnp.exp(-bond_inputs["general_parameters"][17] * (2.0 - sbo2h))
    theta0_rad = (180.0 - theta0 * (1.0 - expsbo)) * (jnp.pi / 180.0)
    thdif = theta0_rad - angle
    exphu = vka * jnp.exp(-vka3 * thdif * thdif)
    exphu2 = jnp.where(vka < 0.0, vka - exphu - vka, vka - exphu)
    safe_boa = jnp.maximum(boa, 1.0e-12)
    safe_bob = jnp.maximum(bob, 1.0e-12)
    exa = jnp.exp(-center_vval1 * safe_boa**vval2)
    exb = jnp.exp(-center_vval1 * safe_bob**vval2)
    expun = jnp.exp(-vkac * exbo)
    expun2 = jnp.exp(bond_inputs["general_parameters"][14] * exbo)
    evboadj2 = center_vval4 - (center_vval4 - 1.0) * (2.0 + expun2) / (1.0 + expun + expun2)
    valence = evboadj2 * (1.0 - exa) * (1.0 - exb) * exphu2

    exbo_penalty = center_abo - center_aval
    expov = jnp.exp(bond_inputs["general_parameters"][21] * exbo_penalty)
    expov2 = jnp.exp(-bond_inputs["general_parameters"][20] * exbo_penalty)
    ecsboadj = (2.0 + expov2) / (1.0 + expov + expov2)
    penalty = (
        vkap
        * ecsboadj
        * jnp.exp(-bond_inputs["general_parameters"][19] * (boa - 2.0) * (boa - 2.0))
        * jnp.exp(-bond_inputs["general_parameters"][19] * (bob - 2.0) * (bob - 2.0))
    )

    unda = abo[:, None, None] - boa
    ovb = center_abo - center_vval3
    undc = abo[None, None, :] - bob
    angle_conjugation = (
        vka8
        * jnp.exp(-bond_inputs["general_parameters"][30] * (boa - 1.5) * (boa - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][30] * (bob - 1.5) * (bob - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][38] * unda * unda)
        * jnp.exp(-bond_inputs["general_parameters"][38] * undc * undc)
        / (1.0 + jnp.exp(bond_inputs["general_parameters"][2] * ovb))
    )

    return {
        "valence": KCAL_TO_EV * jnp.sum(jnp.where(mask, valence, 0.0)),
        "penalty": KCAL_TO_EV * jnp.sum(jnp.where(mask, penalty, 0.0)),
        "angle_conjugation": KCAL_TO_EV * jnp.sum(jnp.where(mask, angle_conjugation, 0.0)),
    }


def reaxff_sparse_valence_angle_terms(
    positions_A,
    bond_inputs,
    topology,
    cell_lengths=None,
):
    """Source-derived `valang` energy terms for sparse topology angle triples."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    total, _, pi, pipi, _ = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    indices = topology["angle_indices"]
    zero = jnp.sum(positions) * 0.0
    if indices.shape[0] == 0:
        return {"valence": zero, "penalty": zero, "angle_conjugation": zero}

    left = indices[:, 0]
    center = indices[:, 1]
    right = indices[:, 2]
    boa = total[left, center] - 0.001
    bob = total[center, right] - 0.001
    mask = (boa >= 0.0) & (bob >= 0.0) & ((total[left, center] * total[center, right]) >= 0.00001)

    vec_a = _minimum_image(positions[center] - positions[left], cell_lengths)
    vec_b = _minimum_image(positions[center] - positions[right], cell_lengths)
    norm_a = jnp.sqrt(jnp.maximum(jnp.sum(vec_a * vec_a, axis=-1), 1.0e-24))
    norm_b = jnp.sqrt(jnp.maximum(jnp.sum(vec_b * vec_b, axis=-1), 1.0e-24))
    cosine = jnp.sum(vec_a * vec_b, axis=-1) / jnp.maximum(norm_a * norm_b, 1.0e-20)
    angle = jnp.arccos(jnp.where(mask, jnp.clip(cosine, -1.0 + 1.0e-12, 1.0 - 1.0e-12), 0.0))

    theta0 = topology["angle_theta0_degrees"]
    vka = topology["angle_vka"]
    vka3 = topology["angle_vka3"]
    vka8 = topology["angle_vka8"]
    vkac = topology["angle_vkac"]
    vkap = topology["angle_vkap"]
    vval2 = topology["angle_vval2"]

    abo = jnp.sum(total, axis=1)
    center_abo = abo[center]
    center_valf = bond_inputs["valf"][center]
    center_vval1 = bond_inputs["vval1"][center]
    center_vval3 = bond_inputs["vval3"][center]
    center_vval4 = bond_inputs["vval4"][center]
    center_aval = bond_inputs["aval"][center]
    center_stlp = bond_inputs["stlp"][center]

    vmbo = jnp.exp(jnp.sum(-(total[center] ** 8), axis=1))
    sbo2_base = jnp.sum(pi[center] + pipi[center], axis=1)
    exbo = center_abo - center_valf
    exlp1 = center_abo - center_stlp
    exlp2 = 2.0 * jnp.trunc(exlp1 / 2.0)
    exlp = exlp1 - exlp2
    vlpnc = jnp.exp(-bond_inputs["general_parameters"][15] * (2.0 + exlp) * (2.0 + exlp)) - jnp.trunc(exlp1 / 2.0)
    vlpadj = jnp.where(exlp < 0.0, vlpnc, 0.0)
    sbo2 = sbo2_base + (1.0 - vmbo) * (-exbo - bond_inputs["general_parameters"][33] * vlpadj)
    powv = bond_inputs["general_parameters"][16]
    sbo2h = jnp.where(
        sbo2 <= 0.0,
        0.0,
        jnp.where(sbo2 <= 1.0, jnp.maximum(sbo2, 0.0) ** powv, jnp.where(sbo2 < 2.0, 2.0 - jnp.maximum(2.0 - sbo2, 0.0) ** powv, 2.0)),
    )
    expsbo = jnp.exp(-bond_inputs["general_parameters"][17] * (2.0 - sbo2h))
    theta0_rad = (180.0 - theta0 * (1.0 - expsbo)) * (jnp.pi / 180.0)
    thdif = theta0_rad - angle
    exphu = vka * jnp.exp(-vka3 * thdif * thdif)
    exphu2 = jnp.where(vka < 0.0, vka - exphu - vka, vka - exphu)
    exa = jnp.exp(-center_vval1 * jnp.maximum(boa, 1.0e-12) ** vval2)
    exb = jnp.exp(-center_vval1 * jnp.maximum(bob, 1.0e-12) ** vval2)
    expun = jnp.exp(-vkac * exbo)
    expun2 = jnp.exp(bond_inputs["general_parameters"][14] * exbo)
    evboadj2 = center_vval4 - (center_vval4 - 1.0) * (2.0 + expun2) / (1.0 + expun + expun2)
    valence = evboadj2 * (1.0 - exa) * (1.0 - exb) * exphu2

    exbo_penalty = center_abo - center_aval
    expov = jnp.exp(bond_inputs["general_parameters"][21] * exbo_penalty)
    expov2 = jnp.exp(-bond_inputs["general_parameters"][20] * exbo_penalty)
    ecsboadj = (2.0 + expov2) / (1.0 + expov + expov2)
    penalty = (
        vkap
        * ecsboadj
        * jnp.exp(-bond_inputs["general_parameters"][19] * (boa - 2.0) * (boa - 2.0))
        * jnp.exp(-bond_inputs["general_parameters"][19] * (bob - 2.0) * (bob - 2.0))
    )

    unda = abo[left] - boa
    ovb = center_abo - center_vval3
    undc = abo[right] - bob
    angle_conjugation = (
        vka8
        * jnp.exp(-bond_inputs["general_parameters"][30] * (boa - 1.5) * (boa - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][30] * (bob - 1.5) * (bob - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][38] * unda * unda)
        * jnp.exp(-bond_inputs["general_parameters"][38] * undc * undc)
        / (1.0 + jnp.exp(bond_inputs["general_parameters"][2] * ovb))
    )

    return {
        "valence": KCAL_TO_EV * jnp.sum(jnp.where(mask, valence, 0.0)),
        "penalty": KCAL_TO_EV * jnp.sum(jnp.where(mask, penalty, 0.0)),
        "angle_conjugation": KCAL_TO_EV * jnp.sum(jnp.where(mask, angle_conjugation, 0.0)),
    }


def reaxff_torsion_terms(
    positions_A,
    atom_type_indices,
    bond_inputs,
    torsion_inputs,
    cell_lengths=None,
):
    """Source-derived `torang` energy terms for dense masked torsion quadruples."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    total, _, pi, _, _ = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    n_atoms = positions.shape[0]
    i_idx = jnp.arange(n_atoms)[:, None, None, None]
    j_idx = jnp.arange(n_atoms)[None, :, None, None]
    k_idx = jnp.arange(n_atoms)[None, None, :, None]
    l_idx = jnp.arange(n_atoms)[None, None, None, :]
    distinct = (i_idx != j_idx) & (i_idx != k_idx) & (i_idx != l_idx) & (j_idx != k_idx) & (j_idx != l_idx) & (k_idx != l_idx)
    central_order = j_idx < k_idx
    boa = total[:, :, None, None] - 0.001
    bob = total[None, :, :, None] - 0.001
    boc = total[None, None, :, :] - 0.001
    graph_mask = (boa >= 0.0) & (bob >= 0.0) & (boc >= 0.0) & (
        total[:, :, None, None] * total[None, :, :, None] * total[None, None, :, :] > 0.001
    )

    type_i = atom_type_indices[:, None, None, None]
    type_j = atom_type_indices[None, :, None, None]
    type_k = atom_type_indices[None, None, :, None]
    type_l = atom_type_indices[None, None, None, :]
    wildcard = torsion_inputs["wildcard"]
    exact_mask = torsion_inputs["mask"][type_i, type_j, type_k, type_l]
    both_wild_mask = torsion_inputs["mask"][wildcard, type_j, type_k, wildcard]
    torsion_mask = exact_mask | both_wild_mask
    mask = distinct & central_order & graph_mask & torsion_mask

    def choose(name):
        values = torsion_inputs[name]
        exact = values[type_i, type_j, type_k, type_l]
        both_wild = values[wildcard, type_j, type_k, wildcard]
        return jnp.where(exact_mask, exact, both_wild)

    v1 = choose("v1")
    v2 = choose("v2")
    v3 = choose("v3")
    v4 = choose("v4")
    vconj = choose("vconj")

    rij = _minimum_image(positions[:, None, None, None, :] - positions[None, :, None, None, :], cell_lengths)
    rkj = _minimum_image(positions[None, None, :, None, :] - positions[None, :, None, None, :], cell_lengths)
    rjk = -rkj
    rlk = _minimum_image(positions[None, None, None, :, :] - positions[None, None, :, None, :], cell_lengths)

    rij_norm = jnp.sqrt(jnp.maximum(jnp.sum(rij * rij, axis=-1), 1.0e-24))
    rkj_norm = jnp.sqrt(jnp.maximum(jnp.sum(rkj * rkj, axis=-1), 1.0e-24))
    rlk_norm = jnp.sqrt(jnp.maximum(jnp.sum(rlk * rlk, axis=-1), 1.0e-24))
    coshd = jnp.sum(rij * rkj, axis=-1) / jnp.maximum(rij_norm * rkj_norm, 1.0e-20)
    coshe = jnp.sum(rjk * rlk, axis=-1) / jnp.maximum(rkj_norm * rlk_norm, 1.0e-20)
    coshd = jnp.clip(coshd, -1.0 + 1.0e-12, 1.0 - 1.0e-12)
    coshe = jnp.clip(coshe, -1.0 + 1.0e-12, 1.0 - 1.0e-12)
    sinhd = jnp.sqrt(jnp.maximum(1.0 - coshd * coshd, 0.0))
    sinhe = jnp.sqrt(jnp.maximum(1.0 - coshe * coshe, 0.0))

    normal_a = jnp.cross(rij, rkj)
    normal_b = jnp.cross(rjk, rlk)
    normal_a_norm = jnp.sqrt(jnp.maximum(jnp.sum(normal_a * normal_a, axis=-1), 1.0e-24))
    normal_b_norm = jnp.sqrt(jnp.maximum(jnp.sum(normal_b * normal_b, axis=-1), 1.0e-24))
    arg = jnp.sum(normal_a * normal_b, axis=-1) / jnp.maximum(normal_a_norm * normal_b_norm, 1.0e-20)
    arg = jnp.where(mask, jnp.clip(arg, -1.0, 1.0), 0.0)
    arg2 = arg * arg

    abo = jnp.sum(total, axis=1)
    center_left_abo = abo[None, :, None, None]
    center_right_abo = abo[None, None, :, None]
    center_left_valf = bond_inputs["valf"][None, :, None, None]
    center_right_valf = bond_inputs["valf"][None, None, :, None]
    htovt = center_left_abo + center_right_abo - center_left_valf - center_right_valf
    expov = jnp.exp(bond_inputs["general_parameters"][25] * htovt)
    expov2 = jnp.exp(-bond_inputs["general_parameters"][24] * htovt)
    etboadj = (2.0 + expov2) / (1.0 + expov + expov2)
    center_pi = pi[None, :, :, None]
    bo2t = 2.0 - center_pi - etboadj
    bocor2 = jnp.exp(v4 * bo2t * bo2t)
    ethhulp = 0.5 * v1 * (1.0 + arg) + v2 * bocor2 * (1.0 - arg2) + v3 * (0.5 + 2.0 * arg2 * arg - 1.5 * arg)
    bocor4 = (
        (1.0 - jnp.exp(-bond_inputs["general_parameters"][23] * boa))
        * (1.0 - jnp.exp(-bond_inputs["general_parameters"][23] * bob))
        * (1.0 - jnp.exp(-bond_inputs["general_parameters"][23] * boc))
    )
    torsion = sinhd * sinhe * ethhulp * bocor4

    sbo = (
        jnp.exp(-bond_inputs["general_parameters"][27] * (boa - 1.5) * (boa - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][27] * (bob - 1.5) * (bob - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][27] * (boc - 1.5) * (boc - 1.5))
    )
    torsion_conjugation = vconj * (((arg2 - 1.0) * sinhd * sinhe) + 1.0) * sbo

    return {
        "torsion": KCAL_TO_EV * jnp.sum(jnp.where(mask, torsion, 0.0)),
        "torsion_conjugation": KCAL_TO_EV * jnp.sum(jnp.where(mask, torsion_conjugation, 0.0)),
    }


def reaxff_sparse_torsion_terms(
    positions_A,
    bond_inputs,
    topology,
    cell_lengths=None,
):
    """Source-derived `torang` energy terms for sparse topology torsion quadruples."""

    import jax.numpy as jnp

    positions = jnp.asarray(positions_A, dtype=jnp.float64)
    total, _, pi, _, _ = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    indices = topology["torsion_indices"]
    zero = jnp.sum(positions) * 0.0
    if indices.shape[0] == 0:
        return {"torsion": zero, "torsion_conjugation": zero}

    left = indices[:, 0]
    center_left = indices[:, 1]
    center_right = indices[:, 2]
    right = indices[:, 3]
    boa = total[left, center_left] - 0.001
    bob = total[center_left, center_right] - 0.001
    boc = total[center_right, right] - 0.001
    mask = (boa >= 0.0) & (bob >= 0.0) & (boc >= 0.0) & (
        total[left, center_left] * total[center_left, center_right] * total[center_right, right] > 0.001
    )

    rij = _minimum_image(positions[left] - positions[center_left], cell_lengths)
    rkj = _minimum_image(positions[center_right] - positions[center_left], cell_lengths)
    rjk = -rkj
    rlk = _minimum_image(positions[right] - positions[center_right], cell_lengths)
    rij_norm = jnp.sqrt(jnp.maximum(jnp.sum(rij * rij, axis=-1), 1.0e-24))
    rkj_norm = jnp.sqrt(jnp.maximum(jnp.sum(rkj * rkj, axis=-1), 1.0e-24))
    rlk_norm = jnp.sqrt(jnp.maximum(jnp.sum(rlk * rlk, axis=-1), 1.0e-24))
    coshd = jnp.sum(rij * rkj, axis=-1) / jnp.maximum(rij_norm * rkj_norm, 1.0e-20)
    coshe = jnp.sum(rjk * rlk, axis=-1) / jnp.maximum(rkj_norm * rlk_norm, 1.0e-20)
    coshd = jnp.clip(coshd, -1.0 + 1.0e-12, 1.0 - 1.0e-12)
    coshe = jnp.clip(coshe, -1.0 + 1.0e-12, 1.0 - 1.0e-12)
    sinhd = jnp.sqrt(jnp.maximum(1.0 - coshd * coshd, 0.0))
    sinhe = jnp.sqrt(jnp.maximum(1.0 - coshe * coshe, 0.0))

    normal_a = jnp.cross(rij, rkj)
    normal_b = jnp.cross(rjk, rlk)
    normal_a_norm = jnp.sqrt(jnp.maximum(jnp.sum(normal_a * normal_a, axis=-1), 1.0e-24))
    normal_b_norm = jnp.sqrt(jnp.maximum(jnp.sum(normal_b * normal_b, axis=-1), 1.0e-24))
    arg = jnp.sum(normal_a * normal_b, axis=-1) / jnp.maximum(normal_a_norm * normal_b_norm, 1.0e-20)
    arg = jnp.where(mask, jnp.clip(arg, -1.0, 1.0), 0.0)
    arg2 = arg * arg

    v1 = topology["torsion_v1"]
    v2 = topology["torsion_v2"]
    v3 = topology["torsion_v3"]
    v4 = topology["torsion_v4"]
    vconj = topology["torsion_vconj"]
    abo = jnp.sum(total, axis=1)
    htovt = abo[center_left] + abo[center_right] - bond_inputs["valf"][center_left] - bond_inputs["valf"][center_right]
    expov = jnp.exp(bond_inputs["general_parameters"][25] * htovt)
    expov2 = jnp.exp(-bond_inputs["general_parameters"][24] * htovt)
    etboadj = (2.0 + expov2) / (1.0 + expov + expov2)
    bo2t = 2.0 - pi[center_left, center_right] - etboadj
    bocor2 = jnp.exp(v4 * bo2t * bo2t)
    ethhulp = 0.5 * v1 * (1.0 + arg) + v2 * bocor2 * (1.0 - arg2) + v3 * (0.5 + 2.0 * arg2 * arg - 1.5 * arg)
    bocor4 = (
        (1.0 - jnp.exp(-bond_inputs["general_parameters"][23] * boa))
        * (1.0 - jnp.exp(-bond_inputs["general_parameters"][23] * bob))
        * (1.0 - jnp.exp(-bond_inputs["general_parameters"][23] * boc))
    )
    torsion = sinhd * sinhe * ethhulp * bocor4

    sbo = (
        jnp.exp(-bond_inputs["general_parameters"][27] * (boa - 1.5) * (boa - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][27] * (bob - 1.5) * (bob - 1.5))
        * jnp.exp(-bond_inputs["general_parameters"][27] * (boc - 1.5) * (boc - 1.5))
    )
    torsion_conjugation = vconj * (((arg2 - 1.0) * sinhd * sinhe) + 1.0) * sbo

    return {
        "torsion": KCAL_TO_EV * jnp.sum(jnp.where(mask, torsion, 0.0)),
        "torsion_conjugation": KCAL_TO_EV * jnp.sum(jnp.where(mask, torsion_conjugation, 0.0)),
    }


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
        "energy_terms_eV": {"coulomb": energy},
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
    terms = reaxff_no_inner_shielded_nonbonded_terms(
        positions,
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
    return {
        "kernel": "reaxff_nonbonded_no_inner_shielded",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": energy,
        "energy_terms_eV": terms,
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def evaluate_bond_order_state_case(case: dict[str, Any]) -> dict[str, Any]:
    _enable_jax_x64()
    import jax.numpy as jnp

    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    cell_lengths = _orthorhombic_cell_lengths(case)
    parameters = _coulomb_parameters_from_case(case)
    bond_inputs = _bond_energy_inputs(case, parameters)
    total, _, pi, pipi, active = corrected_bond_orders(
        positions,
        bond_inputs["rat"],
        bond_inputs["rapt"],
        bond_inputs["vnq"],
        bond_inputs["rob1"],
        bond_inputs["rob2"],
        bond_inputs["rob3"],
        bond_inputs["de1"],
        bond_inputs["pdo"],
        bond_inputs["popi"],
        bond_inputs["pdp"],
        bond_inputs["ptp"],
        bond_inputs["bop1"],
        bond_inputs["bop2"],
        bond_inputs["bond_mask"],
        bond_inputs["cutoff"],
        bond_inputs["aval"],
        bond_inputs["vval3"],
        bond_inputs["bo131"],
        bond_inputs["bo132"],
        bond_inputs["bo133"],
        bond_inputs["ovc"],
        bond_inputs["v13cor"],
        bond_inputs["overcoordination_parameter"],
        bond_inputs["overcoordination_log_parameter"],
        cell_lengths,
    )
    return {
        "fingerprint": geometry_fingerprint(case),
        "abo": jnp.sum(total, axis=1),
        "bo_matrix": total,
        "bopi_matrix": pi,
        "bopi2_matrix": pipi,
        "active_matrix": active,
    }


def _bond_energy_inputs(case: dict[str, Any], parameters: dict[str, Any]):
    amas = _parameter_array(case, parameters["amas_by_symbol"], "atomic mass")
    aval = _parameter_array(case, parameters["aval_by_symbol"], "valency")
    stlp = _parameter_array(case, parameters["stlp_by_symbol"], "lone-pair valency")
    valf = _parameter_array(case, parameters["valf_by_symbol"], "valence-angle valency")
    valp1 = _parameter_array(case, parameters["valp1_by_symbol"], "undercoordination parameter")
    pxel2 = _parameter_array(case, parameters["pxel2_by_symbol"], "overcoordination parameter")
    vlp1 = _parameter_array(case, parameters["vlp1_by_symbol"], "lone-pair parameter")
    vovun = _parameter_array(case, parameters["vovun_by_symbol"], "over/undercoordination exponent")
    vval1 = _parameter_array(case, parameters["vval1_by_symbol"], "valence BO damping")
    rat = _parameter_array(case, parameters["rat_by_symbol"], "covalent radius")
    rapt = _parameter_array(case, parameters["rapt_by_symbol"], "pi covalent radius")
    vnq = _parameter_array(case, parameters["vnq_by_symbol"], "double-pi covalent radius")
    vval3 = _parameter_array(case, parameters["vval3_by_symbol"], "1-3 valency")
    vval4 = _parameter_array(case, parameters["vval4_by_symbol"], "valence undercoordination damping")
    vincr = _parameter_array(case, parameters["vincr_by_symbol"], "electron valency increment")
    bo131 = _parameter_array(case, parameters["bo131_by_symbol"], "1-3 correction parameter 1")
    bo132 = _parameter_array(case, parameters["bo132_by_symbol"], "1-3 correction parameter 2")
    bo133 = _parameter_array(case, parameters["bo133_by_symbol"], "1-3 correction parameter 3")
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
    ovc, _ = _bond_parameter_matrix(case, parameters, "ovc")
    v13cor, _ = _bond_parameter_matrix(case, parameters, "v13cor")
    pxel1, _ = _bond_parameter_matrix(case, parameters, "pxel1")
    vover, _ = _bond_parameter_matrix(case, parameters, "vover")
    return {
        "amas": amas,
        "aval": aval,
        "stlp": stlp,
        "valf": valf,
        "valp1": valp1,
        "pxel2": pxel2,
        "vlp1": vlp1,
        "vovun": vovun,
        "vval1": vval1,
        "rat": rat,
        "rapt": rapt,
        "vnq": vnq,
        "vval3": vval3,
        "vval4": vval4,
        "vincr": vincr,
        "electron_valence_sign": _electron_valence_sign_array(case),
        "bo131": bo131,
        "bo132": bo132,
        "bo133": bo133,
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
        "ovc": ovc,
        "v13cor": v13cor,
        "pxel1": pxel1,
        "vover": vover,
        "bond_mask": bond_mask,
        "cutoff": float(parameters["bond_order_cutoff"]),
        "overcoordination_parameter": float(parameters["overcoordination_parameter"]),
        "overcoordination_log_parameter": float(parameters["overcoordination_log_parameter"]),
        "general_parameters": parameters["general_parameters"],
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
        "energy_terms_eV": {"bond_uncorrected": energy},
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def evaluate_reaxff_corrected_bond_case(case: dict[str, Any]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    positions = jnp.asarray(case["positions_A"], dtype=jnp.float64)
    cell_lengths = _orthorhombic_cell_lengths(case)
    parameters = _coulomb_parameters_from_case(case)
    inputs = _bond_energy_inputs(case, parameters)

    def energy_fn(x):
        return reaxff_corrected_bond_energy(
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
            inputs["aval"],
            inputs["vval3"],
            inputs["bo131"],
            inputs["bo132"],
            inputs["bo133"],
            inputs["ovc"],
            inputs["v13cor"],
            inputs["overcoordination_parameter"],
            inputs["overcoordination_log_parameter"],
            cell_lengths,
        )

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    return {
        "kernel": "reaxff_corrected_bond",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": energy,
        "energy_terms_eV": {"bond": energy},
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def evaluate_reaxff_corrected_bond_nonbonded_case(case: dict[str, Any]) -> dict[str, Any]:
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
    bond_inputs = _bond_energy_inputs(case, parameters)

    def energy_fn(x):
        return reaxff_corrected_bond_nonbonded_energy(
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
            bond_inputs,
            cell_lengths,
        )

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    terms = reaxff_corrected_bond_nonbonded_terms(
        positions,
        charges,
        gamma,
        rvdw,
        eps,
        alfa,
        gammavdw,
        vdw_power,
        coefficients,
        upper_taper,
        bond_inputs,
        cell_lengths,
    )
    return {
        "kernel": "reaxff_corrected_bond_nonbonded",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": energy,
        "energy_terms_eV": terms,
        "forces_eV_per_A": -grad,
        "status": 0,
    }


def evaluate_reaxff_corrected_bond_nonbonded_valence_case(case: dict[str, Any]) -> dict[str, Any]:
    import jax
    import jax.numpy as jnp

    jax.config.update("jax_enable_x64", True)
    _assert_electron_only_case(case)
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
    bond_inputs = _bond_energy_inputs(case, parameters)
    atom_type_indices = _atom_type_index_array(case, parameters)
    angle_inputs = _angle_parameter_tensors(parameters)
    torsion_inputs = _torsion_parameter_tensors(parameters)
    electron_inputs = _electron_energy_inputs(case, parameters)
    topology = _sparse_topology_from_bond_orders(case, parameters, bond_inputs, cell_lengths)

    def energy_fn(x):
        return reaxff_corrected_bond_nonbonded_valence_energy(
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
            bond_inputs,
            atom_type_indices,
            angle_inputs,
            torsion_inputs,
            topology,
            electron_inputs,
            cell_lengths=cell_lengths,
        )

    energy, grad = jax.value_and_grad(energy_fn)(positions)
    terms = reaxff_corrected_bond_nonbonded_valence_terms(
        positions,
        charges,
        gamma,
        rvdw,
        eps,
        alfa,
        gammavdw,
        vdw_power,
        coefficients,
        upper_taper,
        bond_inputs,
        atom_type_indices,
        angle_inputs,
        torsion_inputs,
        topology,
        electron_inputs,
        cell_lengths=cell_lengths,
    )
    return {
        "kernel": "reaxff_corrected_bond_nonbonded_valence",
        "fingerprint": geometry_fingerprint(case),
        "energy_eV": energy,
        "energy_terms_eV": terms,
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
    if kernel == "reaxff_corrected_bond":
        return evaluate_reaxff_corrected_bond_case(case)
    if kernel == "reaxff_corrected_bond_nonbonded":
        return evaluate_reaxff_corrected_bond_nonbonded_case(case)
    if kernel == "reaxff_corrected_bond_nonbonded_valence":
        return evaluate_reaxff_corrected_bond_nonbonded_valence_case(case)
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
        "energy_terms_eV": {
            key: jnp.asarray(value, dtype=jnp.float64)
            for key, value in labels.get("energy_terms_eV", {"total": labels["energy_eV"]}).items()
        },
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
        "reaxff_corrected_bond",
        "reaxff_corrected_bond_nonbonded",
        "reaxff_corrected_bond_nonbonded_valence",
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

        if kernel in (
            "reaxff_shielded_coulomb",
            "reaxff_nonbonded_no_inner_shielded",
            "reaxff_corrected_bond_nonbonded",
            "reaxff_corrected_bond_nonbonded_valence",
        ):
            parameters = [_coulomb_parameters_from_case(case) for case in cases]
            first_parameters = json.dumps(parameters[0], sort_keys=True)
            if any(json.dumps(item, sort_keys=True) != first_parameters for item in parameters):
                raise ValueError("mixed source-derived nonbonded parameter batches are not supported")
            bond_inputs = None
            atom_type_indices = None
            angle_inputs = None
            torsion_inputs = None
            topologies = None
            electron_inputs_by_case = None
            if kernel in ("reaxff_corrected_bond_nonbonded", "reaxff_corrected_bond_nonbonded_valence"):
                first_atom_types = cases[0]["atom_types"]
                if any(case["atom_types"] != first_atom_types for case in cases):
                    raise ValueError("mixed atom-type bond/nonbonded batches are not supported")
                bond_inputs = _bond_energy_inputs(cases[0], parameters[0])
                if kernel == "reaxff_corrected_bond_nonbonded_valence":
                    atom_type_indices = _atom_type_index_array(cases[0], parameters[0])
                    angle_inputs = _angle_parameter_tensors(parameters[0])
                    torsion_inputs = _torsion_parameter_tensors(parameters[0])
                    topologies = [
                        _sparse_topology_from_bond_orders(case, parameters[0], bond_inputs, cell_lengths_by_case[idx])
                        for idx, case in enumerate(cases)
                    ]
                    first_topology_shapes = {
                        key: value.shape
                        for key, value in topologies[0].items()
                    }
                    if any(
                        {key: value.shape for key, value in topology.items()} != first_topology_shapes
                        for topology in topologies
                    ):
                        raise ValueError("mixed sparse topology sizes are not supported in a compiled batch yet")
                    electron_inputs_by_case = [_electron_energy_inputs(case, parameters[0]) for case in cases]
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
            stacked_topology = None
            stacked_electron_inputs = None
            if topologies is not None:
                stacked_topology = {
                    key: jnp.stack([topology[key] for topology in topologies], axis=0)
                    for key in topologies[0]
                }
                stacked_electron_inputs = {
                    key: jnp.stack([inputs[key] for inputs in electron_inputs_by_case], axis=0)
                    for key in electron_inputs_by_case[0]
                }

            if all(has_cell):
                cells = jnp.stack(cell_lengths_by_case, axis=0)

                if kernel == "reaxff_shielded_coulomb":

                    def energy_fn(x, q, gam, _rvdw, _eps, _alfa, _gammavdw, cell):
                        return reaxff_shielded_coulomb_energy(x, q, gam, coefficients, upper_taper, cell)

                elif kernel == "reaxff_nonbonded_no_inner_shielded":

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

                elif kernel == "reaxff_corrected_bond_nonbonded":

                    def energy_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i, cell):
                        return reaxff_corrected_bond_nonbonded_energy(
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
                            bond_inputs,
                            cell,
                        )

                else:

                    def energy_fn(
                        x,
                        q,
                        gam,
                        rvdw_i,
                        eps_i,
                        alfa_i,
                        gammavdw_i,
                        cell,
                        angle_indices,
                        torsion_indices,
                        angle_theta0_degrees,
                        angle_vka,
                        angle_vka3,
                        angle_vka8,
                        angle_vkac,
                        angle_vkap,
                        angle_vval2,
                        torsion_v1,
                        torsion_v2,
                        torsion_v3,
                        torsion_v4,
                        torsion_vconj,
                        torsion_v2bo,
                        torsion_v3bo,
                        is_electron,
                        is_hole,
                        core_charge,
                        betael,
                        alphag,
                    ):
                        topology = {
                            "angle_indices": angle_indices,
                            "torsion_indices": torsion_indices,
                            "angle_theta0_degrees": angle_theta0_degrees,
                            "angle_vka": angle_vka,
                            "angle_vka3": angle_vka3,
                            "angle_vka8": angle_vka8,
                            "angle_vkac": angle_vkac,
                            "angle_vkap": angle_vkap,
                            "angle_vval2": angle_vval2,
                            "torsion_v1": torsion_v1,
                            "torsion_v2": torsion_v2,
                            "torsion_v3": torsion_v3,
                            "torsion_v4": torsion_v4,
                            "torsion_vconj": torsion_vconj,
                            "torsion_v2bo": torsion_v2bo,
                            "torsion_v3bo": torsion_v3bo,
                        }
                        electron_inputs = {
                            "is_electron": is_electron,
                            "is_hole": is_hole,
                            "core_charge": core_charge,
                            "betael": betael,
                            "alphag": alphag,
                        }
                        return reaxff_corrected_bond_nonbonded_valence_energy(
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
                            bond_inputs,
                            atom_type_indices,
                            angle_inputs,
                            torsion_inputs,
                            topology,
                            electron_inputs,
                            cell_lengths=cell,
                        )

                energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                if kernel == "reaxff_corrected_bond_nonbonded_valence":
                    energies, grads = jax.jit(jax.vmap(energy_and_grad))(
                        positions,
                        charges,
                        gamma,
                        rvdw,
                        eps,
                        alfa,
                        gammavdw,
                        cells,
                        stacked_topology["angle_indices"],
                        stacked_topology["torsion_indices"],
                        stacked_topology["angle_theta0_degrees"],
                        stacked_topology["angle_vka"],
                        stacked_topology["angle_vka3"],
                        stacked_topology["angle_vka8"],
                        stacked_topology["angle_vkac"],
                        stacked_topology["angle_vkap"],
                        stacked_topology["angle_vval2"],
                        stacked_topology["torsion_v1"],
                        stacked_topology["torsion_v2"],
                        stacked_topology["torsion_v3"],
                        stacked_topology["torsion_v4"],
                        stacked_topology["torsion_vconj"],
                        stacked_topology["torsion_v2bo"],
                        stacked_topology["torsion_v3bo"],
                        stacked_electron_inputs["is_electron"],
                        stacked_electron_inputs["is_hole"],
                        stacked_electron_inputs["core_charge"],
                        stacked_electron_inputs["betael"],
                        stacked_electron_inputs["alphag"],
                    )
                else:
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

                elif kernel == "reaxff_nonbonded_no_inner_shielded":

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

                elif kernel == "reaxff_corrected_bond_nonbonded":

                    def energy_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i):
                        return reaxff_corrected_bond_nonbonded_energy(
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
                            bond_inputs,
                            None,
                        )

                else:

                    def energy_fn(
                        x,
                        q,
                        gam,
                        rvdw_i,
                        eps_i,
                        alfa_i,
                        gammavdw_i,
                        angle_indices,
                        torsion_indices,
                        angle_theta0_degrees,
                        angle_vka,
                        angle_vka3,
                        angle_vka8,
                        angle_vkac,
                        angle_vkap,
                        angle_vval2,
                        torsion_v1,
                        torsion_v2,
                        torsion_v3,
                        torsion_v4,
                        torsion_vconj,
                        torsion_v2bo,
                        torsion_v3bo,
                        is_electron,
                        is_hole,
                        core_charge,
                        betael,
                        alphag,
                    ):
                        topology = {
                            "angle_indices": angle_indices,
                            "torsion_indices": torsion_indices,
                            "angle_theta0_degrees": angle_theta0_degrees,
                            "angle_vka": angle_vka,
                            "angle_vka3": angle_vka3,
                            "angle_vka8": angle_vka8,
                            "angle_vkac": angle_vkac,
                            "angle_vkap": angle_vkap,
                            "angle_vval2": angle_vval2,
                            "torsion_v1": torsion_v1,
                            "torsion_v2": torsion_v2,
                            "torsion_v3": torsion_v3,
                            "torsion_v4": torsion_v4,
                            "torsion_vconj": torsion_vconj,
                            "torsion_v2bo": torsion_v2bo,
                            "torsion_v3bo": torsion_v3bo,
                        }
                        electron_inputs = {
                            "is_electron": is_electron,
                            "is_hole": is_hole,
                            "core_charge": core_charge,
                            "betael": betael,
                            "alphag": alphag,
                        }
                        return reaxff_corrected_bond_nonbonded_valence_energy(
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
                            bond_inputs,
                            atom_type_indices,
                            angle_inputs,
                            torsion_inputs,
                            topology,
                            electron_inputs,
                            cell_lengths=None,
                        )

                energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                if kernel == "reaxff_corrected_bond_nonbonded_valence":
                    energies, grads = jax.jit(jax.vmap(energy_and_grad))(
                        positions,
                        charges,
                        gamma,
                        rvdw,
                        eps,
                        alfa,
                        gammavdw,
                        stacked_topology["angle_indices"],
                        stacked_topology["torsion_indices"],
                        stacked_topology["angle_theta0_degrees"],
                        stacked_topology["angle_vka"],
                        stacked_topology["angle_vka3"],
                        stacked_topology["angle_vka8"],
                        stacked_topology["angle_vkac"],
                        stacked_topology["angle_vkap"],
                        stacked_topology["angle_vval2"],
                        stacked_topology["torsion_v1"],
                        stacked_topology["torsion_v2"],
                        stacked_topology["torsion_v3"],
                        stacked_topology["torsion_v4"],
                        stacked_topology["torsion_vconj"],
                        stacked_topology["torsion_v2bo"],
                        stacked_topology["torsion_v3bo"],
                        stacked_electron_inputs["is_electron"],
                        stacked_electron_inputs["is_hole"],
                        stacked_electron_inputs["core_charge"],
                        stacked_electron_inputs["betael"],
                        stacked_electron_inputs["alphag"],
                    )
                else:
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
            if kernel in ("reaxff_uncorrected_bond", "reaxff_corrected_bond"):
                parameters = [_coulomb_parameters_from_case(case) for case in cases]
                first_parameters = json.dumps(parameters[0], sort_keys=True)
                if any(json.dumps(item, sort_keys=True) != first_parameters for item in parameters):
                    raise ValueError("mixed bond parameter batches are not supported")
                first_atom_types = cases[0]["atom_types"]
                if any(case["atom_types"] != first_atom_types for case in cases):
                    raise ValueError("mixed atom-type bond batches are not supported")
                bond_inputs = _bond_energy_inputs(cases[0], parameters[0])
                cutoff = float(bond_inputs["cutoff"])

                def bond_energy(x, cell):
                    if kernel == "reaxff_uncorrected_bond":
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
                    return reaxff_corrected_bond_energy(
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
                        bond_inputs["aval"],
                        bond_inputs["vval3"],
                        bond_inputs["bo131"],
                        bond_inputs["bo132"],
                        bond_inputs["bo133"],
                        bond_inputs["ovc"],
                        bond_inputs["v13cor"],
                        bond_inputs["overcoordination_parameter"],
                        bond_inputs["overcoordination_log_parameter"],
                        cell,
                    )

                if all(has_cell):
                    cells = jnp.stack(cell_lengths_by_case, axis=0)

                    def energy_fn(x, cell):
                        return bond_energy(x, cell)

                    energy_and_grad = jax.value_and_grad(energy_fn, argnums=0)
                    energies, grads = jax.jit(jax.vmap(energy_and_grad))(positions, cells)
                else:

                    def energy_fn(x):
                        return bond_energy(x, None)

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

        if kernel == "fixed_charge_coulomb":
            energy_terms = {"coulomb_unshielded": energies}
        elif kernel == "reaxff_shielded_coulomb":
            energy_terms = {"coulomb": energies}
        elif kernel == "reaxff_nonbonded_no_inner_shielded":
            if all(has_cell):

                def terms_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i, cell):
                    return reaxff_no_inner_shielded_nonbonded_terms(
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

                energy_terms = jax.jit(jax.vmap(terms_fn))(positions, charges, gamma, rvdw, eps, alfa, gammavdw, cells)
            else:

                def terms_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i):
                    return reaxff_no_inner_shielded_nonbonded_terms(
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

                energy_terms = jax.jit(jax.vmap(terms_fn))(positions, charges, gamma, rvdw, eps, alfa, gammavdw)
        elif kernel == "reaxff_uncorrected_bond":
            energy_terms = {"bond_uncorrected": energies}
        elif kernel == "reaxff_corrected_bond":
            energy_terms = {"bond": energies}
        elif kernel == "reaxff_corrected_bond_nonbonded":
            if all(has_cell):

                def terms_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i, cell):
                    return reaxff_corrected_bond_nonbonded_terms(
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
                        bond_inputs,
                        cell,
                    )

                energy_terms = jax.jit(jax.vmap(terms_fn))(positions, charges, gamma, rvdw, eps, alfa, gammavdw, cells)
            else:

                def terms_fn(x, q, gam, rvdw_i, eps_i, alfa_i, gammavdw_i):
                    return reaxff_corrected_bond_nonbonded_terms(
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
                        bond_inputs,
                        None,
                    )

                energy_terms = jax.jit(jax.vmap(terms_fn))(positions, charges, gamma, rvdw, eps, alfa, gammavdw)
        elif kernel == "reaxff_corrected_bond_nonbonded_valence":
            if all(has_cell):

                def terms_fn(
                    x,
                    q,
                    gam,
                    rvdw_i,
                    eps_i,
                    alfa_i,
                    gammavdw_i,
                    cell,
                    angle_indices,
                    torsion_indices,
                    angle_theta0_degrees,
                    angle_vka,
                    angle_vka3,
                    angle_vka8,
                    angle_vkac,
                    angle_vkap,
                    angle_vval2,
                    torsion_v1,
                    torsion_v2,
                    torsion_v3,
                    torsion_v4,
                    torsion_vconj,
                    torsion_v2bo,
                    torsion_v3bo,
                    is_electron,
                    is_hole,
                    core_charge,
                    betael,
                    alphag,
                ):
                    topology = {
                        "angle_indices": angle_indices,
                        "torsion_indices": torsion_indices,
                        "angle_theta0_degrees": angle_theta0_degrees,
                        "angle_vka": angle_vka,
                        "angle_vka3": angle_vka3,
                        "angle_vka8": angle_vka8,
                        "angle_vkac": angle_vkac,
                        "angle_vkap": angle_vkap,
                        "angle_vval2": angle_vval2,
                        "torsion_v1": torsion_v1,
                        "torsion_v2": torsion_v2,
                        "torsion_v3": torsion_v3,
                        "torsion_v4": torsion_v4,
                        "torsion_vconj": torsion_vconj,
                        "torsion_v2bo": torsion_v2bo,
                        "torsion_v3bo": torsion_v3bo,
                    }
                    electron_inputs = {
                        "is_electron": is_electron,
                        "is_hole": is_hole,
                        "core_charge": core_charge,
                        "betael": betael,
                        "alphag": alphag,
                    }
                    return reaxff_corrected_bond_nonbonded_valence_terms(
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
                        bond_inputs,
                        atom_type_indices,
                        angle_inputs,
                        torsion_inputs,
                        topology,
                        electron_inputs,
                        cell_lengths=cell,
                    )

                energy_terms = jax.jit(jax.vmap(terms_fn))(
                    positions,
                    charges,
                    gamma,
                    rvdw,
                    eps,
                    alfa,
                    gammavdw,
                    cells,
                    stacked_topology["angle_indices"],
                    stacked_topology["torsion_indices"],
                    stacked_topology["angle_theta0_degrees"],
                    stacked_topology["angle_vka"],
                    stacked_topology["angle_vka3"],
                    stacked_topology["angle_vka8"],
                    stacked_topology["angle_vkac"],
                    stacked_topology["angle_vkap"],
                    stacked_topology["angle_vval2"],
                    stacked_topology["torsion_v1"],
                    stacked_topology["torsion_v2"],
                    stacked_topology["torsion_v3"],
                    stacked_topology["torsion_v4"],
                    stacked_topology["torsion_vconj"],
                    stacked_topology["torsion_v2bo"],
                    stacked_topology["torsion_v3bo"],
                    stacked_electron_inputs["is_electron"],
                    stacked_electron_inputs["is_hole"],
                    stacked_electron_inputs["core_charge"],
                    stacked_electron_inputs["betael"],
                    stacked_electron_inputs["alphag"],
                )
            else:

                def terms_fn(
                    x,
                    q,
                    gam,
                    rvdw_i,
                    eps_i,
                    alfa_i,
                    gammavdw_i,
                    angle_indices,
                    torsion_indices,
                    angle_theta0_degrees,
                    angle_vka,
                    angle_vka3,
                    angle_vka8,
                    angle_vkac,
                    angle_vkap,
                    angle_vval2,
                    torsion_v1,
                    torsion_v2,
                    torsion_v3,
                    torsion_v4,
                    torsion_vconj,
                    torsion_v2bo,
                    torsion_v3bo,
                    is_electron,
                    is_hole,
                    core_charge,
                    betael,
                    alphag,
                ):
                    topology = {
                        "angle_indices": angle_indices,
                        "torsion_indices": torsion_indices,
                        "angle_theta0_degrees": angle_theta0_degrees,
                        "angle_vka": angle_vka,
                        "angle_vka3": angle_vka3,
                        "angle_vka8": angle_vka8,
                        "angle_vkac": angle_vkac,
                        "angle_vkap": angle_vkap,
                        "angle_vval2": angle_vval2,
                        "torsion_v1": torsion_v1,
                        "torsion_v2": torsion_v2,
                        "torsion_v3": torsion_v3,
                        "torsion_v4": torsion_v4,
                        "torsion_vconj": torsion_vconj,
                        "torsion_v2bo": torsion_v2bo,
                        "torsion_v3bo": torsion_v3bo,
                    }
                    electron_inputs = {
                        "is_electron": is_electron,
                        "is_hole": is_hole,
                        "core_charge": core_charge,
                        "betael": betael,
                        "alphag": alphag,
                    }
                    return reaxff_corrected_bond_nonbonded_valence_terms(
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
                        bond_inputs,
                        atom_type_indices,
                        angle_inputs,
                        torsion_inputs,
                        topology,
                        electron_inputs,
                        cell_lengths=None,
                    )

                energy_terms = jax.jit(jax.vmap(terms_fn))(
                    positions,
                    charges,
                    gamma,
                    rvdw,
                    eps,
                    alfa,
                    gammavdw,
                    stacked_topology["angle_indices"],
                    stacked_topology["torsion_indices"],
                    stacked_topology["angle_theta0_degrees"],
                    stacked_topology["angle_vka"],
                    stacked_topology["angle_vka3"],
                    stacked_topology["angle_vka8"],
                    stacked_topology["angle_vkac"],
                    stacked_topology["angle_vkap"],
                    stacked_topology["angle_vval2"],
                    stacked_topology["torsion_v1"],
                    stacked_topology["torsion_v2"],
                    stacked_topology["torsion_v3"],
                    stacked_topology["torsion_v4"],
                    stacked_topology["torsion_vconj"],
                    stacked_topology["torsion_v2bo"],
                    stacked_topology["torsion_v3bo"],
                    stacked_electron_inputs["is_electron"],
                    stacked_electron_inputs["is_hole"],
                    stacked_electron_inputs["core_charge"],
                    stacked_electron_inputs["betael"],
                    stacked_electron_inputs["alphag"],
                )
        else:
            energy_terms = {"total": energies}

        return {
            "kernel": [kernel],
            "fingerprints": [geometry_fingerprint(case) for case in cases],
            "energy_eV": energies,
            "energy_terms_eV": energy_terms,
            "forces_eV_per_A": -grads,
            "status": [0 for _ in cases],
        }

    outputs = [evaluate_force_case(case, allow_replay=allow_replay, kernel=kernel) for case in cases]
    energies = jnp.stack([output["energy_eV"] for output in outputs], axis=0)
    forces = jnp.stack([output["forces_eV_per_A"] for output in outputs], axis=0)
    identity = jax.jit(lambda energy_array, force_array: (energy_array, force_array))
    energies, forces = identity(energies, forces)
    term_keys = sorted({key for output in outputs for key in output.get("energy_terms_eV", {})})
    energy_terms = {
        key: jnp.stack(
            [
                output.get("energy_terms_eV", {}).get(key, jnp.asarray(0.0, dtype=jnp.float64))
                for output in outputs
            ],
            axis=0,
        )
        for key in term_keys
    }
    return {
        "kernel": sorted({output["kernel"] for output in outputs}),
        "fingerprints": [str(output["fingerprint"]) for output in outputs],
        "energy_eV": energies,
        "energy_terms_eV": energy_terms,
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
