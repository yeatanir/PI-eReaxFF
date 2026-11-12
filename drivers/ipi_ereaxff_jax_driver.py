"""i-PI Python-driver PES: eReaxFF through the JAX kernel of the old prototype (interim engine, plan phase 2).

Use with:  EREAXFF_FFIELD=<ffield> EREAXFF_GEO=<geo> \
           <venv-with-jax>/bin/python <i-pi>/bin/i-pi-py_driver -m custom -P this_file.py -u -a <address>

Physics: the kernel is the central-cell (non-periodic) fixed-charge `icharg 8` eReaxFF of the Jamil target;
the rebuilt `ereac` (non-periodic, electron taper 10 A) reproduces it to machine precision on the V2 reference geometry,
so this driver is *flake* physics (no periodic images), good for integrator/stability work, not for science claims.

Atom types and fixed charges (MOLCHARGE) come from the geo file; positions come from i-PI.
Units: i-PI works in atomic units (bohr, Hartree); the kernel in A and eV.
"""
import os
import sys
from pathlib import Path

import numpy as np

from ipi.pes.dummy import Dummy_driver

__DRIVER_NAME__ = "ereaxff_jax"
__DRIVER_CLASS__ = "EReaxFF_JAX_driver"

BOHR_A = 0.529177210903
HARTREE_EV = 27.211386245988

# read-only reference to the old prototype's kernel (not copied, not modified)
KERNEL_DIR = Path(os.environ.get(
    "EREAXFF_JAX_KERNEL_DIR",
    "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_pimd_modular/projects/jamil_replication/fast_backend"))


def read_geo(path):
    """Return atom symbols, positions (A), fixed charges from MOLCHARGE, crystx (or None)."""
    symbols, pos, groups, crystx = [], [], [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith("HETATM"):
            w = line.split()
            symbols.append(w[2])
            pos.append([float(w[4]), float(w[5]), float(w[6])])
        elif line.startswith("MOLCHARGE"):
            w = line.split()
            groups.append((int(w[1]), int(w[2]), float(w[3])))
        elif line.startswith("CRYSTX"):
            crystx = [float(x) for x in line.split()[1:7]]
    q = np.zeros(len(symbols))
    for a, b, c in groups:
        q[a - 1:b] = c / (b - a + 1)
    return symbols, np.array(pos), q, crystx


class EReaxFF_JAX_driver(Dummy_driver):
    def __init__(self, *args, **kwargs):
        ffield = os.environ["EREAXFF_FFIELD"]
        geo = os.environ["EREAXFF_GEO"]
        self.kernel = os.environ.get("EREAXFF_KERNEL", "reaxff_corrected_bond_nonbonded_valence")
        self.use_crystx = os.environ.get("EREAXFF_PERIODIC", "0") == "1"
        sys.path.insert(0, str(KERNEL_DIR))
        import jax
        jax.config.update("jax_enable_x64", True)
        from jax_force_kernel import evaluate_force_batch, parse_ffield_coulomb_parameters
        self._evaluate = evaluate_force_batch
        self.params = parse_ffield_coulomb_parameters(Path(ffield))
        self.symbols, _, self.q, crystx = read_geo(geo)
        self.crystx = crystx if self.use_crystx else None
        self.natoms = len(self.symbols)
        super().__init__(*args, **kwargs)

    def _cases(self, positions_bohr):
        cases = []
        for p in positions_bohr:
            case = {
                "atom_types": self.symbols,
                "positions_A": (np.asarray(p).reshape(-1, 3) * BOHR_A).tolist(),
                "charge_constraints": {"fixed_charges": self.q.tolist()},
                "charge_in": self.q.tolist(),
                "solver_guess_in": [0.0] * self.natoms,
                "coulomb_parameters": self.params,
            }
            if self.crystx is not None:
                case["crystx"] = self.crystx
            cases.append(case)
        return cases

    def evaluate_positions(self, positions_bohr):
        out = self._evaluate(self._cases(positions_bohr), kernel=self.kernel)
        e = np.asarray(out["energy_eV"], float) / HARTREE_EV
        f = np.asarray(out["forces_eV_per_A"], float) * (BOHR_A / HARTREE_EV)
        if not (np.all(np.isfinite(e)) and np.all(np.isfinite(f))):
            raise ValueError("non-finite energy or force from the JAX kernel")
        return e, f

    def compute_structure(self, cell, pos):
        e, f = self.evaluate_positions([pos])
        return float(e[0]), f[0].reshape(pos.shape), cell * 0.0, "nada"

    def compute(self, cell, pos):
        if not isinstance(cell, list):
            return self.compute_structure(cell, pos)
        e, f = self.evaluate_positions(pos)
        return [(float(e[i]), f[i].reshape(np.asarray(pos[i]).shape), np.asarray(cell[i]) * 0.0, "nada")
                for i in range(len(pos))]
