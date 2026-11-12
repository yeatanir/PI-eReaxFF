"""i-PI custom PES for ffdirect: ONE in-process LAMMPS (python module) evaluating all beads of a step, no sockets, plus the Gaussian barrier
(same formula as ipi/pes/gaussian_barrier.py) on the quantum electron. i-PI units in (bohr, Hartree), LAMMPS real units inside.
Parameters (from <parameters>): datadir, plugin, v0 (Ha), sigma (bohr), atom (0-based index of the quantum electron), mode (run1 | run0)."""
import ctypes, os, time
import numpy as np
from ipi.pes.dummy import Dummy_driver

__DRIVER_NAME__ = "lmpdrv"
__DRIVER_CLASS__ = "LmpBatch"
BOHR = 0.529177210903
KCAL_HA = 627.5094740631


class LmpBatch(Dummy_driver):
    def __init__(self, datadir="", plugin="", v0=0.0, sigma=1.0, atom=-1, mode="run1", elems="Av El", *args, **kwargs):
        import lammps
        self.v0, self.sigma, self.atom, self.mode = float(v0), float(sigma), int(atom), mode
        self.t_lmp = 0.0; self.n = 0
        L = lammps.lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"])
        d = os.path.abspath(datadir); old = os.getcwd(); os.chdir(d)
        for c in ["plugin load plug.so", "units real", "atom_style charge", "atom_modify map array", "boundary f f f", "newton on", "read_data data.lmp",
                  "pair_style reaxff/metal NULL ereaxff yes electron_taper 5.7", "pair_coeff * * ffield %s" % elems, "comm_modify cutoff 12.0",
                  "thermo_style custom step pe", "thermo_modify format float %.16g", "thermo 1000000000", "run 0"]:
            L.command(c)
        os.chdir(old); self.L = L; self.nat = L.get_natoms()
        super().__init__(*args, **kwargs)

    def lammps_eval(self, p_ang):
        L = self.L; t = time.perf_counter()
        L.scatter_atoms("x", 1, 3, (ctypes.c_double * (3 * self.nat))(*p_ang.ravel()))
        L.command("run 0" if self.mode == "run0" else "run 1 pre no post no")
        f = np.array(L.gather_atoms("f", 1, 3)).reshape(-1, 3); e = L.get_thermo("pe")
        self.t_lmp += time.perf_counter() - t; self.n += 1
        return e, f

    def one(self, cell, pos):
        p = np.asarray(pos, float).reshape(-1, 3)
        e, f = self.lammps_eval(p * BOHR)
        e_ha = e / KCAL_HA; f_ha = f * (BOHR / KCAL_HA)                    # kcal/mol/A -> Ha/bohr
        if self.v0:
            s = p[self.atom, 2]; g = self.v0 * np.exp(-0.5 * (s / self.sigma) ** 2)
            e_ha += g; f_ha[self.atom, 2] += g * s / self.sigma ** 2
        return float(e_ha), f_ha.reshape(np.asarray(pos).shape), np.asarray(cell) * 0.0, "nada"

    def compute(self, cell, pos):
        if isinstance(cell, list): return [self.one(c, p) for c, p in zip(cell, pos)]
        return self.one(cell, pos)
