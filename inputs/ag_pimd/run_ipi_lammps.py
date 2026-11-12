#!/usr/bin/env python3
"""i-PI server + LAMMPS client (fix ipi, pair_style reaxff/metal ereaxff) for an Ag/El cluster; P beads = P LAMMPS clients.
P = 1 is the classical limit: the i-PI trajectory is compared with the same system run by LAMMPS fix nve from identical positions and velocities.

usage: run_ipi_lammps.py --case <phase4_frames_ag/Ag7_1_md> --lmp BIN --plugin SO --out DIR [--beads 1] [--steps 100] [--dt-fs 0.25] [--temp 300]
Template for the server input: project/inputs/jamil_pimd/run_ipi_jax.py. Env: PYTHONPATH is set to the i-PI fork here.
Velocities (Gaussian at --temp, fixed seed) are written to the i-PI input in m/s and to the LAMMPS data file in A/fs, so both codes start identically.
"""
import argparse, json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve()
PROJ = HERE.parents[2]
IPI = PROJ.parent / "forks" / "i-pi"
LMPT = PROJ.parent / "forks" / "metal-reax-lammps" / "tests" / "lammps"
sys.path.insert(0, str(LMPT))
import run_ereaxff as R  # noqa: E402

KB = 0.0019872041          # kcal/mol/K
AMU_A2_FS2 = 2390.057      # kcal/mol per (amu A^2/fs^2)

ap = argparse.ArgumentParser()
ap.add_argument("--case", required=True); ap.add_argument("--lmp", required=True); ap.add_argument("--plugin", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--beads", type=int, default=1); ap.add_argument("--steps", type=int, default=100)
ap.add_argument("--dt-fs", type=float, default=0.25); ap.add_argument("--temp", type=float, default=300.0); ap.add_argument("--taper", type=float, default=5.7)
ap.add_argument("--zero-charges", action="store_true"); ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--classical-atoms", action="store_true", help="declare all atoms classical (pointless at P=1, exercises the option)")
a = ap.parse_args()
case = Path(os.path.abspath(a.case)); out = Path(os.path.abspath(a.out)); out.mkdir(parents=True, exist_ok=True)
sym, pos, q, _ = R.read_geo(case / "geo")
pos = json.loads((case / "frames.json").read_text())["frames"][0]["pos"]
if a.zero_charges: q = [0.0] * len(q)
masses = R.ffield_masses(str(case / "ffield")); m = np.array([masses[s.upper()] for s in sym])
rng = np.random.default_rng(a.seed)
v = rng.normal(size=(len(sym), 3)) * np.sqrt(KB * a.temp / (m[:, None] * AMU_A2_FS2))    # A/fs
v -= (m[:, None] * v).sum(0) / m.sum()                                                    # no net momentum

# --- common starting files -------------------------------------------------------------------------------------------------------------
n = len(sym)
(out / "init.xyz").write_text("\n".join([str(n), "# CELL(abcABC): 200 200 200 90 90 90 positions{angstrom} cell{angstrom}"] + ["%s %.12f %.12f %.12f" % (s, *p) for s, p in zip(sym, pos)]) + "\n")
(out / "vel.xyz").write_text("\n".join([str(n), "# CELL(abcABC): 200 200 200 90 90 90 velocity{m/s} cell{angstrom}"] + ["%s %.12e %.12e %.12e" % (s, *(1e5 * x for x in vv)) for s, vv in zip(sym, v)]) + "\n")
import shutil; shutil.copy(case / "ffield", out / "ffield")   # relative path: LAMMPS cannot take a path containing a space
R.write_data(out / "data.lmp", sym, pos, q, None, False, masses)
with open(out / "data.lmp", "a") as f:
    f.write("\nVelocities\n\n" + "\n".join("%d %.12e %.12e %.12e" % (i + 1, *v[i]) for i in range(n)) + "\n")
elems = R.elements_in_order(sym)
head = f"""plugin load {os.path.relpath(a.plugin, out) if ' ' in a.plugin else a.plugin}
units real
atom_style charge
boundary f f f
newton on
read_data data.lmp
pair_style reaxff/metal NULL ereaxff yes electron_taper {a.taper}
pair_coeff * * ffield {' '.join(elems)}
comm_modify cutoff 12.0
timestep {a.dt_fs}
"""
# reference: LAMMPS integrates
(out / "in.nve").write_text(head + f"""fix i all nve
dump d all custom {a.steps} nve.dump id x y z
dump_modify d format float %.15g sort id
thermo_style custom step temp pe ke etotal
thermo_modify format float %.12g
run {a.steps}
""")
sock = "ag%d" % (os.getpid() % 100000)
clients = "".join(f"fix ip all ipi {sock} 32123 unix\nrun 100000000\n" for _ in [0])
(out / "in.ipi").write_text(head + clients)
(out / "input.xml").write_text(f"""<simulation verbosity='low'>
  <output prefix='sim'>
    <properties stride='1' filename='out'> [ step, time{{femtosecond}}, conserved{{electronvolt}}, potential{{electronvolt}}, kinetic_md{{electronvolt}}, temperature{{kelvin}} ] </properties>
    <trajectory filename='pos' stride='{a.steps}' format='xyz'> x_centroid{{angstrom}} </trajectory>
  </output>
  <total_steps>{a.steps}</total_steps>
  <prng><seed>{a.seed}</seed></prng>
  <ffsocket name='lmp' mode='unix' pbc='false'><address>{sock}</address></ffsocket>
  <system>
    <initialize nbeads='{a.beads}'>
      <cell mode='abc' units='angstrom'> [200, 200, 200] </cell>
      <file mode='xyz' units='angstrom'> init.xyz </file>
      <masses mode='manual' units='dalton'> [{', '.join(str(x) for x in m)}] </masses>
      <velocities mode='xyz' units='m/s'> vel.xyz </velocities>
    </initialize>
    <forces><force forcefield='lmp'/></forces>
    <motion mode='dynamics'><fixcom> False </fixcom><dynamics mode='nve'><timestep units='femtosecond'> {a.dt_fs} </timestep>
      {"<classical_atoms> [" + ", ".join(map(str, range(n))) + "] </classical_atoms>" if a.classical_atoms else ""}</dynamics></motion>
    <ensemble><temperature units='kelvin'> {a.temp} </temperature></ensemble>
  </system>
</simulation>
""")
env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(Path(a.lmp).resolve().parents[1] / "lib"), PYTHONPATH=str(IPI), PYTHONUNBUFFERED="1")

def lmp(inp, log):
    return subprocess.run([a.lmp, "-in", inp, "-log", "none", "-nocite"], cwd=out, env=env, capture_output=True, text=True)

r = lmp("in.nve", None)
if r.returncode: sys.exit("LAMMPS nve failed:\n" + (r.stdout + r.stderr)[-1200:])
srv = subprocess.Popen([sys.executable, str(IPI / "bin" / "i-pi"), "input.xml"], cwd=out, env=env, stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
time.sleep(4)
clis = [subprocess.Popen([a.lmp, "-in", "in.ipi", "-log", "none", "-nocite"], cwd=out, env=env, stdout=open(out / f"client{k}.log", "w"), stderr=subprocess.STDOUT) for k in range(a.beads)]
try:
    srv.wait(timeout=1800)
finally:
    for c in clis: c.terminate()
if srv.returncode: sys.exit("i-PI failed, see ipi.log")

def last_frame(path):
    L = Path(path).read_text().splitlines(); k = n + 2
    return np.array([[float(x) for x in l.split()[1:4]] for l in L[-n:]])
def lmp_dump(path):
    L = Path(path).read_text().splitlines(); i = max(j for j, l in enumerate(L) if l.startswith("ITEM: ATOMS")) + 1
    return np.array([[float(x) for x in l.split()[1:4]] for l in L[i:i + n]])
import re
rs = (out / "RESTART").read_text()          # full-precision (9 digit) bead positions in bohr; the xyz output has only 6 significant digits
qb = np.array([float(x) for x in re.search(r"<q shape='\(%d, %d\)'>\s*\[(.*?)\]" % (a.beads, 3 * n), rs, re.S).group(1).replace("\n", " ").split(",")])
xi = qb.reshape(a.beads, n, 3).mean(axis=0) * 0.529177210903
xl = lmp_dump(out / "nve.dump")
d = np.abs(xi - xl)
print("P=%d, %d steps of %.3f fs: max |x_iPI - x_LAMMPS-nve| = %.3e A (largest displacement from start %.3f A)" % (a.beads, a.steps, a.dt_fs, d.max(), np.abs(xl - np.array(pos)).max()))
