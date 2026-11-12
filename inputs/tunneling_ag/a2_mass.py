#!/usr/bin/env python3
"""A2: the El dynamical mass (i-PI) and the LAMMPS data-file El mass must not change the engine's energy/forces.
Part 1 (offline, LAMMPS plugin): same geometries, data-file El mass 1.0 / 0.5 / 0.1 / 25.0 (25 is above the legacy mass>21 threshold: a stress test).
Part 2 (i-PI): three short P=2 runs with i-PI El mass 1.0 / 0.5 / 0.1: the step-0 client forces (identical positions) must be identical, and the mass must have reached i-PI.
usage: a2_mass.py"""
import re, subprocess, sys, shutil, os
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve(); PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ.parent / "forks/metal-reax-lammps/tests/lammps")); sys.path.insert(0, str(HERE.parent))
import run_ereaxff as R
import importlib.util
spec = importlib.util.spec_from_file_location("a1", HERE.parent / "a1_compare.py")
src = (HERE.parent / "a1_compare.py").read_text().rsplit("main()", 1)[0]
ns = {"__file__": str(HERE.parent / "a1_compare.py"), "__name__": "a1"}; exec(compile(src, "a1_compare.py", "exec"), ns)
sym, _, q, _ = R.read_geo(ns["CASE"] / "geo"); masses = R.ffield_masses(str(ns["CASE"] / "ffield"))
dump = ns["read_dump"](PROJ.parent / "scratch/tun_A1/p8/client0/fd.dump")
geos = [dump[j][1][:, :3] for j in (3, 300, 600)]
print("Part 1: offline LAMMPS, data-file El mass varied (ffield mass untouched)")
ref = None
for m in (1.0, 0.5, 0.1, 25.0):
    mm = dict(masses); mm["EL"] = m
    res = [ns["lammps_point"](g, sym, q, mm, "m") for g in geos]
    if ref is None: ref = res
    de = max(abs(r[0] - b[0]) for r, b in zip(res, ref)); df = max(np.abs(r[1] - b[1]).max() for r, b in zip(res, ref))
    print(f"  data-file El mass {m:5.2f}: max |dE| vs mass 1.0 = {de:.3e} kcal/mol, max |dF| = {df:.3e} kcal/mol/A; bitwise equal: {all((r[1] == b[1]).all() and r[0] == b[0] for r, b in zip(res, ref))}")
print("Part 2: i-PI El mass reaches i-PI; step-0 client forces identical")
RUN = HERE.parent / "run_pimd_lammps.py"; out0 = PROJ.parent / "scratch/tun_A2"; out0.mkdir(exist_ok=True)
f0 = []
for m in (1.0, 0.5, 0.1):
    d = out0 / f"m{m}"
    r = subprocess.run([sys.executable, str(RUN), "--out", str(d), "--beads", "2", "--clients", "2", "--steps", "3", "--mass-el", str(m), "--dump-forces", "--seed", "5"], capture_output=True, text=True)
    if r.returncode: sys.exit("run failed: " + r.stdout[-500:] + r.stderr[-500:])
    dumpm = ns["read_dump"](d / "client0/fd.dump"); f0.append(dumpm[0][1])
    rs = (d / "sim.restart").read_text() if (d / "sim.restart").exists() else (d / "RESTART").read_text()
    mm_ = re.search(r"<m shape='\((\d+)\)'[^>]*>\s*\[(.*?)\]", rs, re.S); mv = [float(x) for x in mm_.group(2).split(",")]
    print(f"  i-PI El mass {m}: restart mass vector (a.u.) El entries {mv[-1]:.5f} (= {mv[-1] / 1822.888486:.4f} amu), Ag entry {mv[0] / 1822.888486:.4f} amu")
print("  step-0 client forces/positions identical across masses: x:", all((f0[0][:, :3] == f[:, :3]).all() for f in f0[1:]), " F:", all((f0[0][:, 3:] == f[:, 3:]).all() for f in f0[1:]))
