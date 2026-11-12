#!/usr/bin/env python3
"""Time per force call, legacy vs smooth control, on a fixed geometry (tiny timestep so positions barely move). usage: time_ctrl.py <geo> <ffield> <taper> <periodic 0|1> <steps> <out>"""
import sys, os, subprocess, shutil, re
from pathlib import Path
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build")); P = B / "project"   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
sys.path.insert(0, str(B / "forks/metal-reax-lammps/tests/lammps")); import run_ereaxff as R
geo, ffp, taper, per, steps, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4] == "1", int(sys.argv[5]), Path(sys.argv[6]).resolve()
sym, pos, q, cell = R.read_geo(geo); masses = R.ffield_masses(str(ffp)); elems = R.elements_in_order(sym)
plug = B / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"; lmp = B / "external/install_lammps/bin/lmp"
env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(lmp.parents[1] / "lib"))
for mode in ("legacy", "smooth"):
    d = out / mode; d.mkdir(parents=True, exist_ok=True); shutil.copy(ffp, d / "ffield"); shutil.copy(P / f"inputs/control/{mode}.ctrl", d / "ctrl.txt")
    R.write_data(d / "data.lmp", sym, pos, q, cell, per, masses)
    (d / "in.lmp").write_text(f"""plugin load {os.path.relpath(plug, d)}
units real
atom_style charge
boundary {'p p p' if per else 'f f f'}
newton on
read_data data.lmp
pair_style reaxff/metal ctrl.txt ereaxff yes electron_taper {taper}
pair_coeff * * ffield {' '.join(elems)}
comm_modify cutoff 12.0
timestep 0.00001
fix i all nve
thermo {steps}
run {steps}
""")
    r = subprocess.run([str(lmp), "-in", "in.lmp", "-log", "none", "-nocite"], cwd=d, env=env, capture_output=True, text=True)
    if r.returncode: sys.exit(r.stdout[-600:] + r.stderr[-300:])
    m = re.search(r"Loop time of ([\d.e+-]+) on 1 procs for (\d+) steps with (\d+) atoms", r.stdout)
    pe = [l.split() for l in r.stdout.splitlines() if l.split() and l.split()[0] == "0" and len(l.split()) >= 3]
    print(f"CONTROL_USED {mode}: {len(sym)} atoms, {float(m.group(1))/int(m.group(2))*1e3:.3f} ms per force call ({m.group(2)} steps)  pe0={pe[0][2] if pe else '?'}")
