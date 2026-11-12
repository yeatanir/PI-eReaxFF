#!/usr/bin/env python3
"""Per-step NVE through LAMMPS+Metal-Reax (read-only use of shared lmp/plugin). Dumps x,f every step; thermo every step with 15 pair terms."""
import sys, os, json, subprocess, shutil
from pathlib import Path
R_DIR = Path(__file__).resolve().parents[3] / "forks/metal-reax-lammps/tests/lammps"; sys.path.insert(0, str(R_DIR))   # relative to the tree (was a hard-coded MacBook path)
import run_ereaxff as R
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build"))   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
case, out, steps, dt, T, seed = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
frame = int(sys.argv[7]) if len(sys.argv) > 7 else 0
taper = 10.0 if case == "GRAPHENE" else 5.7
cd = B / "project/reference_data/phase4_frames" if case == "GRAPHENE" else B / "project/reference_data/phase4_frames_ag" / case
geo = cd / "geo"; sym, pos, q, cell = R.read_geo(geo)
if q == [0.0] * len(q): q = [1.0 if s == "Av" else -1.0 for s in sym]
fr = json.loads((cd / "frames.json").read_text())["frames"]
p0 = fr[frame]["pos"]
ff = (cd / "ffield") if case == "GRAPHENE" else B / "project/inputs/ag_el/ffield_Ag_El"
out.mkdir(parents=True, exist_ok=True); shutil.copy(ff, out / "ffield")
masses = R.ffield_masses(str(ff)); elems = R.elements_in_order(sym)
(out / "ctrl.txt").write_text(os.environ.get("CTRLTXT") or ("nbrhood_cutoff %s\n" % os.environ.get("BCUT","5.0")))
PER = case == "GRAPHENE"
R.write_data(out / "data.lmp", sym, p0, q, cell if PER else None, PER, masses)
plug = B / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"
lmp = B / "external/install_lammps/bin/lmp"
cols = " ".join("c_pr[%d]" % k for k in range(1, 16))
(out / "in.lmp").write_text(f"""plugin load {os.path.relpath(plug, out)}
units real
atom_style charge
boundary {"p p p" if PER else "f f f"}
newton on
read_data data.lmp
pair_style reaxff/metal {os.environ.get('CTRL','NULL')} ereaxff yes electron_taper {taper}
pair_coeff * * ffield {' '.join(elems)}
comm_modify cutoff 12.0
compute pr all pair reaxff/metal
velocity all create {T} {seed} dist gaussian mom yes rot no
timestep {dt}
fix i all nve
dump fd all custom 1 fd.dump id type x y z fx fy fz
dump_modify fd format float %.16g sort id
thermo_style custom step pe ke etotal {cols}
thermo_modify format float %.12g
thermo 1
run {steps}
""")
env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(lmp.parents[1] / "lib"))
r = subprocess.run([str(lmp), "-in", "in.lmp", "-log", "none", "-nocite"], cwd=out, env=env, capture_output=True, text=True)
(out / "stdout.txt").write_text(r.stdout + r.stderr)
(out / "control_used.txt").write_text("control=%s\n%s" % (os.environ.get("CTRL", "NULL"), (out / "ctrl.txt").read_text() if os.environ.get("CTRL") else ""))
print("CONTROL_USED", os.environ.get("CTRL", "NULL"), (out / "ctrl.txt").read_text().replace("\n", " ") if os.environ.get("CTRL") else "")
print("rc", r.returncode, r.stdout[-300:] if r.returncode else "ok")
