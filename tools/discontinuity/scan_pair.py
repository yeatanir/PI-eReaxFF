"""Static scan of one Ag-Ag pair through a distance: take frame `step` of run dir, move atom a along (a-b) in fine increments, single points."""
import sys, os, subprocess, shutil, numpy as np
from pathlib import Path
sys.path.insert(0, "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build/forks/metal-reax-lammps/tests/lammps")
import run_ereaxff as R
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build"))   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
rd, step, ia, ib = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
lo, hi, n = float(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7]); out = Path(sys.argv[8]); taper = float(sys.argv[9]) if len(sys.argv) > 9 else 5.7
blocks = open(rd / "fd.dump").read().split("ITEM: TIMESTEP\n")[1:]
a = np.array([[float(x) for x in s.split()] for s in blocks[step].split("ITEM: ATOMS")[1].split("\n")[1:] if s.strip()])
pos = a[:, 2:5]; typ = a[:, 1].astype(int)
sym = ["Av" if t == 1 else "El" for t in typ]; q = [1.0 if s == "Av" else -1.0 for s in sym]
out.mkdir(parents=True, exist_ok=True); ff = B / "project/inputs/ag_el/ffield_Ag_El"; shutil.copy(ff, out / "ffield")
masses = R.ffield_masses(str(ff)); elems = R.elements_in_order(sym)
u = pos[ia] - pos[ib]; d0 = np.linalg.norm(u); u /= d0
R.write_data(out / "data.lmp", sym, pos, q, None, False, masses)
plug = B / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"; lmp = B / "external/install_lammps/bin/lmp"
cols = " ".join("c_pr[%d]" % k for k in range(1, 16))
dist = np.linspace(lo, hi, n)
lines = [f"plugin load {os.path.relpath(plug, out)}", "units real", "atom_style charge", "boundary f f f", "newton on", "atom_modify map array", "read_data data.lmp",
         f"pair_style reaxff/metal NULL ereaxff yes electron_taper {taper}", f"pair_coeff * * ffield {' '.join(elems)}", "comm_modify cutoff 12.0",
         "compute pr all pair reaxff/metal", "group mv id %d" % (ia + 1), "velocity all set 0 0 0",
         f"thermo_style custom step pe c_pr[1] c_pr[2] c_pr[3] c_pr[4] c_pr[5] c_pr[6] c_pr[7] c_pr[8] c_pr[9] c_pr[10] c_pr[11] c_pr[12] c_pr[13] c_pr[14] c_pr[15] fx[{ia+1}] fy[{ia+1}] fz[{ia+1}]".replace("fx[", "v_fx").replace("fy[", "v_fy").replace("fz[", "v_fz").replace("]", ""),
         ]
# simple: atom-property variables for force
lines = lines[:-1] + [f"variable fx equal fx[{ia+1}]", f"variable fy equal fy[{ia+1}]", f"variable fz equal fz[{ia+1}]",
         f"thermo_style custom step pe {cols} v_fx v_fy v_fz", "thermo_modify format float %.14g", "thermo 1"]
prev = d0
for k, d in enumerate(dist):
    s = d - prev; prev = d
    lines.append(f"displace_atoms mv move {u[0]*s:.15g} {u[1]*s:.15g} {u[2]*s:.15g} units box")
    lines.append("run 0")
(out / "in.lmp").write_text("\n".join(lines) + "\n")
env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(lmp.parents[1] / "lib"))
r = subprocess.run([str(lmp), "-in", "in.lmp", "-log", "none", "-nocite"], cwd=out, env=env, capture_output=True, text=True)
(out / "stdout.txt").write_text(r.stdout + r.stderr)
rows = [[float(x) for x in l.split()] for l in r.stdout.splitlines() if l.split() and l.split()[0] == "0" and len(l.split()) == 20]
rows = np.array(rows); print("rc", r.returncode, "rows", len(rows), "expected", n)
np.save(out / "scan.npy", np.column_stack([dist[:len(rows)], rows]))
names = "pe eb ea elp emol ev epen ecoa ehb et eco ew ep efi eqeq eelect fx fy fz".split()
for j, nm in enumerate(names):
    col = rows[:, 1 + j]; dd = np.diff(col)
    k = np.argmax(abs(dd))
    print("%-5s max step jump %.3e at d=%.5f (median |step| %.1e)" % (nm, dd[k], dist[k], np.median(abs(dd))))
