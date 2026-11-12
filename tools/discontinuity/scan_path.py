"""Single points along the straight line between two consecutive MD frames; locate the energy jump and the pair that crosses a cutoff there.
usage: scan_path.py run_dir step N out_dir [BCUT]"""
import sys, os, subprocess, shutil, numpy as np
from pathlib import Path
sys.path.insert(0, "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build/forks/metal-reax-lammps/tests/lammps")
import run_ereaxff as R
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build"))   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
rd, step, N, out = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]); bcut = sys.argv[5] if len(sys.argv) > 5 else "5.0"
bl = open(rd / "fd.dump").read().split("ITEM: TIMESTEP\n")[1:]
def get(i):
    a = np.array([[float(x) for x in s.split()] for s in bl[i].split("ITEM: ATOMS")[1].split("\n")[1:] if s.strip()]); return a[:, 2:5], a[:, 1].astype(int)
p0, typ = get(step); p1, _ = get(step + 1)
sym = ["Av" if t == 1 else "El" for t in typ]; q = [1.0 if s == "Av" else -1.0 for s in sym]
out.mkdir(parents=True, exist_ok=True); ff = B / "project/inputs/ag_el/ffield_Ag_El"; shutil.copy(ff, out / "ffield"); (out / "ctrl.txt").write_text("nbrhood_cutoff %s\n" % bcut)
masses = R.ffield_masses(str(ff)); elems = R.elements_in_order(sym)
R.write_data(out / "data.lmp", sym, p0.tolist(), q, None, False, masses)
plug = B / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"; lmp = B / "external/install_lammps/bin/lmp"
cols = " ".join("c_pr[%d]" % k for k in range(1, 16))
S = np.linspace(0, 1, N + 1)
lines = [f"plugin load {os.path.relpath(plug, out)}", "units real", "atom_style charge", "boundary f f f", "newton on", "atom_modify map array", "read_data data.lmp",
         f"pair_style reaxff/metal ctrl.txt ereaxff yes electron_taper 5.7", f"pair_coeff * * ffield {' '.join(elems)}", "comm_modify cutoff 12.0", "compute pr all pair reaxff/metal",
         f"thermo_style custom step pe {cols}", "thermo_modify format float %.14g", "thermo 1"]
prev = p0.copy()
for s in S:
    p = p0 + s * (p1 - p0); d = p - prev; prev = p
    for i in range(len(p)):
        lines.append(f"set atom {i+1} x {p[i,0]:.15g} y {p[i,1]:.15g} z {p[i,2]:.15g}")
    lines.append("run 0")
(out / "in.lmp").write_text("\n".join(lines) + "\n")
env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(lmp.parents[1] / "lib"))
r = subprocess.run([str(lmp), "-in", "in.lmp", "-log", "none", "-nocite"], cwd=out, env=env, capture_output=True, text=True)
rows = np.array([[float(x) for x in l.split()] for l in r.stdout.splitlines() if len(l.split()) == 17 and l.split()[0] == "0"])
print("rc", r.returncode, "rows", len(rows))
names = "pe eb ea elp emol ev epen ecoa ehb et eco ew ep efi eqeq eelect".split()
dE = np.diff(rows[:, 1]); k = np.argmax(abs(dE)); med = np.median(abs(dE))
print("E(s=0)=%.6f E(s=1)=%.6f total change %.5f ; biggest single-substep change %.5f at s in [%.4f,%.4f] (median substep %.1e)" % (rows[0, 1], rows[-1, 1], rows[-1, 1] - rows[0, 1], dE[k], S[k], S[k + 1], med))
for j, nm in enumerate(names):
    dd = np.diff(rows[:, 1 + j]); kk = np.argmax(abs(dd))
    if abs(dd[kk]) > 5 * max(np.median(abs(dd)), 1e-9) and abs(dd[kk]) > 1e-4: print("  %-5s jump %+.5f at s-bin %d (median %.1e)" % (nm, dd[kk], kk, np.median(abs(dd))))
# which pairs cross, in the bin k
ps = p0 + S[k] * (p1 - p0); pe_ = p0 + S[k + 1] * (p1 - p0)
D0 = np.linalg.norm(ps[:, None] - ps[None], axis=-1); D1 = np.linalg.norm(pe_[:, None] - pe_[None], axis=-1)
nm_ = {1: "Ag", 2: "El"}
print("pair distances in jump bin (moving pairs with |dd| > 1e-6, sorted by nearest 'round' value):")
for a in range(len(p0)):
    for b in range(a + 1, len(p0)):
        if abs(D1[a, b] - D0[a, b]) > 1e-7:
            print("   %s%d-%s%d  %.6f -> %.6f" % (nm_[typ[a]], a, nm_[typ[b]], b, D0[a, b], D1[a, b]))
