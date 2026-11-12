#!/usr/bin/env python3
"""Energy/force/term differences, smooth vs legacy control file, over all phase4_frames_ag frames (charges +1/-1, taper 5.7), one LAMMPS run per control file.
Also verifies that legacy.ctrl == NULL (default). usage: frames_compare.py <outdir>"""
import sys, os, json, glob, subprocess, shutil
from pathlib import Path
import numpy as np
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build")); P = B / "project"   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
sys.path.insert(0, str(B / "forks/metal-reax-lammps/tests/lammps")); import run_ereaxff as R
out = Path(sys.argv[1]).resolve(); out.mkdir(parents=True, exist_ok=True)
ff = P / "inputs/ag_el/ffield_Ag_El"; masses = R.ffield_masses(str(ff))
plug = B / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"; lmp = B / "external/install_lammps/bin/lmp"
env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(lmp.parents[1] / "lib"))
names = "eb ea elp emol ev epen ecoa ehb et eco ew ep efi eqeq eelect".split()
def run(case_dir, frames, ctrl):
    sym, pos, q, cell = R.read_geo(case_dir / "geo"); q = [1.0 if s == "Av" else -1.0 for s in sym]
    elems = R.elements_in_order(sym); d = out / ("w_" + case_dir.name + "_" + Path(ctrl).stem); d.mkdir(exist_ok=True)
    (d / "fd.dump").unlink(missing_ok=True); shutil.copy(ff, d / "ffield"); R.write_data(d / "data.lmp", sym, frames[0]["pos"], q, None, False, masses)
    shutil.copy(ctrl, d / "ctrl.txt") if ctrl != "NULL" else None
    L = [f"plugin load {os.path.relpath(plug, d)}", "units real", "atom_style charge", "boundary f f f", "newton on", "atom_modify map array", "read_data data.lmp",
         f"pair_style reaxff/metal {'ctrl.txt' if ctrl != 'NULL' else 'NULL'} ereaxff yes electron_taper 5.7", f"pair_coeff * * ffield {' '.join(elems)}", "comm_modify cutoff 12.0",
         "compute pr all pair reaxff/metal", 
         "thermo_style custom step pe " + " ".join("c_pr[%d]" % k for k in range(1, 16)), "thermo_modify format float %.14g", "thermo 1"]
    for k, fr in enumerate(frames):
        for i, p in enumerate(fr["pos"]): L.append(f"set atom {i+1} x {p[0]:.15g} y {p[1]:.15g} z {p[2]:.15g}")
        L.append("run 0")
        L.append("write_dump all custom fd.dump id fx fy fz modify sort id format float %.16g append yes")
    (d / "in.lmp").write_text("\n".join(L) + "\n")
    r = subprocess.run([str(lmp), "-in", "in.lmp", "-log", "none", "-nocite"], cwd=d, env=env, capture_output=True, text=True)
    if r.returncode: sys.exit(r.stdout[-800:] + r.stderr[-800:])
    rows = np.array([[float(x) for x in l.split()] for l in r.stdout.splitlines() if len(l.split()) == 17 and l.split()[0].isdigit()])
    blocks = (d / "fd.dump").read_text().split("ITEM: TIMESTEP\n")[1:]
    F = np.array([[[float(x) for x in s.split()[1:4]] for s in b.split("ITEM: ATOMS")[1].split("\n")[1:] if s.strip()] for b in blocks])
    assert len(rows) == len(frames) == len(F), (len(rows), len(frames), len(F))
    return rows[:, 1], rows[:, 2:], F
res = {}
for cd in sorted((P / "reference_data/phase4_frames_ag").glob("Ag*_md")):
    fr = json.loads((cd / "frames.json").read_text())["frames"]
    for c in ("NULL", str(P / "inputs/control/legacy.ctrl"), str(P / "inputs/control/smooth.ctrl")):
        res[(cd.name, Path(c).stem)] = run(cd, fr, c)
print("legacy.ctrl vs NULL (max |dE|, max |dF|, max |d term|):")
tot = {"dE": [], "dF": [], "dT": []}
print("%-9s %5s | %-30s | %-30s" % ("case", "n", "legacy vs NULL: dE dF dTerm", "smooth-legacy: max|dE| rms dE | max|dF| rms dF"))
allT = []
for cd in sorted(set(k[0] for k in res)):
    e0, t0, f0 = res[(cd, "NULL")]; e1, t1, f1 = res[(cd, "legacy")]; e2, t2, f2 = res[(cd, "smooth")]
    dE = e2 - e1; dF = (f2 - f1).reshape(len(e1), -1); dT = t2 - t1
    tot["dE"] += list(dE); tot["dF"] += list(np.sqrt((dF ** 2).sum(1))); allT.append(dT)
    print("%-9s %5d | %.1e %.1e %.1e | dE max %.4f rms %.4f | dF max %.4f rms %.4f" % (cd, len(e1), abs(e1 - e0).max(), abs(f1 - f0).max(), abs(t1 - t0).max(), abs(dE).max(), np.sqrt((dE ** 2).mean()), abs(dF).max(), np.sqrt((dF ** 2).mean())))
dE = np.array(tot["dE"]); T = np.vstack(allT)
nF = sum(res[(c, 'smooth')][2].size for c in set(k[0] for k in res))
dFall = np.concatenate([(res[(c, "smooth")][2] - res[(c, "legacy")][2]).ravel() for c in set(k[0] for k in res)])
print("ALL %d frames: smooth-legacy dE max %.5f rms %.5f kcal/mol ; dF component max %.5f rms %.5f kcal/mol/A" % (len(dE), abs(dE).max(), np.sqrt((dE ** 2).mean()), abs(dFall).max(), np.sqrt((dFall ** 2).mean())))
for j, n in enumerate(names):
    if abs(T[:, j]).max() > 0: print("  term %-5s max %.5f rms %.5f" % (n, abs(T[:, j]).max(), np.sqrt((T[:, j] ** 2).mean())))

np.savez(out / "frames_compare.npz", **{f"{k[0]}__{k[1]}__{n}": v for k, r in res.items() for n, v in zip(("E", "T", "F"), r)})
