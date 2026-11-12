#!/usr/bin/env python3
"""Detailed exact-ereac landscape of the transferring electron for the symmetric Ag3...Ag3 pair (revised Part B).
Usage: landscape_gap.py <gap_A> <out_json>
 A  full axis profile (step 0.025), both sides (mirror check)
 B  barrier-top refinement (step 0.005) on both sides
 C  curvature at the top (axial, x, y), from the refined top
 D  2-D (x,y) cuts at several z: does the electron have lower-energy paths off axis?
 E  relaxation of the OTHER four electrons (12 dof, L-BFGS on exact ereac energy+force) with the transferring electron held at z = pit, top, mid
x = out-of-plane, y = in the Ag3 plane (Ag3 units lie in the x = 0 plane), z = donor-acceptor axis.
"""
import json, os, sys, pathlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; PROJ = HERE.parents[1]; ROOT = PROJ.parent
sys.path.insert(0, str(ROOT / "forks/i-pi")); sys.path.insert(0, str(PROJ / "tools")); sys.path.insert(0, str(HERE))
from build_system import write_case
gap, out = float(sys.argv[1]), sys.argv[2]
case = ROOT / "scratch/tunneling" / ("land_%.2f" % gap)
sym, pos, q, ia, zl, zr = write_case(case, "pair", gap)
os.environ.update(EREAC_BIN=str(ROOT / "scratch/phase0/ereac_export_build/ereac"), EREAC_CASE=str(case), EREAC_TMP=str(case / "tmp"), EREAC_WORKERS="2")
import ipi_ereac_driver as d
drv = d.EReac_driver(); KC = d.KCAL_PER_HARTREE; B = d.BOHR_A

def energies(points):                      # points: (n,3) positions of the transferring electron (A, system frame without the 40 A shift)
    xs = []
    for p in points:
        q_ = pos.copy(); q_[ia] = p; xs.append(((q_ + 40.0) / B).flatten())
    return np.array([r[0] for r in drv.compute([np.eye(3)] * len(xs), xs)]) * KC

R = {"gap": gap, "zl": zl, "zr": zr}
# A
zs = np.round(np.arange(-3.0, 3.0 + 1e-9, 0.025), 4)
R["A_z"] = zs.tolist(); R["A_E"] = energies([[0.001, 0, z] for z in zs]).tolist(); print("A done", flush=True)
# B
topz = {}
for side, (a, b) in (("L", (zl + 0.05, zl + 0.8)), ("R", (zr - 0.8, zr - 0.05))):
    zz = np.round(np.arange(a, b + 1e-9, 0.005), 4); E = energies([[0.001, 0, z] for z in zz])
    sel = (zz > zl + 0.1) & (zz < zl + 0.7) if side == "L" else (zz > zr - 0.7) & (zz < zr - 0.1)
    i = np.where(sel)[0][np.argmax(E[sel])]; topz[side] = float(zz[i]); R["B_" + side] = {"z": zz.tolist(), "E": E.tolist()}
R["top_z"] = topz; print("B done, tops", topz, flush=True)
# C
disp = np.round(np.arange(-0.08, 0.08 + 1e-9, 0.02), 4); R["C"] = {}
for side in ("L", "R"):
    z0 = topz[side]
    for ax, v in (("z", [0, 0, 1]), ("x", [1, 0, 0]), ("y", [0, 1, 0])):
        E = energies([np.array([0.001, 0, z0]) + u * np.array(v) for u in disp]); R["C"][side + "_" + ax] = E.tolist()
R["C_disp"] = disp.tolist(); print("C done", flush=True)
# D
g = np.round(np.arange(-0.6, 0.6 + 1e-9, 0.15), 3); R["D_grid"] = g.tolist(); R["D"] = {}
for z in (zl, zl + 0.1, zl + 0.2, topz["L"], zl + 0.5, zl + 0.7, zl + 1.0, 0.0):
    pts = [[0.001 + x, y, z] for x in g for y in g]
    R["D"]["%.3f" % z] = np.array(energies(pts)).reshape(len(g), len(g)).tolist()
print("D done", flush=True)
# E: relax the other electrons (all El except ia) with the transferring electron fixed at (0.001, 0, z)
from scipy.optimize import minimize
others = [k for k, s in enumerate(sym) if s == "El" and k != ia]
base = pos.copy(); R["E"] = {}
def fun(v, zt):
    p = base.copy(); p[ia] = [0.001, 0, zt]; p[others] = v.reshape(-1, 3)
    e, f, _, _ = drv.compute_structure(np.eye(3), ((p + 40.0) / B).flatten())
    return e * KC, (-f.reshape(-1, 3)[others] * (KC / B)).flatten()   # dE/dx (kcal/mol/A) = -F ; F in Ha/bohr -> kcal/mol/A: *KC/B
for name, zt in (("pit", zl), ("top", topz["L"]), ("mid", 0.0)):
    v0 = base[others].flatten(); e0 = fun(v0, zt)[0]
    res = minimize(fun, v0, args=(zt,), jac=True, method="L-BFGS-B", options=dict(maxiter=60, gtol=1e-3))
    mv = np.abs(res.x.reshape(-1, 3) - base[others]).max()
    R["E"][name] = {"z": zt, "E_fixed": e0, "E_relaxed": float(res.fun), "max_move_A": float(mv), "nit": int(res.nit)}
    print("E", name, R["E"][name], flush=True)
json.dump(R, open(out, "w"))
