#!/usr/bin/env python3
"""A3: PILE-L equilibration, El radius-of-gyration stationarity, then NVE 2 ps at dt and dt/2 from the same state, for P in {1, 8, 32}.
usage: a3_pipeline.py <work_root> [--P 1 8 32] [--eq-steps 12000] [--clients 1]
Prints/saves per-P: R_g block means, NVE spans of i-PI's printed conserved quantity (per bead = econs/nbeads) in kcal/mol, also x nbeads."""
import argparse, json, subprocess, sys, re
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve(); RUN = HERE.parent / "run_pimd_lammps.py"
EV = 23.060548   # kcal/mol per eV
NEL = 7          # El atoms are the last 7 (Ag7: 7 Ag then 7 El)

def run(args):
    r = subprocess.run([sys.executable, str(RUN)] + [str(x) for x in args], capture_output=True, text=True)
    if r.returncode: sys.exit("FAILED: %s\n%s" % (args, (r.stdout + r.stderr)[-800:]))
    return r.stdout.strip().splitlines()[-1]

def read_xyz_multi(path, natom):
    L = Path(path).read_text().splitlines(); n = natom + 2
    return np.array([[[float(x) for x in L[i + 2 + j].split()[1:4]] for j in range(natom)] for i in range(0, len(L) - n + 1, n)])

def rg_series(d, P, natom=14):
    files = sorted(d.glob("sim.pos_*.xyz"), key=lambda f: int(f.stem.split("_")[1]))   # i-PI zero-pads bead indices when P > 10
    assert len(files) == P, (len(files), P)
    pos = np.array([read_xyz_multi(f, natom) for f in files])   # (P, frames, atoms, 3)
    el = pos[:, :, natom - NEL:, :]
    cen = el.mean(axis=0, keepdims=True)
    return np.sqrt(((el - cen) ** 2).sum(-1).mean(axis=0)).mean(axis=1)   # per frame: mean over El of sqrt(mean_beads |x-c|^2)

def nve_stats(d, P):
    a = np.atleast_2d(np.loadtxt(d / "sim.out")); t, c = a[:, 1], a[:, 2] * EV
    h = len(c) // 2
    return dict(steps=len(c), t_fs=float(t[-1] - t[0]), span_per_bead=float(np.ptp(c)), span_total=float(np.ptp(c) * P),
                span_second_half=float(np.ptp(c[h:])), drift_per_ps=float(np.polyfit(t, c, 1)[0] * 1000), rms=float(c.std()), T_mean=float(a[:, 6].mean()))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("root"); ap.add_argument("--P", type=int, nargs="+", default=[1, 8, 32])
    ap.add_argument("--eq-steps", type=int, default=12000); ap.add_argument("--clients", type=int, default=1); ap.add_argument("--dt", type=float, default=0.25)
    ap.add_argument("--nve-ps", type=float, default=2.0); ap.add_argument("--mass", type=float, default=1.0); a = ap.parse_args()
    root = Path(a.root).resolve(); root.mkdir(parents=True, exist_ok=True); res = {}
    for P in a.P:
        eq = root / f"P{P}_eq"
        if not (eq / "RESTART").exists():
            print(f"P={P} equilibration:", run(["--out", eq, "--beads", P, "--clients", a.clients, "--steps", a.eq_steps, "--dt-fs", a.dt, "--mode", "nvt", "--thermo", "pile_l",
                                              "--tau-fs", 100, "--mass-el", a.mass, "--pos-stride", 10, "--seed", 11]), flush=True)
        r = {}
        if P > 1:
            rg = rg_series(eq, P); k = len(rg) // 6
            r["rg_blocks_A"] = [float(rg[i * k:(i + 1) * k].mean()) for i in range(6)]
        for tag, dt in (("dt", a.dt), ("half", a.dt / 2)):
            d = root / f"P{P}_nve_{tag}"; steps = int(round(a.nve_ps * 1000 / dt))
            if not (d / "sim.out").exists() or len(np.atleast_2d(np.loadtxt(d / "sim.out"))) < steps:
                print(f"P={P} NVE dt={dt}:", run(["--out", d, "--beads", P, "--clients", a.clients, "--steps", steps, "--dt-fs", dt, "--mode", "nve", "--from-restart", eq,
                                                 "--mass-el", a.mass, "--pos-stride", 100, "--stride", 1]), flush=True)
            r[tag] = nve_stats(d, P); r[tag]["dt_fs"] = dt
        res[P] = r
        print(f"P={P}: {json.dumps(r)}", flush=True)
        (root / "a3_results.json").write_text(json.dumps(res, indent=1))
main()
