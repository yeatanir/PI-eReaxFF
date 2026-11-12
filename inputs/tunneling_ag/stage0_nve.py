#!/usr/bin/env python3
"""Stage 0 c/d: conserved quantity of the production setup (Ag fixed, 4 classical El 1 amu, 1 quantum El mass m, P beads) in NVE from a PILE-L-equilibrated state.
usage: stage0_nve.py <root> --P 64 --dt 0.125 [--eq-ps 4] [--nve-ps 2] [--mass 1.0] [--start pit]  -> NVE at dt and dt/2 from the same state
Reports the span of i-PI's printed conserved quantity (per bead = econs/nbeads, kcal/mol; x P = total), the step-to-step ripple, the slow span (25 fs block means), the linear drift."""
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np
RUN = Path(__file__).resolve().parent / "run_barrier_lammps.py"; EV = 23.060548

def run(args):
    r = subprocess.run([sys.executable, str(RUN)] + [str(x) for x in args], capture_output=True, text=True)
    if r.returncode: sys.exit("FAILED %s\n%s" % (args, (r.stdout + r.stderr)[-800:]))
    return r.stdout.strip().splitlines()[-1]

def stats(d, P):
    a = np.atleast_2d(np.loadtxt(d / "sim.out")); t, c = a[:, 1], a[:, 2] * EV
    nb = max(1, int(round(25.0 / (t[1] - t[0])))); k = len(c) // nb; blk = c[:k * nb].reshape(k, nb).mean(1)
    return dict(steps=len(c), span_per_bead=float(np.ptp(c)), span_total=float(np.ptp(c) * P), ripple_rms=float(np.diff(c).std()), slow_span_25fs_blocks=float(np.ptp(blk)),
                drift_kcal_per_ps_per_bead=float(np.polyfit(t, c, 1)[0] * 1000), T_mean=float(a[:, 6].mean()), T_std=float(a[:, 6].std()))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("root"); ap.add_argument("--P", type=int, required=True); ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--eq-ps", type=float, default=4.0); ap.add_argument("--nve-ps", type=float, default=2.0); ap.add_argument("--mass", type=float, default=1.0)
    ap.add_argument("--start", default="pit"); ap.add_argument("--extra", nargs="*", default=[]); ap.add_argument("--seed", type=int, default=11); a = ap.parse_args()
    root = Path(a.root).resolve(); root.mkdir(parents=True, exist_ok=True); eq = root / "eq"
    common = ["--beads", a.P, "--el-mass", a.mass, "--start", a.start] + a.extra
    if not (eq / "RESTART").exists():
        print("equilibration:", run(["--out", eq, "--steps", int(a.eq_ps * 1000 / a.dt), "--dt-fs", a.dt, "--mode", "nvt", "--thermo", "pile_l", "--tau-fs", 100, "--seed", a.seed, "--pos-stride", 40] + common), flush=True)
    res = {}
    for tag, dt in (("dt", a.dt), ("half", a.dt / 2)):
        d = root / f"nve_{tag}"; steps = int(round(a.nve_ps * 1000 / dt))
        if not (d / "sim.out").exists() or len(np.atleast_2d(np.loadtxt(d / "sim.out"))) < steps:
            print(f"NVE dt={dt}:", run(["--out", d, "--steps", steps, "--dt-fs", dt, "--mode", "nve", "--from-restart", eq, "--pos-stride", 200, "--stride", 1] + common), flush=True)
        res[tag] = stats(d, a.P); res[tag]["dt_fs"] = dt; res[tag]["omega_max_dt"] = 2 * a.P * 0.0392865 * dt
        print(tag, json.dumps(res[tag]), flush=True)
    (root / "nve_results.json").write_text(json.dumps(res, indent=1))
main()
