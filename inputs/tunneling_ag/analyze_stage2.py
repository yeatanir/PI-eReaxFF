#!/usr/bin/env python3
"""Stage 2 analysis: z-density of the quantum electron (atom index 10, all beads of sim.pos_*.xyz) of one or several runs, pooled, versus the 3-D Schroedinger
prediction (data/stage2/*.json from schrod3d_barrier.py: thermal z-marginal for the SAME barrier and mass).
usage: analyze_stage2.py --runs dir:P:dt:stride [dir:P:dt:stride ...] --solver json --mass M [--skip-ps 1] [--nb 8] [--out json] [--label name]
Observables (bead density = exact ring-polymer estimator of the quantum density): rms z, mean z, P(|z|<0.1), P(|z|>0.5), P(|z|>0.75), rho(0)/rho_max
(from 0.05 A histogram, rho(0) = mean of the two bins around 0, rho_max = max of the 3-bin-smoothed histogram), chi2 per bin for |z|<0.8 with block errors.
Errors: standard error over blocks (each run split into nb time blocks; blocks of different runs pooled).
"""
import argparse, json
from pathlib import Path
import numpy as np
IA = 10; EDGES = np.arange(-3.2, 3.2 + 1e-9, 0.05); ZMAX = 3.2     # --zmax restricts data and solver to |z|<zmax (renormalised): classical runs that fall into a pit and cannot leave

def zread(path, natom=11):
    L = Path(path).read_text().splitlines(); n = natom + 2
    return np.array([float(L[i + 2 + IA].split()[3]) for i in range(0, len(L) - n + 1, n)])

def load(d, P):
    d = Path(d); fs = sorted(d.glob("sim.pos_*.xyz"), key=lambda f: int(f.stem.split("_")[1])); assert len(fs) == P, (d, len(fs), P)
    zl = [zread(f) for f in fs]; n = min(len(x) for x in zl); return np.array([x[:n] for x in zl])      # (P, frames)

def obs(z):
    if len(z) < 50: return np.full(6, np.nan)
    h, _ = np.histogram(z, bins=EDGES); h = h / h.sum() / 0.05; c = 0.5 * (EDGES[1:] + EDGES[:-1]); i0 = np.argmin(abs(c))
    sm = np.convolve(h, np.ones(3) / 3, "same"); r0 = 0.5 * (h[i0] + h[i0 - 1]) / sm.max()
    return np.array([np.sqrt((z**2).mean()), z.mean(), (abs(z) < 0.1).mean(), (abs(z) > 0.5).mean(), (abs(z) > 0.75).mean(), r0])
NAMES = ["rms z", "mean z", "P(|z|<0.1)", "P(|z|>0.5)", "P(|z|>0.75)", "rho(0)/rho_max"]

def solver_obs(js, m, key="rho_q"):
    r = [x for x in json.load(open(js)) if abs(x["m"] - m) < 1e-9][0]; zs = np.array(r["zs"]); rho = np.array(r[key]); dz = zs[1] - zs[0]
    # dense resample: bins of 0.05 A by integrating the cell densities (piecewise constant)
    F = np.concatenate([[0], np.cumsum(rho * dz)]); zedge = np.concatenate([zs - dz / 2, [zs[-1] + dz / 2]]); F /= F[-1]
    cdf0 = lambda x: np.interp(x, zedge, F); cdf = lambda x: (cdf0(x) - cdf0(-ZMAX)) / (cdf0(ZMAX) - cdf0(-ZMAX)); h = np.diff(cdf(EDGES)) / 0.05; c = 0.5 * (EDGES[1:] + EDGES[:-1]); i0 = np.argmin(abs(c))
    sm = np.convolve(h, np.ones(3) / 3, "same")
    o = np.array([np.sqrt((h * c**2).sum() * 0.05 / (h.sum() * 0.05)), (h * c).sum() / h.sum(), cdf(0.1) - cdf(-0.1), 1 - (cdf(0.5) - cdf(-0.5)), 1 - (cdf(0.75) - cdf(-0.75)),
                  0.5 * (h[i0] + h[i0 - 1]) / sm.max()])
    return o, h, r

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+", required=True); ap.add_argument("--solver", required=True); ap.add_argument("--mass", type=float, required=True)
    ap.add_argument("--skip-ps", type=float, default=1.0); ap.add_argument("--nb", type=int, default=8); ap.add_argument("--out", default=""); ap.add_argument("--label", default=""); ap.add_argument("--classical", action="store_true", help="compare with the classical 3-D Boltzmann z-marginal (P=1 runs)")
    ap.add_argument("--zmax", type=float, default=3.2); a = ap.parse_args(); global EDGES, ZMAX; ZMAX = a.zmax; EDGES = np.arange(-a.zmax, a.zmax + 1e-9, 0.05); blocks = []; allz = []; info = []
    for spec in a.runs:
        d, P, dt, st = spec.split(":"); P = int(P); dt = float(dt); st = int(st); z = load(d, P); t = np.arange(z.shape[1]) * dt * st / 1000.0
        zz = z[:, t >= a.skip_ps]; k = zz.shape[1] // a.nb; allz.append(zz.ravel()[abs(zz.ravel()) < a.zmax])
        blocks += [obs(zz[:, i * k:(i + 1) * k].ravel()[abs(zz[:, i * k:(i + 1) * k].ravel()) < a.zmax]) for i in range(a.nb)]; info.append((Path(d).name, P, float(t[-1]), obs(zz.ravel()[abs(zz.ravel()) < a.zmax])))
    z = np.concatenate(allz); o = obs(z); B = np.array(blocks); B = B[~np.isnan(B).any(1)]; se = B.std(0, ddof=1) / np.sqrt(len(B)); so, sh, r = solver_obs(a.solver, a.mass, "rho_c" if a.classical else "rho_q")
    print("%s  m=%.1f  runs=%d  samples(all beads)=%d" % (a.label, a.mass, len(a.runs), z.size))
    for n_, x, e, s in zip(NAMES, o, se, so): print("  %-16s PIMD %.5f +- %.5f   solver %.5f   diff %+.1f%% (%+.1f sigma)" % (n_, x, e, s, 100 * (x - s) / abs(s) if s else 0, (x - s) / e if e else 0))
    for nme, P, tt, ob in info: print("  run %-28s P=%d %.1f ps: rms %.4f mean %+.4f P(|z|<.1) %.4f r0 %.3f" % (nme, P, tt, ob[0], ob[1], ob[2], ob[5]))
    hs = [np.histogram(z[:0], bins=EDGES)[0]]
    res = dict(label=a.label, mass=a.mass, n=int(z.size), names=NAMES, pimd=o.tolist(), se=se.tolist(), solver=so.tolist(), runs=[(i[0], i[1], i[2], i[3].tolist()) for i in info],
               hist=(np.histogram(z, bins=EDGES, density=True)[0]).tolist(), solver_hist=sh.tolist(), edges=EDGES.tolist())
    if a.out: Path(a.out).write_text(json.dumps(res))

if __name__ == "__main__": main()
