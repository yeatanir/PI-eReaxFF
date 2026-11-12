#!/usr/bin/env python3
"""Higher barrier (V0 = 8 kcal/mol, sigma = 0.10 A) PIMD analysis: runs vs 3-D solver, equilibration (left vs right start), straddling, crossings.
usage: analyze_higher_barrier.py <rundir_with_hb_*_dirs> <out_json> <out_npz>
Per-run: bead-pooled z density (0.05 A bins, frames with centroid |z| > 1.2 A excluded = the pit rule), rms z, rho(0)/rho_max, P(|z|<0.1),
fraction of frames with beads on both sides of z = 0 (straddle), centroid time on each side, centroid sign changes per ps, first centroid crossing,
conserved quantity std/drift (sim.out col 3, per bead). Run-to-run spread = std/sqrt(n) over independent runs (not block errors)."""
import json, re, sys
from pathlib import Path
import numpy as np

KT = 0.0258520          # eV at 300 K
SKIP_PS = 2.0; ZMAX = 1.5; PIT = 1.2
EDGES = np.arange(-ZMAX, ZMAX + 1e-9, 0.05); C = 0.5 * (EDGES[1:] + EDGES[:-1])
# name -> (mass, P, dt_fs, pos_stride, start sign)
SETUP = {"m1": (1.0, 72, 0.125, 200), "m05": (0.5, 100, 0.09, 200), "cl": (1.0, 1, 0.125, 400)}
SOLVER = Path(__file__).resolve().parents[2] / "validation/tunneling_ag/data/higher_barrier.json"


def solver_hist(m, key):
    r = [x for x in json.load(open(SOLVER)) if x["V0"] == 8.0 and abs(x["sigma"] - 0.10) < 1e-9 and abs(x["m"] - m) < 1e-9][0]
    zs = np.array(r["zs"]); rho = np.array(r[key]); dz = zs[1] - zs[0]
    F = np.concatenate([[0], np.cumsum(rho * dz)]); ze = np.concatenate([zs - dz / 2, [zs[-1] + dz / 2]]); F /= F[-1]
    cdf0 = lambda x: np.interp(x, ze, F); cdf = lambda x: (cdf0(x) - cdf0(-ZMAX)) / (cdf0(ZMAX) - cdf0(-ZMAX))
    return np.diff(cdf(EDGES)) / 0.05, cdf


def obs_from_hist(h):
    i0 = np.argmin(abs(C)); sm = np.convolve(h, np.ones(3) / 3, "same")
    p01 = h[abs(C) < 0.1].sum() * 0.05
    return dict(rms=float(np.sqrt((h * C**2).sum() / h.sum())), p01=float(p01), r0=float(0.5 * (h[i0] + h[i0 - 1]) / sm.max()))


def one_run(d, mass, P, dt_fs, stride):
    ring = np.load(d / "ring.npy")[..., 2].astype(float)             # (P, F)
    assert ring.shape[0] == P, (d, ring.shape)
    dtf = stride * dt_fs / 1000.0; t = np.arange(ring.shape[1]) * dtf; k = t >= SKIP_PS
    z = ring[:, k]; c = z.mean(0); tt = t[k]; keep = abs(c) <= PIT
    zz = z[:, keep].ravel(); zz = zz[abs(zz) < ZMAX]
    h = np.histogram(zz, bins=EDGES, density=True)[0]
    o = obs_from_hist(h) if zz.size > 100 else {}
    sgn = np.sign(c); sgn[sgn == 0] = 1; ch = np.nonzero(np.diff(sgn))[0]
    first = float(tt[ch[0] + 1] - tt[0]) if len(ch) else None
    straddle = float(((z.min(0) < 0) & (z.max(0) > 0)).mean())
    # conserved quantity
    so = np.loadtxt(d / "sim.out"); ts = so[:, 1] / 1000.0; cons = so[:, 2][ts >= SKIP_PS]; tsc = ts[ts >= SKIP_PS]
    sl = np.polyfit(tsc, cons, 1)[0] if len(cons) > 10 else np.nan
    log = (d.parent / (d.name + ".log")).read_text() if (d.parent / (d.name + ".log")).exists() else ""
    m = re.search(r"after ([0-9.]+) s \(P=\d+, \d+ clients, (\d+) steps\)", log)
    sps = float(m.group(2)) / float(m.group(1)) if m else None
    return dict(run=d.name, mass=mass, P=P, dt_fs=dt_fs, frame_ps=dtf, T_ps=float(t[-1]), n_frames=int(k.sum()), frac_pit_excluded=float(1 - keep.mean()),
                obs=o, straddle=straddle, frac_left=float((c < 0).mean()), frac_right=float((c > 0).mean()), start_side=("L" if "_L_" in d.name else "R"),
                mean_c=float(c.mean()), crossings=int(len(ch)), crossings_per_ps=float(len(ch) / (tt[-1] - tt[0])), first_crossing_ps=first,
                cons_std_kT=float(cons.std() / KT), cons_drift_kT_per_ns=float(sl * 1000 / KT), steps_per_s=sps, hist=h.tolist()), c, tt


def verdict(runs):
    """Equilibration rule fixed in advance: equilibrated only if (a) every run has >= 10 centroid crossings after the 2 ps skip, (b) every run first crosses
    within 25 % of its length, (c) left-start and right-start means of rho(0)/rho_max and of the straddle fraction agree within 20 % (relative).
    Otherwise: NOT equilibrated; report each start side separately (leakage into the barrier), not as an equilibrium profile."""
    why = []
    for r in runs:
        if r["crossings"] < 10: why.append("%s: only %d centroid crossings" % (r["run"], r["crossings"]))
        f = r["first_crossing_ps"]
        if f is None or f > 0.25 * (r["T_ps"] - SKIP_PS): why.append("%s: first crossing %s ps (run %.1f ps)" % (r["run"], "never" if f is None else "%.2f" % f, r["T_ps"]))
    for name, g in (("r0", lambda r: r["obs"].get("r0", np.nan)), ("straddle", lambda r: r["straddle"])):
        L = [g(r) for r in runs if r["start_side"] == "L"]; R = [g(r) for r in runs if r["start_side"] == "R"]
        if L and R:
            a, b = float(np.mean(L)), float(np.mean(R))
            if not (abs(a - b) <= 0.2 * max(abs(a), abs(b))): why.append("%s: L-start %.3f vs R-start %.3f differ by > 20 %%" % (name, a, b))
    return dict(equilibrated=not why, reasons=why)


def main():
    root = Path(sys.argv[1]); res = {}; npz = {}
    for d in sorted(root.glob("hb_*")):
        if not d.is_dir() or not (d / "ring.npy").exists(): continue
        key = d.name.split("_")[1]; mass, P, dt, st = SETUP[key]
        r, c, tt = one_run(d, mass, P, dt, st); res.setdefault(key, []).append(r); npz[d.name + "_t"] = tt.astype(np.float32); npz[d.name + "_c"] = c.astype(np.float32)
    out = {}
    for key, runs in res.items():
        mass = runs[0]["mass"]; hs = np.array([r["hist"] for r in runs]); hm = hs.mean(0)
        sh, cdf = solver_hist(mass, "rho_c" if key == "cl" else "rho_q"); so = obs_from_hist(sh)
        per = {n: np.array([r["obs"][n] for r in runs if r["obs"]]) for n in ("rms", "p01", "r0")}
        summ = {n: dict(mean=float(v.mean()), spread=float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None, per_run=v.tolist(), solver=so[n]) for n, v in per.items() if len(v)}
        out[key] = dict(mass=mass, runs=runs, verdict=verdict(runs), summary=summ, hist_mean=hm.tolist(), hist_per_run=hs.tolist(), solver_hist=sh.tolist(), edges=EDGES.tolist(),
                        solver_obs=so, solver_P_left=float(cdf(0.0)))
    Path(sys.argv[2]).write_text(json.dumps(out)); np.savez_compressed(sys.argv[3], **npz)
    for key, o in out.items():
        print(key, "mass", o["mass"], "runs", len(o["runs"]))
        for r in o["runs"]:
            print("  %-18s T=%.1f ps start %s straddle %.3f  L/R %.2f/%.2f  crossings %d (%.2f/ps) first %s  rms %.4f r0 %.3f  cons std %.3f kT drift %.4f kT/ns  %s steps/s" %
                  (r["run"], r["T_ps"], r["start_side"], r["straddle"], r["frac_left"], r["frac_right"], r["crossings"], r["crossings_per_ps"], r["first_crossing_ps"], r["obs"].get("rms", np.nan), r["obs"].get("r0", np.nan), r["cons_std_kT"], r["cons_drift_kT_per_ns"], r["steps_per_s"]))
        print("  VERDICT:", "equilibrated" if o["verdict"]["equilibrated"] else "NOT equilibrated", o["verdict"]["reasons"])
        for n, s in o["summary"].items(): print("  %s: PIMD %.4f +- %s (run-to-run) solver %.4f" % (n, s["mean"], ("%.4f" % s["spread"]) if s["spread"] is not None else "n/a", s["solver"]))


if __name__ == "__main__": main()
