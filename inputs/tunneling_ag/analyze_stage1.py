#!/usr/bin/env python3
"""Stage 1 analysis: z-density of the transferring electron (atom index 10, all beads) from i-PI bead trajectories (sim.pos_<k>.xyz).
usage: analyze_stage1.py <run_dir> <P> <dt_fs> <pos_stride> [--skip-ps X] [--out json]
Zones along the donor-acceptor axis: left z < -1.70 | middle | right z > +1.70 (barrier tops at +-1.70 A; start is the left pit, z = -1.975).
Reports: bead-density zone fractions over the production part (with 5-block standard error), first time any bead / the centroid is beyond the top,
fraction of frames with a ring straddling the left top, density histogram (0.05 A bins), and a coarse time series of the left-zone fraction.
"""
import argparse, json
from pathlib import Path
import numpy as np
TOP = 1.70; IA = 10

def zread(path, natom=11):
    out = []; L = Path(path).read_text().splitlines(); n = natom + 2
    for i in range(0, len(L) - n + 1, n): out.append(float(L[i + 2 + IA].split()[3]))
    return np.array(out)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("P", type=int); ap.add_argument("dt", type=float); ap.add_argument("stride", type=int)
    ap.add_argument("--skip-ps", type=float, default=0.0); ap.add_argument("--out", default=""); a = ap.parse_args()
    d = Path(a.run); files = sorted(d.glob("sim.pos_*.xyz"), key=lambda f: int(f.stem.split("_")[1])); assert len(files) == a.P, (len(files), a.P)
    zl = [zread(f) for f in files]; n = min(len(x) for x in zl); z = np.array([x[:n] for x in zl])   # (P, frames); truncated to the shortest file if the run is still going
    t = np.arange(z.shape[1]) * a.dt * a.stride / 1000.0          # ps  (frame 0 = initial structure)
    cen = z.mean(0); sel = t >= a.skip_ps; zz = z[:, sel]
    zone = lambda x: np.array([(x < -TOP).mean(), ((x >= -TOP) & (x <= TOP)).mean(), (x > TOP).mean()])
    fr = zone(zz.ravel()); nb = 5; k = zz.shape[1] // nb
    blocks = np.array([zone(zz[:, i * k:(i + 1) * k].ravel()) for i in range(nb)]); se = blocks.std(0, ddof=1) / np.sqrt(nb)
    first = lambda x: (float(t[np.argmax(x > -TOP)]) if (x > -TOP).any() else None)
    straddle = float(((z.min(0) < -TOP) & (z.max(0) > -TOP))[sel].mean())
    h, e = np.histogram(zz.ravel(), bins=np.arange(-3.2, 3.2 + 1e-9, 0.05), density=True)
    w = max(1, int(round(0.5 / (a.dt * a.stride / 1000.0)))); ts = [float((z[:, i:i + w] < -TOP).mean()) for i in range(0, z.shape[1] - w + 1, w)]
    res = dict(run=str(d), P=a.P, frames=int(z.shape[1]), t_total_ps=float(t[-1]), skip_ps=a.skip_ps, frac_left=float(fr[0]), frac_mid=float(fr[1]), frac_right=float(fr[2]), blocks=blocks.tolist(), se=se.tolist(),
               first_bead_beyond_top_ps=first(z.max(0)), first_centroid_beyond_top_ps=first(cen), straddle_fraction=straddle, mean_z=float(zz.mean()), hist_edges=e.tolist(), hist=h.tolist(), left_frac_vs_time_0p5ps=ts)
    print("%s  P=%d  t=%.1f ps (skip %.1f): left/mid/right = %.4f / %.4f / %.4f  (+-%.4f / %.4f / %.4f)  first bead beyond top %s ps, centroid %s ps, straddling frames %.3f" % (
        d.name, a.P, t[-1], a.skip_ps, *fr, *se, res["first_bead_beyond_top_ps"], res["first_centroid_beyond_top_ps"], straddle))
    print("left-zone fraction vs time (0.5 ps windows):", np.round(ts, 2).tolist())
    if a.out: Path(a.out).write_text(json.dumps(res))
if __name__ == '__main__': main()
