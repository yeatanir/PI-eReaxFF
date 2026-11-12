#!/usr/bin/env python3
"""Stage 1 table: bead-density moments of the transferring electron along z for several runs vs the 3-D Schroedinger prediction (frozen other electrons).
usage: stage1_report.py <schrod3d json> <mass> <skip_ps> <run_dir:P:dt:stride:label> ...   (5-block standard errors over the production part)"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent)); from analyze_stage1 import zread   # noqa  (module runs main on import: guard below)
sc = json.load(open(sys.argv[1]))[sys.argv[2]]; zs = np.array(sc["zs"]); rho = np.array(sc["rho_z"]); skip = float(sys.argv[3])
def obs(z):   # z: array of bead z values
    return np.array([np.sqrt((z**2).mean()), (np.abs(z) > 0.5).mean(), (np.abs(z) > 0.75).mean(), (np.abs(z) > 1.0).mean(), (z < -1.7).mean(), (np.abs(z) > 1.7).mean()])
names = ["rms z (A)", "P(|z|>0.5)", "P(|z|>0.75)", "P(|z|>1.0)", "P(z<-1.7) left zone", "P(|z|>1.7)"]
sch = np.array([np.sqrt((rho * zs**2).sum()), rho[np.abs(zs) > 0.5].sum(), rho[np.abs(zs) > 0.75].sum(), rho[np.abs(zs) > 1.0].sum(), rho[zs < -1.7].sum(), rho[np.abs(zs) > 1.7].sum()])
print("%-34s" % "run", "  ".join("%-17s" % n for n in names))
print("%-34s" % ("Schroedinger 3-D, %s amu" % sys.argv[2]), "  ".join("%-17.4g" % v for v in sch))
for spec in sys.argv[4:]:
    d, P, dt, stride, label = spec.split(":"); P = int(P); dt = float(dt); stride = int(stride)
    files = sorted(Path(d).glob("sim.pos_*.xyz"), key=lambda f: int(f.stem.split("_")[1])); zl = [zread(f) for f in files]; n = min(len(x) for x in zl); z = np.array([x[:n] for x in zl]); t = np.arange(z.shape[1]) * dt * stride / 1000
    z = z[:, t >= skip]; k = z.shape[1] // 5; B = np.array([obs(z[:, i * k:(i + 1) * k].ravel()) for i in range(5)]); m = obs(z.ravel()); se = B.std(0, ddof=1) / np.sqrt(5)
    print("%-34s" % label, "  ".join("%-8.4g+-%-7.2g" % (a, b) for a, b in zip(m, se)), "  [%.1f ps used]" % (t[-1] - skip))
