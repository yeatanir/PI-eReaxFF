#!/usr/bin/env python3
"""A5 fine scan: displace the transferring electron from its paired position at the LEFT apex (and from mid-gap) along z, x, y.
Usage: scan_fine.py <kind> <gap> <out_json> [half_range_A=0.3] [step_A=0.02]"""
import json, os, sys, pathlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; PROJ = HERE.parents[1]; ROOT = PROJ.parent
sys.path.insert(0, str(ROOT / "forks/i-pi")); sys.path.insert(0, str(PROJ / "tools")); sys.path.insert(0, str(HERE))
from build_system import write_case
kind, gap, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
hr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3; st = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02
case = ROOT / "scratch/tunneling" / ("fine_%s_%.2f" % (kind, gap))
sym, pos, q, ia, zl, zr = write_case(case, kind, gap)
os.environ.update(EREAC_BIN=str(ROOT / "scratch/phase0/ereac_export_build/ereac"), EREAC_CASE=str(case), EREAC_TMP=str(case / "tmp"), EREAC_WORKERS="2")
import ipi_ereac_driver as d
drv = d.EReac_driver()
disp = np.arange(-hr, hr + 1e-9, st); res = {}
for name, origin in (("apex", np.array([0.001, 0.0, zl])), ("mid", np.array([0.001, 0.0, 0.0]))):
    for ax, dirv in (("z", [0, 0, 1]), ("x", [1, 0, 0]), ("y", [0, 1, 0])):
        xs = []
        for u in disp:
            p = pos.copy(); p[ia] = origin + u * np.array(dirv); xs.append(((p + 40.0) / d.BOHR_A).flatten())
        E = np.array([r[0] for r in drv.compute([np.eye(3)] * len(xs), xs)]) * d.KCAL_PER_HARTREE
        res[name + "_" + ax] = E.tolist()
json.dump({"kind": kind, "gap": gap, "disp": disp.tolist(), "E_kcal": res}, open(out, "w")); print("done", out)
