#!/usr/bin/env python3
"""A5: static scan of the transferring electron along the donor-acceptor axis with the exact rebuilt ereac.
Usage: scan_landscape.py <kind pair|ion> <gap_A> <out_json> [zmargin_A=0.6] [step_A=0.1]
Everything else fixed. Needs i-PI importable (for the driver base class) - PYTHONPATH set below.
"""
import json, os, sys, pathlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
PROJ = HERE.parents[1]; ROOT = PROJ.parent
sys.path.insert(0, str(ROOT / "forks/i-pi")); sys.path.insert(0, str(PROJ / "tools")); sys.path.insert(0, str(HERE))
from build_system import write_case
kind, gap, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
zm = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6; step = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
case = ROOT / "scratch/tunneling" / ("scan_%s_%.2f" % (kind, gap))
sym, pos, q, ia, zl, zr = write_case(case, kind, gap)
os.environ.update(EREAC_BIN=str(ROOT / "scratch/phase0/ereac_export_build/ereac"), EREAC_CASE=str(case),
                  EREAC_TMP=str(case / "tmp"), EREAC_WORKERS="2")
import ipi_ereac_driver as d
drv = d.EReac_driver()
zs = np.arange(zl - zm, zr + zm + 1e-9, step)
xs = []
for z in zs:
    p = pos.copy(); p[ia] = [0.001, 0.0, z]            # keep the +0.001 x offset convention
    xs.append(((p + 40.0) / d.BOHR_A).flatten())
res = drv.compute([np.eye(3)] * len(xs), xs)
E = np.array([r[0] for r in res]) * d.KCAL_PER_HARTREE   # kcal/mol
json.dump({"kind": kind, "gap": gap, "z": zs.tolist(), "E_kcal": E.tolist(), "zl": zl, "zr": zr}, open(out, "w"))
i = int(np.argmin(E))
print("%s gap %.2f: %d points, Emin %.3f at z=%.2f; E(left apex) %.3f ; E(right apex) %.3f" % (kind, gap, len(zs), E.min(), zs[i],
      E[int(round((zl - (zl - zm)) / step))], E[int(round((zr - (zl - zm)) / step))]))
