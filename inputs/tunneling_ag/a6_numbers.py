#!/usr/bin/env python3
"""A5/A6 numbers from scan_fine.py JSON: curvature (quadratic fit, window +-win A) -> k (eV/A^2) -> omega per mass -> P, dt, free-ring radius.
Usage: a6_numbers.py <fine.json> [window_A=0.08]"""
import json, sys
import numpy as np
from analyze_landscape import KCAL_EV, AMU, HBAR, KB, T, omega, req_P, rfree
d = json.load(open(sys.argv[1])); win = float(sys.argv[2]) if len(sys.argv) > 2 else 0.08
u = np.array(d["disp"]); print("gap %.2f  (electron displacement scans, quadratic fit |u|<=%.2f A)" % (d["gap"], win))
rows = {}
for name, E in d["E_kcal"].items():
    E = np.array(E); m = np.abs(u) <= win + 1e-9; c = np.polyfit(u[m], E[m], 2)
    k = 2 * c[0] * KCAL_EV; ub = -c[1] / (2 * c[0]) if c[0] != 0 else np.nan
    rows[name] = k
    asym = (E[np.argmin(np.abs(u - 0.2))] - E[np.argmin(np.abs(u + 0.2))])
    print("  %-8s k = %8.3f eV/A^2  (minimum offset %+.3f A, E(+0.2)-E(-0.2) = %+.2f kcal/mol)" % (name, k, ub, asym))
for site in ("apex", "mid"):
    ks = [rows[site + "_" + a] for a in "zxy"]
    kmax = max(ks)
    print("  %s: stiffest axis k = %.2f eV/A^2" % (site, kmax))
    if kmax > 0:
        for m in (1.0, 0.5, 0.1):
            w = omega(kmax, m); P = req_P(w)
            print("     m=%.1f amu: omega=%.3f /fs (%.0f cm^-1), beta*hbar*omega=%.1f -> P >= %d for <1%%; dt <= 0.2/omega = %.3f fs" % (m, w, w * 5308.8, HBAR * w / (KB * T), P, 0.2 / w))
print("  free-ring rms radius about the centroid at 300 K (3-D): " + ", ".join("%.2f A (%.1f amu)" % (rfree(m), m) for m in (1.0, 0.5, 0.1)))
