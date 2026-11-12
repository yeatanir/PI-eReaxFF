#!/usr/bin/env python3
"""1-D Schroedinger solve on the measured axis profile V(z) (LAMMPS + Metal-Reax fine scan, transferring electron at (0.001, 0, z), everything else frozen).
Per mass: lowest levels, tunnel splitting of the two pit ground states, 300 K thermal density and zone fractions (left z<-1.7 | middle | right z>1.7),
ground-state / pit-state probability beyond the barrier top, harmonic-pit comparison. Usage: schrod1d.py <E_F.npy> <zf.npy> <out.json>
"""
import json, sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
HB2_2M = 0.0481741; KT = 0.0019872041 * 300.0
E = np.load(sys.argv[1]); z = np.load(sys.argv[2]); V = E - E.min(); h = z[1] - z[0]; top = 1.70
res = {}
# harmonic comparison: pit stiffness 45 eV/A^2 = 1037.7 kcal/mol/A^2 at 0.275 A from the minimum
print("harmonic pit energy at 0.275 A: %.1f kcal/mol (measured barrier %.2f)" % (0.5 * 45 * 23.0605 * 0.275**2, V[(z > -1.9) & (z < -1.0)].max()))
for m in (4.0, 1.0, 0.5, 0.1):
    t = HB2_2M / m / h**2
    w, v = eigh_tridiagonal(2 * t + V, -t * np.ones(len(z) - 1), select="i", select_range=(0, 59))
    p = v**2 / h; dens = p; E0 = w[0]
    outer = np.array([(p[:, k][(z < -top) | (z > top)].sum() * h) for k in range(len(w))])
    left = np.array([(p[:, k][z < -top].sum() * h) for k in range(len(w))])
    pit = [k for k in range(len(w)) if outer[k] > 0.5][:2]
    wt = np.exp(-(w - E0) / KT); rho = (p * wt).sum(1) / wt.sum()
    zone = lambda lo, hi: float(rho[(z >= lo) & (z < hi)].sum() * h)
    d = dict(E0=float(E0), levels=(w[:6] - E0).tolist(), n_states_60=int(len(w)), highest_level_kT=float((w[-1] - E0) / KT),
             frac_left=zone(-9, -top), frac_mid=zone(-top, top), frac_right=zone(top, 9))
    if pit:
        d["pit_states_E_above_ground"] = [float(w[k] - E0) for k in pit]
        if len(pit) == 2: d["pit_splitting_kcal"] = float(w[pit[1]] - w[pit[0]])
        # beyond-the-top probability of the (symmetric) pit ground state: weight of the left-localised combination on z > -top
        k = pit[0]; d["pit_state_weight_|z|>top"] = float(outer[k])
        a, b = p[:, pit[0]], p[:, pit[1]] if len(pit) > 1 else None
    res[m] = d
    print("m=%.1f amu: levels above ground %s ; zones L/mid/R = %.4f/%.4f/%.4f ; pit states %s ; splitting %s ; highest level %.1f kT" % (m, np.round(w[:6] - E0, 2).tolist(), d["frac_left"], d["frac_mid"], d["frac_right"],
          np.round(d.get("pit_states_E_above_ground", []), 2).tolist(), "%.2e" % d["pit_splitting_kcal"] if "pit_splitting_kcal" in d else "-", d["highest_level_kT"]))
wc = np.exp(-V / KT); rc = wc / wc.sum()
zc = lambda lo, hi: float(rc[(z >= lo) & (z < hi)].sum()); res["classical"] = dict(frac_left=zc(-9, -top), frac_mid=zc(-top, top), frac_right=zc(top, 9))
print("classical 1-D Boltzmann: L/mid/R = %.4f/%.4f/%.4f" % (res["classical"]["frac_left"], res["classical"]["frac_mid"], res["classical"]["frac_right"]))
json.dump(res, open(sys.argv[3], "w"))
