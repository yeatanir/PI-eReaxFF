#!/usr/bin/env python3
"""A5/A6: wells, barriers, curvatures and frequencies from scan_landscape.py JSON files.
Usage: analyze_landscape.py <png_out> <json> [<json> ...]"""
import json, sys
import numpy as np
KCAL_EV = 0.0433641153; AMU = 103.6427; HBAR = 0.6582119569; KB = 8.617333262e-5; T = 300.0

def minima(z, E, fitw=0.25):
    out = []
    for i in range(1, len(z) - 1):
        if E[i] <= E[i - 1] and E[i] <= E[i + 1]:
            m = np.abs(z - z[i]) <= fitw + 1e-9
            c = np.polyfit(z[m] - z[i], E[m], 2)
            zmin = z[i] - c[1] / (2 * c[0]) if c[0] > 0 else z[i]
            out.append(dict(z=float(zmin), E=float(np.polyval(c, zmin - z[i])), k_kcal=float(2 * c[0]), k_eV=float(2 * c[0] * KCAL_EV)))
    return out

def omega(k_eV, m_amu): return float(np.sqrt(k_eV / (m_amu * AMU)))
def req_P(w): return float(np.ceil((HBAR * w / (KB * T)) / np.sqrt(0.24)))      # (beta hbar w / P)^2 / 24 < 1 %
def rfree(m_amu): return float(np.sqrt(3 * HBAR ** 2 / (12 * m_amu * AMU * KB * T)))   # sqrt(3 d hbar^2/(12 m kT)) free ring, large P (3-D rms about centroid)

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for f in sys.argv[2:]:
        d = json.load(open(f)); z = np.array(d["z"]); E = np.array(d["E_kcal"])
        ms = minima(z, E); Emin = min(m["E"] for m in ms)
        ax.plot(z, E - E.min(), label="%s gap %.1f" % (d["kind"], d["gap"]))
        wells = [m for m in ms]
        s = "%-5s gap %.2f A: " % (d["kind"], d["gap"])
        if len(wells) >= 2:
            zb = np.linspace(wells[0]["z"], wells[-1]["z"], 400); Eb = np.interp(zb, z, E); ib = np.argmax(Eb)
            s += "wells at z=" + ", ".join("%.2f (E-Emin %.2f)" % (m["z"], m["E"] - Emin) for m in wells) + \
                 "; barrier top z=%.2f, height above lower well %.2f kcal/mol (%.1f kT)" % (zb[ib], Eb[ib] - Emin, (Eb[ib] - Emin) / 0.5962)
        else:
            s += "single minimum z=%.2f" % wells[0]["z"]
        print(s)
        for m in wells: print("      well z=%.2f  k=%.1f eV/A^2  omega(1/0.5/0.1 amu)= %.3f/%.3f/%.3f fs^-1" % (m["z"], m["k_eV"], *(omega(m["k_eV"], x) for x in (1, 0.5, 0.1))))
    ax.set_xlabel("transferring electron z (A)"); ax.set_ylabel("E - Emin (kcal/mol)"); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(sys.argv[1], dpi=130)
