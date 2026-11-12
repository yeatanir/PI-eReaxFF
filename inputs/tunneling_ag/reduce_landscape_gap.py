#!/usr/bin/env python3
"""Reduce landscape_gap.py output: pit/mid/top energies, curvatures, T_c per mass, off-axis minima, relaxation effect.
Usage: reduce_landscape_gap.py <landscape.json> <png_out>"""
import json, sys
import numpy as np
KCAL_EV = 0.0433641153; AMU = 103.6427; HBAR = 0.6582119569; KB = 8.617333262e-5; T0 = 300.0
R = json.load(open(sys.argv[1])); z = np.array(R["A_z"]); E = np.array(R["A_E"]); zl, zr = R["zl"], R["zr"]
def at(zq): return float(np.interp(zq, z, E))
ipit_l = np.argmin(np.where(np.abs(z - zl) < 0.3, E, 1e9)); ipit_r = np.argmin(np.where(np.abs(z - zr) < 0.3, E, 1e9))
Epit = min(E[ipit_l], E[ipit_r]); print("pits: left z=%.3f E=%.3f ; right z=%.3f E=%.3f (kcal/mol); |diff| = %.4f" % (z[ipit_l], E[ipit_l], z[ipit_r], E[ipit_r], abs(E[ipit_l] - E[ipit_r])))
imid = np.argmin(np.abs(z)); mid_seg = (z > R["top_z"]["L"] + 0.1) & (z < R["top_z"]["R"] - 0.1); im = np.where(mid_seg)[0][np.argmin(E[mid_seg])]
print("middle region: E(z=0) - E(pit) = %+.3f ; lowest point of the middle region z=%.3f, E - E(pit) = %+.3f ; highest point of the middle region z=%.3f, %+.3f" % (
    E[imid] - Epit, z[im], E[im] - Epit, z[mid_seg][np.argmax(E[mid_seg])], E[mid_seg].max() - Epit))
for s in ("L", "R"):
    b = R["B_" + s]; zb = np.array(b["z"]); Eb = np.array(b["E"]); i = np.argmax(np.where((zb > (zl + 0.1 if s == "L" else zr - 0.7)) & (zb < (zl + 0.7 if s == "L" else zr - 0.1)), Eb, -1e9))
    print("barrier top %s: z = %.3f A (%.3f A from the pit), height above pit = %.3f kcal/mol (%.1f kT at 300 K), above middle-region minimum = %.3f" % (
        s, zb[i], abs(zb[i] - (zl if s == "L" else zr)), Eb[i] - Epit, (Eb[i] - Epit) / (KB * T0 / KCAL_EV), Eb[i] - E[im]))
u = np.array(R["C_disp"]); print("curvature at the top (quadratic fit |u|<=0.08 A), k in eV/A^2 (negative = unstable):")
for s in ("L", "R"):
    ks = {}
    for ax in "zxy":
        c = np.polyfit(u, np.array(R["C"][s + "_" + ax]), 2); ks[ax] = 2 * c[0] * KCAL_EV
    print("  %s: k_z = %+.2f, k_x = %+.2f, k_y = %+.2f" % (s, ks["z"], ks["x"], ks["y"]))
    if ks["z"] < 0:
        for m in (4.0, 1.0, 0.5, 0.1):
            wb = np.sqrt(-ks["z"] / (m * AMU)); Tc = HBAR * wb / (2 * np.pi * KB)
            print("     m = %.1f amu: omega_b = %.3f /fs (%.0f cm^-1), T_c = %.0f K (T_c/300K = %.2f)" % (m, wb, wb * 5308.8, Tc, Tc / 300))
print("off-axis: for each z, lowest energy over the (x,y) grid (+-0.6 A, step 0.15) minus the on-axis value; location of the minimum")
g = np.array(R["D_grid"])
for zk, M in R["D"].items():
    M = np.array(M); i, j = np.unravel_index(np.argmin(M), M.shape); onaxis = M[len(g) // 2, len(g) // 2]
    print("  z = %6.3f: on-axis %+.2f vs pit ; grid min %+.2f (x offset %+.2f, y offset %+.2f) -> gain %.2f kcal/mol" % (float(zk), onaxis - Epit, M[i, j] - Epit, g[i], g[j], onaxis - M[i, j]))
print("relaxation of the other four electrons (12 dof, L-BFGS):")
for k, v in R["E"].items(): print("  %-3s z=%6.3f: E_fixed %.4f -> E_relaxed %.4f (lowering %.4f kcal/mol), max move %.3f A, %d iterations" % (k, v["z"], v["E_fixed"], v["E_relaxed"], v["E_fixed"] - v["E_relaxed"], v["max_move_A"], v["nit"]))
if all(k in R["E"] for k in ("pit", "top", "mid")):
    fx = {k: R["E"][k]["E_fixed"] for k in R["E"]}; rl = {k: R["E"][k]["E_relaxed"] for k in R["E"]}
    print("  barrier height (top - pit): fixed %.3f, relaxed %.3f ; middle - pit: fixed %.3f, relaxed %.3f" % (fx["top"] - fx["pit"], rl["top"] - rl["pit"], fx["mid"] - fx["pit"], rl["mid"] - rl["pit"]))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(z, E - Epit); ax[0].set_xlabel("z (A)"); ax[0].set_ylabel("E - E(pit) (kcal/mol)"); ax[0].set_title("on-axis profile, gap %.1f A" % R["gap"]); ax[0].set_ylim(-1, 30)
for s in ("L", "R"): ax[0].axvline(R["top_z"][s], color="r", lw=0.5)
zk = sorted(R["D"], key=float); M = np.array(R["D"][zk[3]]); im_ = ax[1].imshow((M - Epit).T, origin="lower", extent=[g[0], g[-1], g[0], g[-1]], vmax=25); ax[1].set_xlabel("x offset (A)"); ax[1].set_ylabel("y offset (A)"); ax[1].set_title("E - E(pit), plane z = %s" % zk[3]); fig.colorbar(im_)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=130)
