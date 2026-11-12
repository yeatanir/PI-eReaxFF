#!/usr/bin/env python3
"""Exact 3-D Schroedinger problem of ONE quantum electron (mass m) in the frozen eReaxFF landscape of the 4.0 A donor-barrier-acceptor system
(Ag fixed, other four El frozen on their cations). V(x,y,z) comes from LAMMPS+Metal-Reax (fix ipi with a mini server, see probe_server.py);
V is even in x, y, z (mirror symmetry of the Ag3 pair), so only x,y >= 0 is solved: four parity sectors (even/odd in x and y), cell-centred grid,
second-order finite differences (convergence checked by comparing two grid spacings). Thermal density rho(z) = sum_n exp(-beta E_n) |psi_n|^2 / Z.
Usage: schrod3d.py <h> <masses comma> [--lx 1.2 --lz 3.0 --emax-kt 14]
"""
import argparse, json, os, subprocess, sys, time
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HB2_2M = 0.0481741            # hbar^2/(2 * 1 amu) in kcal/mol A^2
KT = 0.0019872041 * 300.0
HERE = os.path.dirname(os.path.abspath(__file__)); STAGE = os.path.abspath(os.path.join(HERE, "../../../scratch/tun_stage"))

def potential(h, lx, lz):
    f = os.path.join(STAGE, "V3d_h%.4f_lx%.2f_lz%.2f.npz" % (h, lx, lz))
    nx = int(round(lx / h)); nz = int(round(2 * lz / h))
    xs = (np.arange(nx) + 0.5) * h; zs = -lz + (np.arange(nz) + 0.5) * h
    if not os.path.exists(f):
        X, Y, Z = np.meshgrid(xs, xs, zs, indexing="ij"); pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
        np.save(os.path.join(STAGE, "pts_tmp.npy"), pts)
        t0 = time.time(); subprocess.run([os.path.join(STAGE, "probe.sh"), os.path.join(STAGE, "pts_tmp.npy"), os.path.join(STAGE, "E_tmp.npy"), "s3d"], check=True)
        E = np.load(os.path.join(STAGE, "E_tmp.npy")).reshape(nx, nx, nz); np.savez(f, xs=xs, zs=zs, E=E); print("potential evaluated", E.size, "points in %.0f s" % (time.time() - t0), flush=True)
    d = np.load(f); return d["xs"], d["zs"], d["E"]

def lap1(n, sgn):                         # 1-D second difference, cell-centred, mirror at index -1/2 (sgn=+1 even, -1 odd), Dirichlet beyond n
    L = sp.diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1], format="lil"); L[0, 0] += sgn
    return L.tocsr()
def lap1z(n):
    return sp.diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1], format="csr")

def solve(h, m, V, nk, sx, sy):
    nx, ny, nz = V.shape; I = lambda n: sp.identity(n, format="csr")
    T = sp.kron(sp.kron(lap1(nx, sx), I(ny)), I(nz)) + sp.kron(sp.kron(I(nx), lap1(ny, sy)), I(nz)) + sp.kron(sp.kron(I(nx), I(ny)), lap1z(nz))
    H = (-HB2_2M / m / h**2) * T + sp.diags(V.ravel(), 0)
    w, v = eigsh(H.tocsc(), k=nk, which="SA", tol=1e-7, ncv=min(max(3 * nk, 60), H.shape[0] - 1))
    return w, v

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("h", type=float); ap.add_argument("masses"); ap.add_argument("--lx", type=float, default=1.2)
    ap.add_argument("--lz", type=float, default=3.0); ap.add_argument("--emax-kt", type=float, default=14.0); ap.add_argument("--nk", type=int, default=24); ap.add_argument("--out", default="")
    ap.add_argument("--lift", type=float, default=0.0, help="design probe: add H*exp(-z^2/(2 s^2)) (kcal/mol) in the middle, H = this value")
    ap.add_argument("--lift-sigma", type=float, default=0.8)
    ap.add_argument("--pit-wall", type=float, default=None, help="also solve the left pit alone with a hard wall at z = this value (A): lowest 3-D pit state energy")
    a = ap.parse_args()
    xs, zs, E = potential(a.h, a.lx, a.lz); V = E - E.min()
    if a.lift: V = V + a.lift * np.exp(-zs**2 / (2 * a.lift_sigma**2))[None, None, :]
    nx, nz = len(xs), len(zs); res = {}
    Vpit = None
    if a.pit_wall is not None:
        Vpit = V.copy(); Vpit[:, :, zs > a.pit_wall] = 1e3
    for m in map(float, a.masses.split(",")):
        t0 = time.time(); levels = []; dens = []; tail = []
        for sx in (1, -1):
            for sy in (1, -1):
                nk = a.nk
                while True:
                    w, v = solve(a.h, m, V, nk, sx, sy)
                    if (w.max() - 0.0) > 0 and w.max() > (levels_min := min([l[0] for l in levels] + [w.min()])) + a.emax_kt * KT or nk >= 400: break
                    nk *= 2
                for k in range(len(w)):
                    p = (v[:, k].reshape(nx, nx, nz) ** 2); levels.append((w[k], sx, sy)); dens.append(p.sum((0, 1)) / p.sum())
                    tail.append(max(p[-3:, :, :].sum(), p[:, -3:, :].sum(), p[:, :, :3].sum()) / p.sum())
        lev = np.array([l[0] for l in levels]); order = np.argsort(lev); E0 = lev[order[0]]
        outer = np.array([d[(zs < -1.7) | (zs > 1.7)].sum() for d in dens])            # weight of each state beyond |z| > 1.7 (the two pit regions)
        pit_states = [(float(lev[i] - E0), float(outer[i])) for i in order if outer[i] > 0.5][:3]
        wt = np.exp(-(lev - E0) / KT); sel = (lev - E0) < a.emax_kt * KT
        rho = sum(wt[i] * dens[i] for i in np.where(sel)[0]) / wt[sel].sum()          # normalised per cell (sum over z = 1) on the half-space (x,y>=0), full density in z
        zone = lambda lo, hi: float(rho[(zs >= lo) & (zs < hi)].sum())
        top = 1.70
        res[m] = dict(E0=float(E0), levels=[float(x) for x in np.sort(lev)[:12] - E0], n_states=int(sel.sum()), lowest_pit_states_E_above_ground_and_weight=pit_states, highest_included_kT=float(((lev[sel] - E0) / KT).max()),
                      frac_left=zone(-9, -top), frac_mid=zone(-top, top), frac_right=zone(top, 9), P_beyond_top_from_left_pit=zone(-top, 9),
                      worst_wavefunction_edge_weight=float(max(tail[i] for i in np.where(sel)[0])), rho_z=rho.tolist(), zs=zs.tolist(), time_s=time.time() - t0)
        if Vpit is not None:
            wp, _ = solve(a.h, m, Vpit, 3, 1, 1); res[m]["pit_ground_E_above_global_ground"] = float(wp[0] - E0)
            print("   left pit alone (hard wall at z=%.2f): ground state %.3f kcal/mol above the global ground state = %.1f kT" % (a.pit_wall, wp[0] - E0, (wp[0] - E0) / KT), flush=True)
        print("m=%.2f amu: E0=%.3f kcal/mol, %d states up to %.1f kT, zones L/mid/R = %.4f / %.4f / %.4f, edge weight %.1e, first levels %s, %.0f s" % (
            m, E0, sel.sum(), res[m]["highest_included_kT"], res[m]["frac_left"], res[m]["frac_mid"], res[m]["frac_right"], res[m]["worst_wavefunction_edge_weight"],
            np.round(res[m]["levels"][:6], 2).tolist(), res[m]["time_s"]), flush=True)
        print("   lowest pit-localised states (E above the ground state in kcal/mol, weight at |z|>1.7):", [(round(e, 2), round(w, 2)) for e, w in pit_states], " P(|z|>1.7) = %.2e" % (res[m]["frac_left"] + res[m]["frac_right"]), flush=True)
    # classical Boltzmann reference on the same grid (full 3-D, includes the transverse volume)
    wc = np.exp(-V / KT); rc = wc.sum((0, 1)); rc /= rc.sum(); zonec = lambda lo, hi: float(rc[(zs >= lo) & (zs < hi)].sum())
    res["classical_boltzmann"] = dict(frac_left=zonec(-9, -1.7), frac_mid=zonec(-1.7, 1.7), frac_right=zonec(1.7, 9), rho_z=rc.tolist())
    print("classical Boltzmann (3-D): L/mid/R = %.4f / %.4f / %.4f" % (res["classical_boltzmann"]["frac_left"], res["classical_boltzmann"]["frac_mid"], res["classical_boltzmann"]["frac_right"]))
    if a.out: json.dump(res, open(a.out, "w"))
if __name__ == "__main__":
    main()
