#!/usr/bin/env python3
"""Stage 2 design solver: 3-D Schroedinger problem (schrod3d.py grid/solver) of ONE quantum electron in the frozen eReaxFF landscape PLUS a z-only Gaussian
    Vb(z) = V0 exp(-z^2 / (2 sigma^2))        (same form as ipi/pes/gaussian_barrier.py with s = z, r_mid = 0, acting on the quantum electron only)
Per mass (300 K) it reports:
  E0                    3-D ground-state energy (kcal/mol, V shifted so min of the bare landscape is 0)
  W0 / margin_kT        adiabatic barrier top W(0) = lowest eigenvalue of the 2-D transverse problem at z = 0 (this includes the transverse zero-point energy);
                        margin = (W(0) - E0)/kT  (how far the ground state lies below the effective barrier top)
  rho_q0 / rho_c0       thermal z-marginal density at z = 0, quantum (all states up to emax kT above E0) and classical Boltzmann (same box), and their ratio
  rho_q0 / rho_qmax     quantum density at z = 0 relative to its maximum (is the barrier region populated enough to be measured?)
  Bcl_PMF_kT            classical barrier = ln( min(max_{z<0} rho_c, max_{z>0} rho_c) / rho_c(0) ) (potential of mean force incl. transverse entropy); V0/kT also printed
  split_kcal            tunnel splitting = lowest odd-z state minus E0 inside the (even x, even y) sector
  f_ground / f_thermal  fraction of the quantum density where V > E (ground state only / thermal average with each state's own energy)
Usage: schrod3d_barrier.py <h> <emax_kT> <out.json> <V0:sigma> [<V0:sigma> ...]   (V0=0 gives the bare landscape; use <V0:sigma>=15:0.8 for the earlier 'middle lift')
Masses 4, 1, 0.5, 0.1 amu. Pass --masses to change.
"""
import argparse, json, os, sys, time
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schrod3d import potential, lap1, lap1z, HB2_2M, KT

def solve_sector(h, m, V, sx, sy, emax, nk0=24):
    """Eigenpairs of one (x,y)-parity sector up to emax above its own lowest level (doubling nk)."""
    nx, ny, nz = V.shape; I = lambda n: sp.identity(n, format="csr")
    T = (sp.kron(sp.kron(lap1(nx, sx), I(ny)), I(nz)) + sp.kron(sp.kron(I(nx), lap1(ny, sy)), I(nz)) + sp.kron(sp.kron(I(nx), I(ny)), lap1z(nz)))
    H = ((-HB2_2M / m / h**2) * T + sp.diags(V.ravel(), 0)).tocsc(); nk = nk0
    while True:
        w, v = eigsh(H, k=nk, which="SA", tol=1e-8, ncv=min(max(3 * nk, 60), H.shape[0] - 1))
        if w.max() - w.min() > emax or nk >= 600: break
        nk *= 2
    o = np.argsort(w); return w[o], v[:, o]

def adiabatic_W(h, m, V):
    """W(z_k): lowest eigenvalue of the transverse 2-D (even,even) problem at each z slice."""
    nx, ny, nz = V.shape; I = lambda n: sp.identity(n, format="csr")
    T = sp.kron(lap1(nx, 1), I(ny)) + sp.kron(I(nx), lap1(ny, 1)); W = np.zeros(nz)
    for k in range(nz):
        H = ((-HB2_2M / m / h**2) * T + sp.diags(V[:, :, k].ravel(), 0)).tocsc()
        W[k] = eigsh(H, k=1, which="SA", tol=1e-9)[0][0]
    return W

def run(h, m, V0, sig, emax_kt, lx=1.2, lz=3.0):
    xs, zs, E = potential(h, lx, lz); V = E - E.min()
    V = V + V0 * np.exp(-zs**2 / (2 * sig**2))[None, None, :]
    nx, nz = len(xs), len(zs); emax = emax_kt * KT; dz = h; t0 = time.time()
    mid = nz // 2                                     # zs[mid-1], zs[mid] = -h/2, +h/2: z = 0 sits between them
    st = []                                           # (E, sx, sy, p_z, parity, f_V>E)
    for sx in (1, -1):
        for sy in (1, -1):
            w, v = solve_sector(h, m, V, sx, sy, emax if (sx, sy) != (1, 1) else emax + 2.0 * KT)
            for k in range(len(w)):
                p = v[:, k].reshape(nx, nx, nz) ** 2; p /= p.sum(); pz = p.sum((0, 1))
                par = float(np.sum(v[:, k].reshape(nx, nx, nz) * v[:, k].reshape(nx, nx, nz)[:, :, ::-1]) / np.sum(v[:, k] ** 2))
                st.append((w[k], sx, sy, pz, par, float(p[V > w[k]].sum())))
    E0 = min(s[0] for s in st); sel = [s for s in st if s[0] - E0 < emax]
    wt = np.array([np.exp(-(s[0] - E0) / KT) for s in sel]); rho_q = sum(w_ * s[3] for w_, s in zip(wt, sel)) / wt.sum()
    wc = np.exp(-V / KT); rho_c = wc.sum((0, 1)); rho_c /= rho_c.sum()
    r0 = lambda r: 0.5 * (r[mid - 1] + r[mid]) / dz                            # density per A at z = 0
    # classical PMF barrier
    cl_left, cl_right = rho_c[:mid].max(), rho_c[mid:].max(); Bcl = float(np.log(min(cl_left, cl_right) / (0.5 * (rho_c[mid - 1] + rho_c[mid]))))
    W = adiabatic_W(h, m, V); W0 = 0.5 * (W[mid - 1] + W[mid])
    g0 = [s for s in sel if abs(s[0] - E0) < 1e-9][0]
    odd = [s for s in st if s[1] == 1 and s[2] == 1 and s[4] < 0]; even = [s for s in st if s[1] == 1 and s[2] == 1 and s[4] > 0]
    split = float(min(s[0] for s in odd) - min(s[0] for s in even)) if odd else None
    fth = float(sum(w_ * s[5] for w_, s in zip(wt, sel)) / wt.sum())
    return dict(m=m, V0=V0, sigma=sig, h=h, E0=float(E0), W0=float(W0), margin_kT=float((W0 - E0) / KT), V0_kT=V0 / KT, Bcl_PMF_kT=Bcl,
                rho_q0=float(r0(rho_q)), rho_c0=float(r0(rho_c)), ratio_q_over_c=float(r0(rho_q) / r0(rho_c)), rho_q0_over_qmax=float(r0(rho_q) / (rho_q.max() / dz)),
                split_kcal=split, f_ground=float(g0[5]), f_thermal=fth, n_states=len(sel), time_s=time.time() - t0,
                rms_z=float(np.sqrt((rho_q * zs**2).sum())), zs=zs.tolist(), rho_q=(rho_q / dz).tolist(), rho_c=(rho_c / dz).tolist(), W=W.tolist(),
                levels=[float(x - E0) for x in sorted(s[0] for s in st)[:8]])

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("h", type=float); ap.add_argument("emax_kt", type=float); ap.add_argument("out"); ap.add_argument("designs", nargs="+")
    ap.add_argument("--masses", default="4,1,0.5,0.1"); a = ap.parse_args(); res = []
    for d in a.designs:
        V0, sig = map(float, d.split(":"))
        for m in map(float, a.masses.split(",")):
            r = run(a.h, m, V0, sig, a.emax_kt); res.append(r)
            print("V0=%.1f sigma=%.2f m=%.1f: E0=%.3f W0=%.3f margin=%.2f kT | Bcl_PMF=%.2f kT (V0/kT %.2f) | rho_q0/rho_c0=%.2f rho_q0/qmax=%.3f | split=%s | f_V>E: gs %.3f th %.3f | %d st %.0f s" % (
                V0, sig, m, r["E0"], r["W0"], r["margin_kT"], r["Bcl_PMF_kT"], r["V0_kT"], r["ratio_q_over_c"], r["rho_q0_over_qmax"],
                "%.3f" % r["split_kcal"] if r["split_kcal"] is not None else "-", r["f_ground"], r["f_thermal"], r["n_states"], r["time_s"]), flush=True)
            json.dump(res, open(a.out, "w"))
if __name__ == "__main__":
    main()
