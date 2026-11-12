#!/usr/bin/env python3
"""A1: forces/energies returned by the P LAMMPS clients (client*/fd.dump, lmp.log) against (a) the same LAMMPS plugin evaluated offline
and (b) the exact rebuilt legacy ereac, at 10 bead geometries picked from the run.
ereac reads positions as f10.5, so (a-rounded) and (b) use positions rounded to 5 decimals; the client force at the UNrounded position is compared too.
usage: a1_compare.py <run_dir> [--n 10]
"""
import argparse, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve(); PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ.parent / "forks/metal-reax-lammps/tests/lammps"))
import run_ereaxff as R
LMP = PROJ.parent / "external/install_lammps/bin/lmp"
PLUGIN = PROJ.parent / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"
EREAC = PROJ.parent / "scratch/phase0/ereac_export_build/ereac"
CASE = PROJ / "reference_data/phase4_frames_ag/Ag7_1_md"
TEMPLATE = PROJ.parent / "scratch/phase5_ag2/ipi_ag7"   # fort.20, fort.35

def read_dump(path):
    L = Path(path).read_text().splitlines(); out = []; i = 0
    while i < len(L):
        if L[i].startswith("ITEM: TIMESTEP"):
            step = int(L[i + 1]); n = int(L[i + 3]); k = i + 9
            out.append((step, np.array([[float(x) for x in L[k + j].split()[1:7]] for j in range(n)]))); i = k + n
        else: i += 1
    return out

def read_pe(path):
    rows, on = [], False
    for l in Path(path).read_text().splitlines():
        if l.lstrip().startswith("Step"): on = True; continue
        if l.startswith("Loop time") or l.startswith("ERROR"): on = False
        if on and len(l.split()) == 2:
            try: rows.append(float(l.split()[1]))
            except ValueError: pass
    return rows

def lammps_point(pos, sym, q, masses, tag):
    d = Path(tempfile.mkdtemp(prefix="a1_", dir=PROJ.parent / "scratch/tun_A1")); shutil.copy(CASE / "ffield", d / "ffield")
    R.write_data(d / "data.lmp", sym, pos, q, None, False, masses)
    plug = os.path.relpath(PLUGIN, d)
    pe, pv, f = R.run_lammps(str(LMP), plug, "ffield", d, False, 5.7, R.elements_in_order(sym))
    return pe, np.array(f)

def ereac_point(pos_rounded, geo_lines):
    d = Path(tempfile.mkdtemp(prefix="a1e_", dir=PROJ.parent / "scratch/tun_A1")); out = []; i = 0
    for l in geo_lines:
        if l.startswith("HETATM"):
            l = l[:30] + "".join("%10.5f" % v for v in pos_rounded[i]) + l[60:]; i += 1
        out.append(l)
    (d / "geo").write_text("\n".join(out)); (d / "fort.3").write_text("\n".join(out))
    ff = (CASE / "ffield").read_text(); (d / "ffield").write_text(ff); (d / "fort.4").write_text(ff)
    for f in ("fort.20", "fort.35"): shutil.copy(TEMPLATE / f, d / f)
    c = (CASE / "control").read_text()
    for k, v in (("imdmet", "3"), ("iravel", "0"), ("nmdit", "1"), ("iout1", "1"), ("imetho", "0")):
        p = re.compile(r"^\s*\S+(\s+%s\b.*)$" % k, re.M)
        c = p.sub(lambda m: "%7s%s" % (v, m.group(1)), c, count=1) if p.search(c) else c.rstrip("\n") + "\n%7s %s\n" % (v, k)
    (d / "control").write_text(c)
    for att in range(80):
        (d / "ereac_export.out").unlink(missing_ok=True)
        r = subprocess.run([str(EREAC)], cwd=d, capture_output=True, text=True, timeout=300)
        if "dyld" in r.stderr: time.sleep(0.3); continue
        break
    e, g, st = None, [], None
    for l in (d / "ereac_export.out").read_text().splitlines():
        w = l.split()
        if w[0] == "STEP": st = int(w[1])
        elif st == 0 and w[0] == "ESTRC": e = float(w[1])
        elif st == 0 and w[0] == "GRAD": g.append([float(x) for x in w[2:5]])
    return e, -np.array(g)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--n", type=int, default=10); a = ap.parse_args()
    run = Path(a.run).resolve()
    sym, _, q, _ = R.read_geo(CASE / "geo"); masses = R.ffield_masses(str(CASE / "ffield"))
    geo_lines = (CASE / "geo").read_text().splitlines()
    clients = sorted(run.glob("client*")); data = []
    for c in clients:
        dump = read_dump(c / "fd.dump"); pe = read_pe(c / "lmp.log")
        m = min(len(dump), len(pe))
        for j in range(m):   # fix ipi repeats each configuration (2 of 3 LAMMPS steps are repeats at identical x): keep only new geometries
            if j == 0 or np.abs(dump[j][1][:, :3] - dump[j - 1][1][:, :3]).max() > 0: data.append((c.name, dump[j][0], dump[j][1], pe[j]))
    print("clients:", len(clients), " evaluations dumped:", len(data), " (clients that received at least one request: %d)" % len({d[0] for d in data}))
    rng = np.random.default_rng(3)
    # 10 geometries spread over clients and time; skip nothing
    by_client = {}
    for d in data: by_client.setdefault(d[0], []).append(d)
    names = sorted(by_client); pick = []
    for j in range(a.n):
        c = names[j % len(names)]; lst = by_client[c]; pick.append(lst[int(rng.integers(0, len(lst)))])
    print(f"{'client':8s} {'step':>5s} | {'|F_client - F_lmp(unrounded)|':>30s} | {'|F_lmp - F_ereac| (rounded)':>28s} {'dE(lmp-ereac)':>14s} | {'|F_client - F_ereac|':>20s} {'|F|max':>8s}")
    worst = np.zeros(4)
    for (c, step, arr, pe) in pick:
        x = arr[:, 0:3]; f = arr[:, 3:6]
        pe_off, f_off = lammps_point(x, sym, q, masses, "u")
        xr = np.round(x, 5)
        pe_r, f_r = lammps_point(xr, sym, q, masses, "r")
        e_er, f_er = ereac_point(xr, geo_lines)
        fmax = np.abs(f_er).max()
        d1 = np.abs(f - f_off).max(); d2 = np.abs(f_r - f_er).max(); de = pe_r - e_er; d3 = np.abs(f - f_er).max()
        print(f"{c:8s} {step:5d} | {d1:30.3e} | {d2:28.3e} {de:14.3e} | {d3:20.3e} {fmax:8.3f}")
        worst = np.maximum(worst, [d1 / fmax, d2 / fmax, abs(de) / max(1, abs(e_er)), d3 / fmax])
    print("worst relative to max|F| (or |E|): client-vs-lmp-offline %.2e | lmp(rounded)-vs-ereac %.2e | dE/E %.2e | client(unrounded)-vs-ereac %.2e" % tuple(worst))
main()
