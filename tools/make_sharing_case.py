#!/usr/bin/env python3
"""Purpose-made eReaxFF geometries with electrons close enough to carbons for the sharing term to matter (exp(-vincr r^2), vincr(C) = 181 -> r ~ 0.1 A).

Usage: make_sharing_case.py <out_dir> <tag> --pxel2 <value> --dist d1 d2 ... [--seed N]
Writes control, ffield (C atom pxel2 replaced), fort.4 (= ffield), geo, fort.3 (= geo), fort.20, fort.35 into <out_dir>/<tag>, non-periodic, fixed charges
(El -1 each, the carbons share +n_el), to be run with tools/run_ereac.py (export binary) and compared with test_ereaxff_reference.
"""
import argparse
import math
import pathlib
import random
import shutil

REF = pathlib.Path(__file__).resolve().parents[1] / "reference_data/phase0_jamil_225atom/np_export"
SCRATCH = pathlib.Path(__file__).resolve().parents[2] / "scratch/phase0/jamil_np_stage"


def hetatm(i, name, x, y, z):
    return "%-6s %5d %-5s %-3s %-1s %-5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f" % ("HETATM", i, name, "", "", name, x, y, z, "0", 0, 0, 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("tag")
    ap.add_argument("--pxel2", type=float, default=0.8)
    ap.add_argument("--plp2", type=float, default=None, help="C lone-pair strength (atom line 3, column 2); default: as in the ffield")
    ap.add_argument("--vincr", type=float, default=None, help="C sharing width (atom line 3, column 3)")
    ap.add_argument("--dist", type=float, nargs="+", required=True)
    ap.add_argument("--seed", type=int, default=3)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    d = pathlib.Path(a.out) / a.tag
    d.mkdir(parents=True, exist_ok=True)
    # carbons: ring of 6 plus two attached atoms, slightly distorted
    c = []
    for k in range(6):
        ang = 2 * math.pi * k / 6
        c.append((1.42 * math.cos(ang) + 0.06 * rng.uniform(-1, 1), 1.42 * math.sin(ang) + 0.06 * rng.uniform(-1, 1), 10.0 + 0.06 * rng.uniform(-1, 1)))
    c.append((2.85 + 0.05 * rng.uniform(-1, 1), 0.05 * rng.uniform(-1, 1), 10.0 + 0.05 * rng.uniform(-1, 1)))
    c.append((-2.75 + 0.05 * rng.uniform(-1, 1), 0.95 + 0.05 * rng.uniform(-1, 1), 10.0 + 0.05 * rng.uniform(-1, 1)))
    els = []
    for k, dist in enumerate(a.dist):
        host = c[k % len(c)]
        v = [rng.uniform(-1, 1) for _ in range(3)]
        n = math.sqrt(sum(t * t for t in v))
        els.append(tuple(host[i] + dist * v[i] / n for i in range(3)))
    nel = len(els)
    lines = ["BIOGRF 200", "DESCRP sharing_" + a.tag, "RUTYPE NORMAL RUN", "REMARK purpose-made sharing check", "REMARK fixed-charge MOLCHARGE branch (icharg=8)",
             "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"]
    for i, p in enumerate(els, 1):
        lines.append(hetatm(i, "El", *p))
    for i, p in enumerate(c, nel + 1):
        lines.append(hetatm(i, "C", *p))
    for i in range(1, nel + 1):
        lines.append("MOLCHARGE %4d %4d  -1.00" % (i, i))
    lines.append("MOLCHARGE %4d %4d %7.2f" % (nel + 1, nel + len(c), float(nel)))
    lines.append("END")
    (d / "geo").write_text("\n".join(lines) + "\n")
    shutil.copy(d / "geo", d / "fort.3")
    # ffield with C pxel2 replaced (atom line 2, column 5)
    ff = (REF / "ffield").read_text().split("\n")
    i = next(k for k, l in enumerate(ff) if l.startswith(" C ") and "1.3417" in l)
    cols = ff[i + 1].split()
    cols[4] = "%.4f" % a.pxel2
    ff[i + 1] = "   " + "".join("%9.4f" % float(t) for t in cols)   # legacy format (3x,10f9.4)
    c3 = ff[i + 2].split()
    if a.plp2 is not None: c3[1] = "%.4f" % a.plp2
    if a.vincr is not None: c3[2] = "%.4f" % a.vincr
    ff[i + 2] = "   " + "".join("%9.4f" % float(t) for t in c3)
    (d / "ffield").write_text("\n".join(ff))
    shutil.copy(d / "ffield", d / "fort.4")
    shutil.copy(REF / "control", d / "control")
    for f in ("fort.20", "fort.35"):
        shutil.copy(SCRATCH / f, d / f)
    print("wrote", d)


if __name__ == "__main__":
    main()
