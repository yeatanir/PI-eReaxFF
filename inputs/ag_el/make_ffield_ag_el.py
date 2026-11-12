#!/usr/bin/env python3
"""Write ffield_Ag_El (legacy fixed-column layout) from the text of the SI of Evangelisti, Fichthorn, van Duin, J. Chem. Phys. 153, 104106 (2020).
Usage: make_ffield_ag_el.py <si_text_from_pdftotext> <out_ffield>   (reads the numbers, never retyped by hand)"""
import re, sys
L = open(sys.argv[1]).read().split("\n")
i0 = [i for i, l in enumerate(L) if l.startswith("Force Field File")][0]
body = [l for l in L[i0 + 1:] if l.strip()]
h = [i for i, l in enumerate(body) if l.lstrip().startswith("[ ereaxff")][0]
title = body[h].strip()
gen, names = [], []
for l in body[h + 2:h + 2 + 39]:
    m = re.match(r"\s*(-?\d+\.\d+)\s*!(.*)$", l); gen.append(float(m.group(1))); names.append(m.group(2).strip())
rest = body[h + 2 + 39:]
nums = lambda s: [float(x) for x in s.split()]
ia = [i for i, l in enumerate(rest) if "Nr of atoms" in l][0]
atoms, j = [], ia + 4
for a in range(2):
    sym, r0 = rest[j].split()[0], rest[j].split()[1:]
    rows = [[float(x) for x in r0]] + [nums(rest[j + k]) for k in (1, 2, 3)]
    assert [len(r) for r in rows] == [8] * 4; atoms.append((sym, rows)); j += 4
ib = [i for i, l in enumerate(rest) if "Nr of bonds" in l][0]
b1, b2 = rest[ib + 2].split(), rest[ib + 3].split()
assert len(b1) == 10 and len(b2) == 8, (b1, b2)
ia2 = [i for i, l in enumerate(rest) if "Nr of angles" in l][0]; ang = rest[ia2 + 1].split()
it = [i for i, l in enumerate(rest) if "Nr of torsions" in l][0]; tor = rest[it + 1].split()
f = lambda xs: "".join("%9.4f" % float(x) for x in xs)
out = [title.ljust(100), " 39       ! Number of general parameters".ljust(60)]
out += [("%10.4f !%s" % (g, n)).ljust(48) for g, n in zip(gen, names)]
out += [" 2    ! Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;#",
        "            alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.               ",
        "            cov r3;Elp;Heat inc.;n.u.;n.u.;n.u.;n.u.                            ",
        "            ov/un;val1;n.u.;val3,vval4                                          "]
for sym, rows in atoms:
    out.append(" %-2s" % sym + f(rows[0]))
    out += ["   " + f(r) for r in rows[1:]]
out += [" 1      ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6      ",
        "                         pbe2;pbo3;pbo4;Etrip;pbo1;pbo2;ovcorr                  ",
        "%3d%3d" % (int(b1[0]), int(b1[1])) + f(b1[2:]), "      " + f(b2),
        " 0    ! Nr of off-diagonal terms; Ediss;Ro;gamma;rsigma;rpi;rpi2   ",
        " 1    ! Nr of angles;at1;at2;at3;Thetao,o;ka;kb;pv1;pv2            ",
        "%3d%3d%3d" % tuple(int(x) for x in ang[:3]) + f(ang[3:]),
        " 1    ! Nr of torsions;at1;at2;at3;at4;;V1;V2;V3;V2(BO);vconj;n.u;n",
        "%3d%3d%3d%3d" % tuple(int(x) for x in tor[:4]) + f(tor[4:]),
        " 0                                                                 "]
open(sys.argv[2], "w").write("\n".join(out) + "\n")
