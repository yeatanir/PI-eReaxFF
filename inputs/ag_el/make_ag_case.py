#!/usr/bin/env python3
"""Build an ereac case directory for a neutral Ag cluster from the SI of Evangelisti et al., JCP 153, 104106.
Each Ag (type Av) is an Ag+ cation (MOLCHARGE +1) paired with an explicit electron El (MOLCHARGE -1) placed at the
cation position + 0.001 A in a rotating coordinate (SI: 'small (+-0.001 A) offset'). Fixed charges (icharg 8).
Usage: make_ag_case.py <si_text> <cluster_name e.g. Ag7_1> <out_dir> [--charge-ions N]
"""
import pathlib, re, sys
si, name, out = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3])
L = open(si).read().split("\n")
i = [k for k, l in enumerate(L) if l.strip() == name][0]
pos = []
for l in L[i + 1:]:
    w = l.split()
    if len(w) == 4 and w[0] == "Ag": pos.append([float(x) for x in w[1:]])
    elif pos: break
pos = pos if pos else None
n = len(pos)
cen = [sum(p[k] for p in pos) / n for k in range(3)]
pos = [[p[k] - cen[k] + 40.0 for k in range(3)] for p in pos]
el = [[p[0] + (0.001 if j % 3 == 0 else 0), p[1] + (0.001 if j % 3 == 1 else 0), p[2] + (0.001 if j % 3 == 2 else 0)] for j, p in enumerate(pos)]
g = ["BIOGRF 200", "DESCRP " + name + "_El", "RUTYPE NORMAL RUN", "REMARK Ag cluster + explicit electrons (Evangelisti et al. 2020)",
     "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"]
idx = 1
for p in pos:
    g.append("HETATM %5d %-5s %-3s %1s %-5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f" % (idx, "Av", "Av", " ", "Av", *p, "", 0, 0, 0.0)); idx += 1
for p in el:
    g.append("HETATM %5d %-5s %-3s %1s %-5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f" % (idx, "El", "El", " ", "El", *p, "", 0, 0, 0.0)); idx += 1
for k in range(1, n + 1): g.append("MOLCHARGE %4d %4d %6.2f" % (k, k, 1.0))
for k in range(n + 1, 2 * n + 1): g.append("MOLCHARGE %4d %4d %6.2f" % (k, k, -1.0))
g.append("END")  # END must be last: the legacy reader stops at END, so MOLCHARGE lines after it are ignored (all charges 0)
out.mkdir(parents=True, exist_ok=True)
(out / "geo").write_text("\n".join(g) + "\n"); (out / "fort.3").write_text("\n".join(g) + "\n")
here = pathlib.Path(__file__).resolve().parent
ff = (here / "ffield_Ag_El").read_text()
(out / "ffield").write_text(ff); (out / "fort.4").write_text(ff)
ref = here.parents[1] / "reference_data/phase0_jamil_225atom/np_export"
ctl = (ref / "control").read_text()
ctl = re.sub(r"^\s*80\.00000 axis(\d)", lambda m: "   80.00000 axis" + m.group(1), ctl, flags=re.M)
(out / "control").write_text(ctl)
(out / "fort.20").write_text((ref / "../np_export/../np_export/ffield").read_text()[:0] or "  0\n  0   0: Normal run\n")
(out / "fort.35").write_text("0\n")
print(name, n, "Ag +", n, "El written to", out)
