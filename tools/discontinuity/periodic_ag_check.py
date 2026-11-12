#!/usr/bin/env python3
"""Periodic Ag/El frames: rebuilt exact ereac (CRYSTX cell, minimum-image/periodic pair list) vs the LAMMPS adapter, in a cubic cell of edge L.
The cluster is placed so that it STRADDLES the cell boundary (atoms wrapped into [0,L)). usage: periodic_ag_check.py <L> <work_dir> [stride] [mode legacy|smooth]"""
import sys, json, re, shutil, subprocess, os
from pathlib import Path
import numpy as np
# the harness needs space-free paths: use a symlink to the tree (ln -s "<tree>" $EB_LINK)
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build")); P = B / "project"   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
L = float(sys.argv[1]); work = Path(sys.argv[2]); stride = int(sys.argv[3]) if len(sys.argv) > 3 else 5
case_src = P / "reference_data/phase4_frames_ag/Ag7_1_md"
work.mkdir(parents=True, exist_ok=True); cd = work / "case"; cd.mkdir(exist_ok=True)
fj = json.loads((case_src / "frames.json").read_text()); frames = fj["frames"][::stride]
for f in frames:
    a = np.array(f["pos"]); a = a - np.array([40.0, 40.0, 40.0]); f["pos"] = (a % L).tolist()       # cluster centre at the origin: straddles the boundary
geo = (case_src / "geo").read_text().split("\n"); het = [i for i, l in enumerate(geo) if l.startswith("HETATM")]
head = [l for l in geo if not l.startswith(("HETATM", "MOLCHARGE", "END"))]
cry = "CRYSTX   " + "%11.5f" * 6 % (L, L, L, 90.0, 90.0, 90.0)     # legacy reader format (8x,6f11.5)
lines = head + [cry] + [geo[i][:30] + "".join("%10.5f" % x for x in frames[0]["pos"][n]) + geo[i][60:] for n, i in enumerate(het)]
lines += [l for l in geo if l.startswith("MOLCHARGE")] + ["END"]
(cd / "geo").write_text("\n".join(lines) + "\n"); shutil.copy(case_src / "ffield", cd / "ffield")
ctl = re.sub(r"^\s*80\.00000 axis(\d)", lambda m: "%10.5f axis%s" % (L, m.group(1)), (case_src / "control").read_text(), flags=re.M)
(cd / "control").write_text(ctl); json.dump({"terms_order": fj["terms_order"], "frames": frames}, open(cd / "frames.json", "w"))
r = subprocess.run([sys.executable, str(P / "tools/make_ag_charged_refs.py"), str(cd), str(B / "scratch/phase0/ereac_export_build/ereac"), str(work / "ereac"), *map(str, range(len(frames)))], capture_output=True, text=True)
print(r.stdout[-200:], r.stderr[-300:])
env = dict(os.environ, DYLD_LIBRARY_PATH=str(B / "external/install_lammps/lib"))
# tolerances = the periodic gates of the earlier 60C/60El frames
h = subprocess.run([sys.executable, str(B / "forks/metal-reax-lammps/tests/lammps/run_ereaxff.py"), "--lmp", str(B / "external/install_lammps/bin/lmp"),
                    "--plugin", str(B / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"), "--ffield", str(cd / "ffield"), "--geo", str(cd / "geo"),
                    "--frames", str(cd / "frames_charged.json"), "--periodic", "--taper", "5.7", "--tol-term", "1e-7", "--tol-torsion", "1e-6", "--tol-force", "2e-6"], capture_output=True, text=True, env=env)
print(h.stdout[-900:], h.stderr[-600:])
