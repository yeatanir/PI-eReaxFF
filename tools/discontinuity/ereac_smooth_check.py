#!/usr/bin/env python3
"""Validate a control-file mode (legacy|smooth) of the adapter against the rebuilt exact ereac run with the SAME thresholds on its control file
(ereac keys: vlbora = bond-order distance range, default 5.0; cutof2 = 3-/4-body BO threshold, default 0.001). usage: ereac_smooth_check.py <mode> <fc_dir> <work_dir> [stride]"""
import os, sys, json, re, shutil, subprocess, time
from pathlib import Path
import numpy as np
B = Path(os.environ.get("EB_LINK", "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build")); P = B / "project"   # EB_LINK: space-free symlink to the tree (LAMMPS cannot load plugins from paths with spaces)
mode, fc, work = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]).resolve(); stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1
binary = (B / "scratch/phase0/ereac_export_build/ereac").resolve(); work.mkdir(parents=True, exist_ok=True)
extra = {"legacy": [], "smooth": ["    5.2 vlbora", "1.0E-09 cutof2"]}[mode]
Z = np.load(fc / "frames_compare.npz")
names = ["bond", "over+under", "lone_pair", "molecule", "valence+penalty", "angle_conj", "hbond", "torsion", "torsion_conj", "vdw", "coulomb", "charge", "field", "electron"]
worst = dict(E=0, F=0, T=0, tors=0); nfr = 0
for cd in sorted((P / "reference_data/phase4_frames_ag").glob("Ag*_md")):
    frames = json.loads((cd / "frames.json").read_text())["frames"]
    geo = (cd / "geo").read_text().split("\n"); het = [i for i, l in enumerate(geo) if l.startswith("HETATM")]
    mol = [l for l in geo if l.startswith("MOLCHARGE")]; head = [l for l in geo if not l.startswith(("HETATM", "MOLCHARGE", "END"))]
    ctrl = re.sub(r"^\s*\d+ nmdit", "      1 nmdit", (cd / "control").read_text(), flags=re.M).rstrip("\n") + "\n" + "\n".join(extra) + ("\n" if extra else "")
    Ez, Tz, Fz = (Z[f"{cd.name}__{mode}__{n}"] for n in ("E", "T", "F"))
    for k in range(0, len(frames), stride):
        d = work / f"{cd.name}_f{k}"; d.mkdir(exist_ok=True)
        lines = list(head)
        for n, i in enumerate(het): lines.append(geo[i][:30] + "".join("%10.5f" % x for x in frames[k]["pos"][n]) + geo[i][60:])
        lines += mol + ["END"]
        (d / "geo").write_text("\n".join(lines) + "\n"); shutil.copy(d / "geo", d / "fort.3")
        shutil.copy(cd / "ffield", d / "ffield"); shutil.copy(cd / "ffield", d / "fort.4"); (d / "control").write_text(ctrl)
        (d / "fort.20").write_text(" 0\n 0   0: Normal run\n"); (d / "fort.35").write_text("0\n"); (d / "ereac_export.out").unlink(missing_ok=True)
        for attempt in range(60):
            r = subprocess.run([str(binary)], cwd=d, capture_output=True, text=True)
            if r.returncode == 0 and (d / "ereac_export.out").exists(): break
            time.sleep(1.0)
        else: sys.exit("ereac failed " + r.stderr[-200:])
        terms, est, grad, step = None, None, [], None
        for l in (d / "ereac_export.out").read_text().splitlines():
            w = l.split()
            if w[0] == "STEP": step = int(w[1])
            elif step == 0 and w[0] == "TERMS": terms = [float(x) for x in w[1:15]]
            elif step == 0 and w[0] == "ESTRC": est = float(w[1])
            elif step == 0 and w[0] == "GRAD": grad.append([float(x) for x in w[2:5]])
        # positions used by ereac are quantised to 5 decimals: evaluate the adapter on the same positions is not available here, so compare with
        # the adapter on the unquantised frame and accept the quantisation error (~|F|*5e-6 ~ 1e-4 kcal/mol)
        dE = est - Ez[k]; dF = np.abs(-np.array(grad) - Fz[k]).max()
        worst["E"] = max(worst["E"], abs(dE)); worst["F"] = max(worst["F"], dF); nfr += 1
print(f"mode={mode}: {nfr} frames; max |E_ereac - E_adapter| = {worst['E']:.2e} kcal/mol; max |F diff| = {worst['F']:.2e} kcal/mol/A  (positions quantised to 1e-5 A in ereac; legacy-mode level is the baseline)")
