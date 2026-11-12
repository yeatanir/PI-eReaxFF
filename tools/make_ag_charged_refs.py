#!/usr/bin/env python3
"""Corrected Ag/El references with the intended fixed charges (Ag+ +1, El -1).
The delivered phase4_frames_ag geo files put END before the MOLCHARGE lines, so the legacy reader sees no charges (Coulomb = 0 in every frame).
This tool rewrites each selected frame as a geo with MOLCHARGE before END, runs the rebuilt ereac (export build) for one step, and writes frames_charged.json
(same layout as frames.json; positions are those quantised to the geo format 3f10.5, as the legacy code reads them).
usage: make_ag_charged_refs.py <case_dir> <ereac_export_binary> <work_dir> <frame_index>... """
import json, shutil, subprocess, sys, time
from pathlib import Path

case, binary, work = Path(sys.argv[1]), Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve()
idx = [int(x) for x in sys.argv[4:]]
frames = json.loads((case / "frames.json").read_text())["frames"]
geo = (case / "geo").read_text().split("\n")
hetatm = [i for i, l in enumerate(geo) if l.startswith("HETATM")]
mol = [l for l in geo if l.startswith("MOLCHARGE")]
head = [l for l in geo if not l.startswith(("HETATM", "MOLCHARGE", "END"))]
out = {"terms_order": json.loads((case / "frames.json").read_text())["terms_order"], "frames": []}
ctrl = (case / "control").read_text()
import re
ctrl = re.sub(r"^\s*\d+ nmdit", "      1 nmdit", ctrl, flags=re.M)
for k in idx:
    fr = frames[k]; d = work / f"{case.name}_f{k}"; d.mkdir(parents=True, exist_ok=True)
    lines = list(head)
    qpos = []
    for n, i in enumerate(hetatm):
        p = fr["pos"][n]
        l = geo[i][:30] + "".join("%10.5f" % x for x in p) + geo[i][60:]
        lines.append(l); qpos.append([float(l[30:40]), float(l[40:50]), float(l[50:60])])
    lines += mol + ["END"]
    (d / "geo").write_text("\n".join(lines) + "\n"); shutil.copy(d / "geo", d / "fort.3")
    shutil.copy(case / "ffield", d / "ffield"); shutil.copy(case / "ffield", d / "fort.4")
    (d / "control").write_text(ctrl)
    for f in ("fort.20", "fort.35"): (d / f).write_text("0\n") if f == "fort.35" else (d / f).write_text(" 0\n 0   0: Normal run\n")
    for attempt in range(40):
        for f in ("ereac_export.out",): (d / f).unlink(missing_ok=True)
        r = subprocess.run([str(binary)], cwd=d, capture_output=True, text=True)
        if r.returncode == 0 and (d / "ereac_export.out").exists(): break
        time.sleep(1.0)   # intermittent dyld start-up failure under load
    else: sys.exit("ereac failed: " + r.stderr[-200:])
    terms, est, grad, step = None, None, [], None
    for l in (d / "ereac_export.out").read_text().splitlines():
        w = l.split()
        if w[0] == "STEP": step = int(w[1])
        elif step == 0 and w[0] == "TERMS": terms = [float(x) for x in w[1:15]]
        elif step == 0 and w[0] == "ESTRC": est = float(w[1])
        elif step == 0 and w[0] == "GRAD": grad.append([float(x) for x in w[2:5]])
    out["frames"].append({"step": fr["step"], "pos": qpos, "terms": terms, "estrc": est, "grad": grad})
(case / "frames_charged.json").write_text(json.dumps(out))
print("wrote", case / "frames_charged.json", len(out["frames"]), "frames; coulomb of first:", out["frames"][0]["terms"][10])
