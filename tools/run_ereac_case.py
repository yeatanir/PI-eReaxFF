#!/usr/bin/env python3
"""Prepare a scratch run directory for the rebuilt ereac from a reference case (control, ffield, geo) and run it via run_ereac.py.

Usage: run_ereac_case.py <ereac_binary> <case_dir_with_control_ffield_geo> <run_dir> [--expect FILE ...]
Copies control, ffield (-> ffield and fort.4), geo (-> geo and fort.3), fort.20 and fort.35 (from the 225-atom Jamil np_export if missing).
Never deletes anything.
"""
import pathlib, shutil, subprocess, sys

here = pathlib.Path(__file__).resolve().parent
binary, case, run = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
expect = sys.argv[sys.argv.index("--expect") + 1:] if "--expect" in sys.argv else ["fort.73", "ereac_export.out"]
run.mkdir(parents=True, exist_ok=True)
ref = here.parent / "reference_data/phase0_jamil_225atom/np_export"
for src, dst in (("control", "control"), ("ffield", "ffield"), ("ffield", "fort.4"), ("geo", "geo"), ("geo", "fort.3")):
    shutil.copy(case / src, run / dst)
for f in ("fort.20", "fort.35"):
    if (case / f).exists():
        shutil.copy(case / f, run / f)
    elif not (run / f).exists():
        (run / f).write_text({"fort.20": "  0\n  0   0: Normal run\n", "fort.35": "0\n"}[f])
sys.exit(subprocess.call([sys.executable, str(here / "run_ereac.py"), binary, str(run), "--expect", *expect]))
