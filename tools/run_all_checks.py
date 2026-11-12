#!/usr/bin/env python3
"""One command that runs every existing check and prints PASS / FAIL / SKIP per check (exit status 1 if any FAIL).

  python3 tools/run_all_checks.py                # quick set (~10-15 min on 1-2 cores)
  python3 tools/run_all_checks.py --full         # also the full classical-atoms pytest file (~10 min more)
  python3 tools/run_all_checks.py --only fd_forces schrodinger
  python3 tools/run_all_checks.py --update-ipi-reference   # (re)store the stored i-PI+LAMMPS output (regression baseline, NOT an independent reference)

Checks (names for --only):
  frames_ag        Ag/El reference frames: adapter vs exact ereac (stored frames_charged.json), legacy thresholds
  frames_graphene  C/El periodic frames (stored ereac references), legacy thresholds
  fd_forces        finite-difference test: adapter forces vs central differences of its own energy (legacy and smooth control files)
  smooth_vs_ereac  smooth control file vs ereac run with the same thresholds (vlbora 5.2, cutof2 1e-9)
  conservation_p1  Ag7_1 NVE, mobile Ag, 1 ps: smooth span falls as dt^2 (dt 0.25 / 0.125)
  classical_atoms  i-PI classical-atoms tests (quick subset; --full for the whole file)
  ipi_regression   i-PI's own regression runner on the CLASSICAL_ATOMS tests
  ipi_lammps_short short i-PI + LAMMPS run (P=4, NVE) against stored output
  schrodinger      1-D solver re-run vs stored result, and stored PIMD density moments vs solver (3 sigma)
  eparams          Metal-Reax unit test + LAMMPS test for the `eparams` species table (needs --mr-build)
  periodic_ag      periodic Ag/El frames vs ereac, cell > 2 x cutoff
Every LAMMPS check prints the control file it used (CONTROL_USED ...)."""
import argparse, json, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path

HERE = Path(__file__).resolve(); PROJ = HERE.parents[1]; TREE = PROJ.parent
ap = argparse.ArgumentParser()
ap.add_argument("--only", nargs="*"); ap.add_argument("--full", action="store_true"); ap.add_argument("--update-ipi-reference", action="store_true")
ap.add_argument("--mr-build", default="", help="build dir of the metal-reax-lammps checkout that has the eparams change (for the `eparams` check)")
ap.add_argument("--ipi-reg-dir", default=str(TREE / "scratch/ipi_reg_diag"), help="i-PI checkout used for the regression runner (needs the committed ref_simulation.out files)")
a = ap.parse_args()

TMP = Path(tempfile.mkdtemp(prefix="allchecks_"))
LINK = TMP / "eb"; os.symlink(TREE, LINK)   # the harnesses need space-free absolute paths
LMP = LINK / "external/install_lammps/bin/lmp"; PLUG = LINK / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"
LT = LINK / "forks/metal-reax-lammps/tests/lammps"; P = LINK / "project"
ENV = dict(os.environ, EB_LINK=str(LINK), OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(LINK / "external/install_lammps/lib"))
sys.path.insert(0, str(LT)); import run_ereaxff as R   # noqa: E402
results = []

def sh(cmd, cwd=None, timeout=3600, env=None):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env or ENV, capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout + r.stderr

def check(name):
    def deco(f):
        f.check_name = name; CHECKS.append(f); return f
    return deco
CHECKS = []

@check("frames_ag")
def frames_ag():
    out = []
    for case in ("Ag3_md", "Ag7_1_md"):
        d = P / "reference_data/phase4_frames_ag" / case; w = TMP / ("ag_" + case); shutil.copytree(d, w)
        rc, o = sh([sys.executable, P / "tools/make_ag_charged_refs.py", w, TREE / "scratch/phase0/ereac_export_build/ereac", TMP / ("agw_" + case), 0, 5, 10, 19], timeout=1800)
        if rc: out.append((False, case + ": make_ag_charged_refs failed " + o[-150:])); continue
        rc, o = sh([sys.executable, LT / "run_ereaxff.py", "--lmp", LMP, "--plugin", PLUG, "--ffield", w / "ffield", "--geo", w / "geo", "--frames", w / "frames_charged.json", "--taper", "5.7"])
        out.append((rc == 0 and "RESULT: PASS" in o, case + ": " + (re.findall(r"WORST.*", o) or [o[-150:]])[0]))
    return all(x[0] for x in out), "CONTROL_USED NULL(legacy); " + "; ".join(x[1] for x in out)

@check("frames_graphene")
def frames_graphene():
    d = P / "reference_data/phase4_frames"
    rc, o = sh([sys.executable, LT / "run_ereaxff.py", "--lmp", LMP, "--plugin", PLUG, "--ffield", d / "ffield", "--geo", d / "geo", "--frames", d / "frames.json", "--periodic", "--taper", "10", "--max-frames", "8", "--tol-term", "1e-7", "--tol-torsion", "1e-6", "--tol-force", "2e-6"])
    return rc == 0 and "RESULT: PASS" in o, "CONTROL_USED NULL(legacy); " + (re.findall(r"WORST.*", o) or [o[-200:]])[0]

@check("fd_forces")
def fd_forces():
    import numpy as np
    case = P / "reference_data/phase4_frames_ag/Ag7_1_md"; sym, _, _, _ = R.read_geo(case / "geo"); q = [1.0 if s == "Av" else -1.0 for s in sym]
    pos = np.array(json.loads((case / "frames.json").read_text())["frames"][3]["pos"]); n = len(sym); h = 2e-4
    masses = R.ffield_masses(str(case / "ffield")); elems = R.elements_in_order(sym); ok = True; det = []
    for mode in ("legacy", "smooth"):
        d = TMP / ("fd_" + mode); d.mkdir(); shutil.copy(case / "ffield", d / "ffield"); shutil.copy(P / f"inputs/control/{mode}.ctrl", d / "ctrl.txt")
        R.write_data(d / "data.lmp", sym, pos.tolist(), q, None, False, masses)
        L = [f"plugin load {PLUG}", "units real", "atom_style charge", "boundary f f f", "newton on", "atom_modify map array", "read_data data.lmp",
             "pair_style reaxff/metal ctrl.txt ereaxff yes electron_taper 5.7", f"pair_coeff * * ffield {' '.join(elems)}", "comm_modify cutoff 12.0", "thermo_style custom step pe", "thermo_modify format float %.14g", "thermo 1",
             "variable fx equal fx[1]", "variable fy equal fy[1]", "variable fz equal fz[1]"]
        sel = [(0, 0), (0, 2), (n - 1, 1), (3, 2)]   # (atom index, axis)
        L += ["run 0", "print \"FORCE ${fx} ${fy} ${fz}\""]
        for i, ax in sel:
            for s in (+1, -1):
                p = pos[i].copy(); p[ax] += s * h
                L += [f"set atom {i+1} x {p[0]:.15g} y {p[1]:.15g} z {p[2]:.15g}", "run 0"]
            L += [f"set atom {i+1} x {pos[i][0]:.15g} y {pos[i][1]:.15g} z {pos[i][2]:.15g}", f"variable f{i}_{ax} equal f{'xyz'[ax]}[{i+1}]", "run 0", f"print \"FREF {i} {ax} ${{f{i}_{ax}}}\""]
        (d / "in.lmp").write_text("\n".join(L) + "\n")
        rc, o = sh([LMP, "-in", "in.lmp", "-log", "none", "-nocite"], cwd=d)
        if rc: return False, o[-300:]
        isnum = lambda x: re.fullmatch(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", x) is not None
        pes = [float(l.split()[1]) for l in o.splitlines() if len(l.split()) == 2 and l.split()[0].isdigit() and isnum(l.split()[1])]
        refs = {(int(l.split()[1]), int(l.split()[2])): float(l.split()[3]) for l in o.splitlines() if l.startswith("FREF")}
        # pes order: base, then per selection: +h, -h, base again
        k = 1; worst = 0.0
        for j, (i, ax) in enumerate(sel):
            ep, em = pes[1 + 3 * j], pes[2 + 3 * j]
            fd = -(ep - em) / (2 * h); worst = max(worst, abs(fd - refs[(i, ax)]))
        det.append(f"{mode}: max|F_fd - F| = {worst:.1e} kcal/mol/A over {len(sel)} coords (CONTROL_USED {mode})"); ok &= worst < 2e-3
    return ok, "; ".join(det)

@check("smooth_vs_ereac")
def smooth_vs_ereac():
    fc = TMP / "fc"; rc, o = sh([sys.executable, P / "tools/discontinuity/frames_compare.py", fc], timeout=1800)
    if rc: return False, o[-300:]
    det = []; ok = True
    for mode in ("legacy", "smooth"):
        rc, o = sh([sys.executable, P / "tools/discontinuity/ereac_smooth_check.py", mode, fc, TMP / ("ereac_" + mode), "25"], timeout=3000)
        m = re.search(r"max \|E_ereac - E_adapter\| = ([\d.e+-]+) kcal/mol; max \|F diff\| = ([\d.e+-]+)", o)
        if not m: return False, o[-300:]
        ok &= float(m.group(1)) < 2e-3 and float(m.group(2)) < 5e-2; det.append(f"{mode}: dE {m.group(1)} dF {m.group(2)} (CONTROL_USED {mode})")
    return ok, "; ".join(det) + " (floor = 1e-5 A geo quantisation in ereac)"

@check("conservation_p1")
def conservation_p1():
    import numpy as np
    spans = {}
    for dt in (0.25, 0.125):
        env = dict(ENV, CTRL="ctrl.txt", CTRLTXT=(P / "inputs/control/smooth.ctrl").read_text())
        d = f"cons_{dt}"; rc, o = sh([sys.executable, P / "tools/discontinuity/run_nve.py", "Ag7_1_md", TMP / d, int(1000 / dt), dt, 300, 4928459, 0], env=env, cwd=TMP)
        if rc: return False, o[-300:]
        rc, o = sh([sys.executable, P / "tools/discontinuity/analyze.py", TMP / d]); spans[dt] = float(re.search(r"etot span ([\d.]+)", o).group(1))
    ratio = spans[0.25] / spans[0.125]
    return 3.0 < ratio < 5.0 and spans[0.25] < 0.05, f"CONTROL_USED smooth; span {spans[0.25]:.4f} (dt 0.25) -> {spans[0.125]:.4f} (dt 0.125), ratio {ratio:.2f} (dt^2 -> 4)"

@check("classical_atoms")
def classical_atoms():
    ipi = TREE / "forks/i-pi"; sel = [] if a.full else ["-k", "L_P1 or L_P2 or restart or exact_ring or L_P4_nve"]
    env = dict(os.environ, OMP_NUM_THREADS="1", IPI_UPSTREAM_ROOT=str(TREE / "scratch/ipi_upstream"), PYTHONPATH=str(ipi), PATH=f"{ipi/'bin'}:{os.environ['PATH']}")
    rc, o = sh([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "ipi_tests/classical_atoms/test_classical_atoms.py", *sel], cwd=ipi, env=env, timeout=3600)
    return rc == 0, ("full file" if a.full else "quick subset") + ": " + (o.strip().splitlines() or ["(no output)"])[-1]

@check("ipi_regression")
def ipi_regression():
    d = Path(a.ipi_reg_dir)
    if not (d / "ipi_tests/regression_tests/tests/CLASSICAL_ATOMS/nvt_pilel_obabo/coupled_harmonic_python/ref_simulation.out").exists():
        return False, f"{d}: ref_simulation.out missing (they are ignored by .gitignore '*.out' unless added with git add -f; branch regression-runner-diag has them)"
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONPATH=str(d), PATH=f"{d/'bin'}:{os.environ['PATH']}")
    rc, o = sh([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-o", "faulthandler_timeout=300", "ipi_tests/regression_tests/test_run.py", "-k", "CLASSICAL"], cwd=d, env=env, timeout=1800)
    return rc == 0, (o.strip().splitlines() or ["(no output)"])[-1]

@check("ipi_lammps_short")
def ipi_lammps_short():
    import numpy as np
    ref = PROJ / "reference_data/ipi_lammps_short/sim.out"; out = TMP / "ipilmp"
    rc, o = sh([sys.executable, P / "tools/run_pimd_lammps_ctrl.py", "--out", out, "--beads", "4", "--clients", "1", "--steps", "200", "--dt-fs", "0.25", "--mode", "nve", "--seed", "5", "--control", "legacy", "--latency", "1e-5", "--stride", "10", "--pos-stride", "100000"], timeout=900)
    if rc: return False, o[-300:]
    new = np.loadtxt(out / "sim.out")
    if a.update_ipi_reference:
        ref.parent.mkdir(parents=True, exist_ok=True); shutil.copy(out / "sim.out", ref); return True, f"reference stored at {ref} ({len(new)} rows)"
    if not ref.exists(): return None, "no stored output yet (run once with --update-ipi-reference)"
    old = np.loadtxt(ref); err = np.abs(new[:, 2:5] - old[:, 2:5]).max() / np.abs(old[:, 2:5]).max()
    return err < 1e-8, f"CONTROL_USED legacy; max rel diff of conserved/potential/kinetic vs stored = {err:.1e} (a regression baseline, not an independent reference)"

@check("schrodinger")
def schrodinger():
    import numpy as np
    ref = PROJ / "reference_data/tunneling_schrod"; det = []
    out = TMP / "s1d.json"; rc, o = sh([sys.executable, P / "inputs/tunneling_ag/schrod1d.py", ref / "E_F.npy", ref / "zf.npy", out])
    if rc: return False, o[-300:]
    new, old = json.loads(out.read_text()), json.loads((PROJ / "validation/tunneling_ag/data/schrod1d.json").read_text())
    worst = max(abs(new[k]["E0"] - old[k]["E0"]) + abs(new[k]["frac_mid"] - old[k]["frac_mid"]) for k in ("4.0", "1.0", "0.5", "0.1"))
    ok = worst < 1e-6; det.append(f"1-D solver re-run vs stored: max diff {worst:.1e}")
    for f in sorted((PROJ / "validation/tunneling_ag/data/stage2").glob("pimd_m*_P*_base*.json")):
        d = json.loads(f.read_text())
        if d.get("names") and d.get("pimd") and d.get("se") and d.get("solver"):
            z = [abs(p - s) / e for p, s, e in zip(d["pimd"], d["solver"], d["se"]) if e > 0]
            ok &= max(z) < 3.0; det.append(f"{f.stem}: max |PIMD - solver|/se = {max(z):.1f}")
    return ok, "; ".join(det)

@check("eparams")
def eparams():
    if not a.mr_build: return None, "pass --mr-build <build dir of the checkout with the eparams change>"
    b = Path(a.mr_build); rc, o = sh(["ctest", "-R", "compat_flags", "--output-on-failure"], cwd=b, timeout=900)
    if rc: return False, "ctest compat_flags: " + o[-200:]
    ag = P / "reference_data/phase4_frames_ag/Ag7_1_md"
    rc, o = sh([sys.executable, LINK / "scratch/mr_merged/tests/lammps/run_eparams.py", "--lmp", LMP, "--plugin", b / "plugin/reaxmetaladapterplugin.so", "--ag-ffield", ag / "ffield", "--ag-geo", ag / "geo", "--ag-frames", ag / "frames.json",
                "--cho-ffield", P / "reference_data/phase4_frames/ffield", "--table", LINK / "scratch/mr_merged/tests/data/legacy_table.eparams"], timeout=900)
    return rc == 0, "run_eparams.py: " + ("; ".join(re.findall(r"(?:PASS|FAIL) .*", o))[:300] or o[-200:])

@check("periodic_ag")
def periodic_ag():
    rc, o = sh([sys.executable, P / "tools/discontinuity/periodic_ag_check.py", "22.0", TMP / "per22", "5"], timeout=1800, env=dict(ENV, EB_LINK=str(LINK)))
    m = re.findall(r"WORST.*", o)
    return rc == 0 and "RESULT: PASS" in o, "cell 22 A (> 2 x 10 A cutoff): " + (m[0] if m else o[-200:])

t0 = time.time(); fails = 0
for f in CHECKS:
    if a.only and f.check_name not in a.only: continue
    t1 = time.time()
    try: ok, detail = f()
    except Exception as e: ok, detail = False, f"exception: {type(e).__name__}: {e}"
    status = "SKIP" if ok is None else ("PASS" if ok else "FAIL"); fails += status == "FAIL"
    print(f"{status:4s}  {f.check_name:16s} ({time.time()-t1:5.0f} s)  {detail}", flush=True)
print(f"\n{'ALL CHECKS PASSED' if not fails else str(fails) + ' CHECK(S) FAILED'} in {time.time()-t0:.0f} s; scratch {TMP}")
sys.exit(1 if fails else 0)
