#!/usr/bin/env python3
"""Donor-barrier-acceptor Ag/El system (gap 4.0 A) through i-PI (fork classical-atoms) + LAMMPS/Metal-Reax clients (fix ipi, unix sockets).
Ag fixed (<fixatoms>); the transferring electron (last atom, id 11) is a P-bead ring of mass --el-mass; the other four El are classical (1 amu)
via <classical_atoms> (default) or frozen (--freeze-others, benchmark for the 1-electron Schroedinger solve) or rings (--all-rings).
--start pit|middle  : initial z of the transferring electron (-1.975 or 0.0 A).   --pos-stride : output stride of the xyz trajectory (bead positions).
--v0-kcal/--fwhm/--barrier-z : optional Gaussian on |z| (not used in stage 1). --dump-forces: each client dumps x,f at every step (A1-style check).
"""
import argparse, os, re, shutil, subprocess, sys, time
from pathlib import Path
HERE = Path(__file__).resolve(); PROJ = HERE.parents[2]; ROOT = PROJ.parent
IPI = Path(os.environ.get("IPI_ROOT", ROOT / "forks/i-pi")); LMP = ROOT / "external/install_lammps/bin/lmp"
PLUGIN = Path(os.environ.get("REAX_PLUGIN", ROOT / "forks/metal-reax-lammps/build-ns/plugin/reaxmetaladapterplugin.so"))
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(ROOT / "forks/metal-reax-lammps/tests/lammps"))
from build_system import build
import run_ereaxff as R

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True)
    ap.add_argument("--beads", type=int, default=8); ap.add_argument("--clients", type=int, default=1)
    ap.add_argument("--steps", type=int, default=1000); ap.add_argument("--dt-fs", type=float, default=0.25); ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--mode", default="nvt"); ap.add_argument("--thermo", default="pile_l"); ap.add_argument("--tau-fs", type=float, default=100.0)
    ap.add_argument("--el-mass", type=float, default=1.0); ap.add_argument("--start", default="pit"); ap.add_argument("--gap", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=7); ap.add_argument("--stride", type=int, default=1); ap.add_argument("--pos-stride", type=int, default=10)
    ap.add_argument("--freeze-others", action="store_true"); ap.add_argument("--all-rings", action="store_true"); ap.add_argument("--no-fix-ag", action="store_true")
    ap.add_argument("--dump-forces", action="store_true"); ap.add_argument("--timeout", type=int, default=86400)
    ap.add_argument("--from-restart", default="", help="continue from a checkpoint dir (RESTART): sets mode/steps/dt")
    ap.add_argument("--checkpoint-stride", type=int, default=0)
    ap.add_argument("--v0-kcal", type=float, default=0.0, help="Stage 2: z-only Gaussian V0*exp(-z^2/(2 sigma^2)) (kcal/mol) on the quantum electron, centred at z = 0 (python gaussian_barrier driver)")
    ap.add_argument("--sigma-z", type=float, default=0.10, help="Gaussian sigma (A)"); ap.add_argument("--keep-xyz", action="store_true", help="keep the bead xyz files (default: reduced to ring.npy)"); ap.add_argument("--latency", type=float, default=None, help="i-PI ffsocket polling latency (s); default: i-PI default 1e-4");  ap.add_argument("--start-z", type=float, default=None, help="initial z of the quantum electron (A); overrides --start")
    a = ap.parse_args()
    out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    sym, pos, q, ia, zl, zr = build("pair", a.gap)
    pos = [list(map(float, p)) for p in pos]
    pos[ia][2] = a.start_z if a.start_z is not None else ({"pit": -1.975, "middle": 0.0}[a.start] if a.gap == 4.0 else zl)
    n = len(sym); ffp = PROJ / "inputs/ag_el/ffield_Ag_El"; masses = R.ffield_masses(str(ffp))
    ag = [i for i, s in enumerate(sym) if s == "Av"]; els = [i for i, s in enumerate(sym) if s == "El"]; others = [i for i in els if i != ia]
    ipi_mass = [masses["AV"] if s == "Av" else (1.0 if i != ia and not a.all_rings else a.el_mass) for i, s in enumerate(sym)]
    (out / "init.xyz").write_text("\n".join([str(n), "# CELL(abcABC): 200 200 200 90 90 90 positions{angstrom} cell{angstrom}"] + ["%s %.12f %.12f %.12f" % ("Ag" if s == "Av" else "El", *p) for s, p in zip(sym, pos)]) + "\n")
    elems = R.elements_in_order(sym); sock = "bl%d" % (os.getpid() % 100000)
    for k in range(a.clients):
        d = out / f"client{k}"; d.mkdir(exist_ok=True); shutil.copy(ffp, d / "ffield")
        R.write_data(d / "data.lmp", sym, pos, q, None, False, masses)
        extra = ""
        if a.dump_forces: extra = "dump fd all custom 1 fd.dump id x y z fx fy fz\ndump_modify fd format float %.16g sort id\nthermo_style custom step pe\nthermo_modify format float %.16g\nthermo 1\n"
        plug = os.path.relpath(PLUGIN, d)
        (d / "in.ipi").write_text(f"plugin load {plug}\nunits real\natom_style charge\nboundary f f f\nnewton on\nread_data data.lmp\n"
            f"pair_style reaxff/metal NULL ereaxff yes electron_taper 5.7\npair_coeff * * ffield {' '.join(elems)}\ncomm_modify cutoff 12.0\n{extra}fix ip all ipi {sock} 32123 unix\nrun 100000000\n")
    fixed = [] if a.no_fix_ag else list(ag)
    classical = [] if a.all_rings else others
    if a.freeze_others: fixed += others; classical = []
    fx = "<fixatoms> [%s] </fixatoms>" % ", ".join(map(str, fixed)) if fixed else ""
    cl = "<classical_atoms> [%s] </classical_atoms>" % ", ".join(map(str, classical)) if classical else ""
    lat = "" if a.latency is None else f"<latency> {a.latency} </latency>"
    bsock = sock + "b"; fb = ffb = ""
    if a.v0_kcal:
        ffb = f"<ffsocket name='barrier' mode='unix' pbc='false'>{lat}<address>{bsock}</address><consolidate_messages> True </consolidate_messages><batch_size> {a.beads} </batch_size></ffsocket>"
        fb = "<force forcefield='barrier'/>"
    props = f"[ step, time{{femtosecond}}, conserved{{electronvolt}}, potential{{electronvolt}}, kinetic_md{{electronvolt}}, kinetic_cv{{electronvolt}}, temperature{{kelvin}}, kinetic_cv{{electronvolt}}({ia}) ]"
    ck = a.checkpoint_stride or a.steps
    xml = f"""<simulation verbosity='low'>
  <output prefix='sim'>
    <properties stride='{a.stride}' filename='out'> {props} </properties>
    <trajectory filename='pos' stride='{a.pos_stride}' format='xyz'> positions{{angstrom}} </trajectory>
    <checkpoint stride='{ck}' overwrite='True'/>
  </output>
  <total_steps>{a.steps}</total_steps>
  <prng><seed>{a.seed}</seed></prng>
  <ffsocket name='lmp' mode='unix' pbc='false'>{lat}<address>{sock}</address></ffsocket>
  {ffb}
  <system>
    <initialize nbeads='{a.beads}'>
      <cell mode='abc' units='angstrom'> [200, 200, 200] </cell>
      <file mode='xyz' units='angstrom'> init.xyz </file>
      <masses mode='manual' units='dalton'> [{', '.join(repr(x) for x in ipi_mass)}] </masses>
      <velocities mode='thermal' units='kelvin'> {a.temp} </velocities>
    </initialize>
    <forces><force forcefield='lmp'/>{fb}</forces>
    <motion mode='dynamics'><fixcom> False </fixcom>{fx}<dynamics mode='{a.mode}'><timestep units='femtosecond'> {a.dt_fs} </timestep>
      <thermostat mode='{a.thermo}'><tau units='femtosecond'> {a.tau_fs} </tau></thermostat>
      {cl}</dynamics></motion>
    <ensemble><temperature units='kelvin'> {a.temp} </temperature></ensemble>
  </system>
</simulation>
"""
    if a.from_restart:
        rs = (Path(a.from_restart).resolve() / "RESTART").read_text()
        step0 = int(re.search(r"<step>\s*(\d+)\s*</step>", rs).group(1))
        rs = re.sub(r"<total_steps>\s*\d+\s*</total_steps>", f"<total_steps>{step0 + a.steps}</total_steps>", rs)
        rs = re.sub(r"(<dynamics mode=')\w+(')", r"\g<1>%s\2" % a.mode, rs, count=1)
        rs = re.sub(r"(<address>)[^<]*(</address>)", r"\g<1>%s\2" % sock, rs)
        rs = re.sub(r"<output prefix='[^']*'>", "<output prefix='sim'>", rs, count=1)
        rs = re.sub(r"<timestep>[^<]*</timestep>", "<timestep> %.12e </timestep>" % (a.dt_fs * 41.341373575751), rs, count=1)
        (out / "input.xml").write_text(rs)
    else:
        (out / "input.xml").write_text(xml)
    env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(LMP.parents[1] / "lib"), PYTHONPATH=str(IPI), PYTHONUNBUFFERED="1")
    srv = subprocess.Popen([sys.executable, str(IPI / "bin" / "i-pi"), "input.xml"], cwd=out, env=env, stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
    time.sleep(4)
    clis = [subprocess.Popen([str(LMP), "-in", "in.ipi", "-log", "none", "-nocite"], cwd=out / f"client{k}", env=env, stdout=open(out / f"client{k}/lmp.log", "w"), stderr=subprocess.STDOUT) for k in range(a.clients)]
    def _git(r):
        g = lambda *x: subprocess.run(["git", "-C", str(r)] + list(x), capture_output=True, text=True).stdout.strip()
        return "%s dirty=%d" % (g("rev-parse", "--short", "HEAD"), len([l for l in g("status", "--short").splitlines() if not l.startswith("??")]))
    import platform
    (out / "run_info.txt").write_text("\n".join(["command: " + " ".join(sys.argv), "host: " + platform.node() + " " + platform.machine(), "date: " + time.strftime("%Y-%m-%d %H:%M:%S"),
        "CONTROL_USED NULL(legacy)   pair_style reaxff/metal NULL ereaxff yes electron_taper 5.7   (Ag fixed: %s)" % (not a.no_fix_ag),
        "seed %d  beads %d  dt_fs %g  el_mass %g  temp %g  thermo %s tau_fs %g  latency %s" % (a.seed, a.beads, a.dt_fs, a.el_mass, a.temp, a.thermo, a.tau_fs, a.latency),
        "design V0_kcal %g sigma_A %g  start_z %s  freeze_others %s  mode %s  from_restart %s" % (a.v0_kcal, a.sigma_z, a.start_z, a.freeze_others, a.mode, a.from_restart or "-"),
        "commits: project " + _git(PROJ) + " | i-pi " + _git(IPI) + " | metal-reax-lammps " + _git(ROOT / "forks/metal-reax-lammps") + " | plugin " + str(PLUGIN)]) + "\n")
    print("CONTROL_USED NULL(legacy); seed %d" % a.seed, flush=True)
    if a.v0_kcal:
        sys.path.insert(0, str(IPI)); from ipi.pes.gaussian_barrier import barrier_params
        prm = barrier_params(a.v0_kcal, a.sigma_z * 2 * (2 * __import__('math').log(2)) ** 0.5, [0, 0, 0], [0, 0, 1], ia, ia)
        clis.append(subprocess.Popen([sys.executable, str(IPI / "bin/i-pi-py_driver"), "-m", "gaussian_barrier", "-u", "-a", bsock, "-o", prm], cwd=out, env=env, stdout=open(out / "driver_barrier.log", "w"), stderr=subprocess.STDOUT))
    t0 = time.time()
    try: srv.wait(timeout=a.timeout)
    finally:
        for c in clis: c.terminate()
        if srv.poll() is None: srv.kill()
    if not a.keep_xyz and srv.returncode == 0:                      # reduce the bead trajectories of the quantum electron to one small array (P, frames, 3), then drop the xyz files
        import numpy as np
        fs = sorted(out.glob("sim.pos_*.xyz"), key=lambda f: int(f.stem.split("_")[1])); arr = []
        for f in fs:
            L = f.read_text().splitlines(); nn = len(sym) + 2
            arr.append([[float(v) for v in L[i + 2 + ia].split()[1:4]] for i in range(0, len(L) - nn + 1, nn)])
        n = min(len(x) for x in arr); np.save(out / "ring.npy", np.array([x[:n] for x in arr], dtype=np.float32))
        for f in fs: f.unlink()
    print("i-PI exit %s after %.1f s (P=%d, %d clients, %d steps)" % (srv.returncode, time.time() - t0, a.beads, a.clients, a.steps)); sys.exit(srv.returncode)
main()
