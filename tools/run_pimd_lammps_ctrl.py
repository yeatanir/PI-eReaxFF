#!/usr/bin/env python3
"""Ag/El path-integral MD: i-PI server (fork classical-atoms) + LAMMPS clients (pair reaxff/metal ereaxff, fix ipi, unix sockets).
Ag (type Av) cations are classical via <classical_atoms>; every El is a P-bead ring polymer. Every client runs in its OWN directory.

Stage 1 (equilibration):  --mode nvt --thermo pile_l  [--init-frame k]
Stage 2 (NVE from the equilibrated state): --from-restart <stage1 dir> --mode nve  (edits the i-PI checkpoint: dynamics mode, step counter)
--clients N : number of LAMMPS clients (default min(P, 4)); with N < P each client serves several beads per step (CPU cap).
--mass-el   : dynamical El mass in i-PI (amu). The LAMMPS data-file El mass is set by --data-mass-el (default: ffield mass) and must not matter.
--dump-forces : every client dumps x,f (%.16g) and thermo pe at every step, for the A1 force check.
"""
import argparse, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).absolute(); PROJ = HERE.parents[1]
IPI = PROJ.parent / "forks" / "i-pi"
LMPT = PROJ.parent / "forks" / "metal-reax-lammps" / "tests" / "lammps"
LMP = PROJ.parent / "external" / "install_lammps" / "bin" / "lmp"
PLUGIN = PROJ.parent / "forks" / "metal-reax-lammps" / "build-ns" / "plugin" / "reaxmetaladapterplugin.so"
sys.path.insert(0, str(LMPT))
import run_ereaxff as R  # noqa: E402

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=str(PROJ / "reference_data/phase4_frames_ag/Ag7_1_md")); ap.add_argument("--out", required=True)
    ap.add_argument("--beads", type=int, default=8); ap.add_argument("--clients", type=int, default=0)
    ap.add_argument("--steps", type=int, default=100); ap.add_argument("--dt-fs", type=float, default=0.25); ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--mode", default="nvt"); ap.add_argument("--thermo", default="pile_l"); ap.add_argument("--tau-fs", type=float, default=100.0)
    ap.add_argument("--mass-el", type=float, default=1.0); ap.add_argument("--data-mass-el", type=float, default=0.0)
    ap.add_argument("--taper", type=float, default=5.7); ap.add_argument("--seed", type=int, default=7); ap.add_argument("--init-frame", type=int, default=0)
    ap.add_argument("--from-restart", default=""); ap.add_argument("--stride", type=int, default=1); ap.add_argument("--pos-stride", type=int, default=10)
    ap.add_argument("--control", default="none", choices=["none", "legacy", "smooth"], help="LAMMPS control file (inputs/control/*.ctrl); none = built-in defaults (= legacy). smooth for mobile atoms.")
    ap.add_argument("--latency", type=float, default=None, help="i-PI ffsocket polling latency (s)")
    ap.add_argument("--dump-forces", action="store_true"); ap.add_argument("--timeout", type=int, default=86400)
    a = ap.parse_args()
    case = Path(a.case).resolve(); out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    nclient = a.clients or min(a.beads, 4)
    sym, _, q, _ = R.read_geo(case / "geo")
    pos = json.loads((case / "frames.json").read_text())["frames"][a.init_frame]["pos"]
    n = len(sym); masses = R.ffield_masses(str(case / "ffield"))
    dmass = dict(masses)
    if a.data_mass_el > 0: dmass["EL"] = a.data_mass_el
    ipi_mass = [a.mass_el if s.upper() == "EL" else masses[s.upper()] for s in sym]
    shutil.copy(case / "ffield", out / "ffield")
    (out / "init.xyz").write_text("\n".join([str(n), "# CELL(abcABC): 200 200 200 90 90 90 positions{angstrom} cell{angstrom}"] + ["%s %.12f %.12f %.12f" % (s, *p) for s, p in zip(sym, pos)]) + "\n")
    elems = R.elements_in_order(sym)
    ctrl_name = "NULL" if a.control == "none" else "ctrl.txt"
    (out / "control_used.txt").write_text(f"control={a.control}\n")
    print(f"CONTROL_USED {a.control}", flush=True)
    sock = "ag%d" % (os.getpid() % 100000)
    for k in range(nclient):
        d = out / f"client{k}"; d.mkdir(exist_ok=True); shutil.copy(case / "ffield", d / "ffield")
        R.write_data(d / "data.lmp", sym, pos, q, None, False, dmass)
        if a.control != "none": shutil.copy(PROJ / f"inputs/control/{a.control}.ctrl", d / "ctrl.txt")
        extra = ""
        if a.dump_forces:
            extra = ("dump fd all custom 1 fd.dump id x y z fx fy fz\ndump_modify fd format float %.16g sort id\n"
                     "thermo_style custom step pe\nthermo_modify format float %.16g\nthermo 1\n")
        plug = os.path.relpath(PLUGIN, d) if ' ' in str(PLUGIN) else str(PLUGIN)
        (d / "in.ipi").write_text(f"""plugin load {plug}
units real
atom_style charge
boundary f f f
newton on
read_data data.lmp
pair_style reaxff/metal {ctrl_name} ereaxff yes electron_taper {a.taper}
pair_coeff * * ffield {' '.join(elems)}
comm_modify cutoff 12.0
{extra}fix ip all ipi {sock} 32123 unix
run 100000000
""")
    ag = [i for i, s in enumerate(sym) if s.upper() == "AV"]
    cl = "<classical_atoms> [%s] </classical_atoms>" % ", ".join(map(str, ag))
    props = "[ step, time{femtosecond}, conserved{electronvolt}, potential{electronvolt}, kinetic_md{electronvolt}, kinetic_cv{electronvolt}, temperature{kelvin} ]"
    xml = f"""<simulation verbosity='low'>
  <output prefix='sim'>
    <properties stride='{a.stride}' filename='out'> {props} </properties>
    <trajectory filename='pos' stride='{a.pos_stride}' format='xyz'> positions{{angstrom}} </trajectory>
    <checkpoint stride='{a.steps}' overwrite='True'/>
  </output>
  <total_steps>{a.steps}</total_steps>
  <prng><seed>{a.seed}</seed></prng>
  <ffsocket name='lmp' mode='unix' pbc='false'>{'' if a.latency is None else f'<latency> {a.latency} </latency>'}<address>{sock}</address></ffsocket>
  <system>
    <initialize nbeads='{a.beads}'>
      <cell mode='abc' units='angstrom'> [200, 200, 200] </cell>
      <file mode='xyz' units='angstrom'> init.xyz </file>
      <masses mode='manual' units='dalton'> [{', '.join(repr(x) for x in ipi_mass)}] </masses>
      <velocities mode='thermal' units='kelvin'> {a.temp} </velocities>
    </initialize>
    <forces><force forcefield='lmp'/></forces>
    <motion mode='dynamics'><fixcom> False </fixcom><dynamics mode='{a.mode}'><timestep units='femtosecond'> {a.dt_fs} </timestep>
      <thermostat mode='{a.thermo}'><tau units='femtosecond'> {a.tau_fs} </tau></thermostat>
      {cl}</dynamics></motion>
    <ensemble><temperature units='kelvin'> {a.temp} </temperature></ensemble>
  </system>
</simulation>
"""
    if a.from_restart:   # NVE (or other) continuation: edit the checkpoint of the previous stage
        rs = (Path(a.from_restart).resolve() / "RESTART").read_text()
        step0 = int(re.search(r"<step>\s*(\d+)\s*</step>", rs).group(1))
        rs = re.sub(r"<total_steps>\s*\d+\s*</total_steps>", f"<total_steps>{step0 + a.steps}</total_steps>", rs)
        rs = re.sub(r"(<dynamics mode=')\w+(')", r"\g<1>%s\2" % a.mode, rs, count=1)
        rs = re.sub(r"(<address>)[^<]*(</address>)", r"\g<1>%s\2" % sock, rs)
        rs = re.sub(r"<output prefix='[^']*'>", "<output prefix='sim'>", rs, count=1)
        rs = re.sub(r"<timestep>[^<]*</timestep>", "<timestep> %.12e </timestep>" % (a.dt_fs * 41.341373575751), rs, count=1)   # checkpoint is in atomic units
        (out / "input.xml").write_text(rs)
    else:
        (out / "input.xml").write_text(xml)
    env = dict(os.environ, OMP_NUM_THREADS="1", DYLD_LIBRARY_PATH=str(LMP.parents[1] / "lib"), PYTHONPATH=str(IPI), PYTHONUNBUFFERED="1")
    srv = subprocess.Popen([sys.executable, str(IPI / "bin" / "i-pi"), "input.xml"], cwd=out, env=env, stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
    time.sleep(4)
    clis = [subprocess.Popen([str(LMP), "-in", "in.ipi", "-log", "none", "-nocite"], cwd=out / f"client{k}", env=env,
                             stdout=open(out / f"client{k}" / "lmp.log", "w"), stderr=subprocess.STDOUT) for k in range(nclient)]
    t0 = time.time()
    try:
        srv.wait(timeout=a.timeout)
    finally:
        for c in clis: c.terminate()
        if srv.poll() is None: srv.kill()
    print("i-PI exit %s after %.0f s (P=%d, %d clients, %d steps)" % (srv.returncode, time.time() - t0, a.beads, nclient, a.steps))
    sys.exit(srv.returncode)

main()
