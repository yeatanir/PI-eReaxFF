#!/usr/bin/env python3
"""i-PI run of the donor-barrier-acceptor Ag/El system: eReaxFF force (exact ereac driver) + optional Gaussian barrier as a
second <force> component (own ffsocket), Ag atoms fixed (<fixatoms>), El atoms P-bead rings.
Usage: run_ipi_barrier.py <out_dir> --kind pair --gap 4.5 --beads 4 --steps 100 [--v0-kcal 3 --fwhm 0.4 | --no-barrier] [--mode nve|nvt]
"""
import argparse, os, pathlib, subprocess, sys, time
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; PROJ = HERE.parents[1]; ROOT = PROJ.parent; IPI = ROOT / "forks/i-pi"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(IPI))
from build_system import write_case
from ipi.pes.gaussian_barrier import barrier_params
EREAC = ROOT / "scratch/phase0/ereac_export_build/ereac"
MASS = {"Av": 107.8682, "El": 1.0}
SHIFT = 40.0
T = """<simulation verbosity='low'>
  <output prefix='sim'>
    <properties stride='{stride}' filename='out'> [ step, time{{femtosecond}}, conserved{{electronvolt}}, potential{{electronvolt}}, kinetic_md{{electronvolt}}, temperature{{kelvin}} ] </properties>
    <trajectory filename='pos' stride='{tstride}' format='xyz'> positions{{angstrom}} </trajectory>
  </output>
  <total_steps>{steps}</total_steps>
  <prng><seed>{seed}</seed></prng>
  <ffsocket name='eRFF' mode='unix' pbc='false'>
    <address>{addr1}</address><consolidate_messages> True </consolidate_messages><batch_size> {beads} </batch_size>
  </ffsocket>
  {ffbarrier}
  <system>
    <initialize nbeads='{beads}'>
      <cell mode='abc' units='angstrom'> [200, 200, 200] </cell>
      <file mode='xyz' units='angstrom'> init.xyz </file>
      <masses mode='manual' units='dalton'> [{masses}] </masses>
      <velocities mode='thermal' units='kelvin'> {temp} </velocities>
    </initialize>
    <forces><force forcefield='eRFF'/>{forcebarrier}</forces>
    <motion mode='dynamics'>
      <fixcom> False </fixcom>
      <fixatoms> [{fixed}] </fixatoms>
      <dynamics mode='{mode}'>
        <timestep units='femtosecond'> {dt} </timestep>
        <thermostat mode='{thermo}'><tau units='femtosecond'> {tau} </tau></thermostat>
        {classical}
      </dynamics>
    </motion>
    <ensemble><temperature units='kelvin'> {temp} </temperature></ensemble>
  </system>
</simulation>
"""
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out")
    ap.add_argument("--kind", default="pair"); ap.add_argument("--gap", type=float, default=4.5)
    ap.add_argument("--beads", type=int, default=4); ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--dt-fs", type=float, default=0.25); ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--mode", default="nve"); ap.add_argument("--thermo", default="pile_l"); ap.add_argument("--tau-fs", type=float, default=100.0)
    ap.add_argument("--v0-kcal", type=float, default=3.0); ap.add_argument("--fwhm", type=float, default=0.4); ap.add_argument("--no-barrier", action="store_true")
    ap.add_argument("--barrier-z", type=float, default=0.0, help="barrier centre along z relative to the donor-acceptor midpoint (A)")
    ap.add_argument("--one-quantum-electron", action="store_true", help="only the transferring electron (last atom) is a ring polymer of mass --el-mass; the other El are classical (1 amu) via <classical_atoms>")
    ap.add_argument("--el-mass", type=float, default=1.0); ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--stride", type=int, default=1); ap.add_argument("--tstride", type=int, default=1); ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    out = pathlib.Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    case = out / "case"; sym, pos, q, ia, zl, zr = write_case(case, a.kind, a.gap, SHIFT)
    nag = sum(1 for s in sym if s == "Av")
    (out / "init.xyz").write_text("\n".join([str(len(sym)), "# CELL(abcABC): 200 200 200 90 90 90 positions{angstrom} cell{angstrom}"] +
        ["%s %.8f %.8f %.8f" % ("Ag" if s == "Av" else s, *(p + SHIFT)) for s, p in zip(sym, pos)]) + "\n")
    addr1, addr2 = "tb%d" % (os.getpid() % 100000), "tc%d" % (os.getpid() % 100000)
    ffb = fb = ""
    if not a.no_barrier:
        ffb = "<ffsocket name='barrier' mode='unix' pbc='false'><address>%s</address><consolidate_messages> True </consolidate_messages><batch_size> %d </batch_size></ffsocket>" % (addr2, a.beads)
        fb = "<force forcefield='barrier'/>"
    masses = [MASS["Av"] if s == "Av" else a.el_mass for s in sym]
    classical = ""
    if a.one_quantum_electron:
        others = [k for k, x in enumerate(sym) if x == "El" and k != len(sym) - 1]
        for k in others: masses[k] = MASS["El"]                      # classical El keep the force-field mass (1 amu)
        classical = "<classical_atoms> [%s] </classical_atoms>" % ", ".join(map(str, others))
    (out / "input.xml").write_text(T.format(stride=a.stride, tstride=a.tstride, steps=a.steps, seed=a.seed, beads=a.beads, addr1=addr1, ffbarrier=ffb,
        forcebarrier=fb, masses=", ".join(map(str, masses)), temp=a.temp, mode=a.mode, dt=a.dt_fs, thermo=a.thermo, tau=a.tau_fs, fixed=", ".join(map(str, range(nag))), classical=classical))
    env = dict(os.environ, PYTHONPATH=str(IPI), PYTHONUNBUFFERED="1", OMP_NUM_THREADS="1", EREAC_BIN=str(EREAC), EREAC_CASE=str(case),
               EREAC_WORKERS=str(a.workers), EREAC_TMP=str(out / "drv_tmp"))
    srv = subprocess.Popen([sys.executable, str(IPI / "bin/i-pi"), "input.xml"], cwd=out, env=env, stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
    time.sleep(3)
    clients = [subprocess.Popen([sys.executable, str(IPI / "bin/i-pi-py_driver"), "-m", "custom", "-P", str(PROJ / "tools/ipi_ereac_driver.py"), "-u", "-a", addr1],
                                cwd=out, env=env, stdout=open(out / "driver_ereac.log", "w"), stderr=subprocess.STDOUT)]
    if not a.no_barrier:
        prm = barrier_params(a.v0_kcal, a.fwhm, [SHIFT, SHIFT, SHIFT + a.barrier_z], [0, 0, 1], nag, len(sym) - 1)
        clients.append(subprocess.Popen([sys.executable, str(IPI / "bin/i-pi-py_driver"), "-m", "gaussian_barrier", "-u", "-a", addr2, "-o", prm],
                                        cwd=out, env=env, stdout=open(out / "driver_barrier.log", "w"), stderr=subprocess.STDOUT))
    t0 = time.time()
    try: srv.wait()
    finally:
        for c in clients: c.terminate()
        if srv.poll() is None: srv.kill()
    print("i-PI exit", srv.returncode, "%.0f s" % (time.time() - t0)); sys.exit(srv.returncode)
main()
