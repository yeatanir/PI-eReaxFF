#!/usr/bin/env python3
"""Run i-PI PIMD on an Ag/El cluster with the exact legacy ereac as force engine (tools/ipi_ereac_driver.py).
Ag (type Av) cations are classical via <classical_atoms> unless --all-ring; El are P-bead ring polymers.
Usage: run_ipi_ereac.py <case_dir (control, ffield, geo, fort.20, fort.35)> <out_dir> --beads P --steps N [--mode nve|nvt] ...
"""
import argparse, os, pathlib, subprocess, sys, time
HERE = pathlib.Path(__file__).resolve(); PROJ = HERE.parents[2]
IPI = PROJ.parent / "forks" / "i-pi"
EREAC = PROJ.parent / "scratch/phase0/ereac_export_build/ereac"
MASS = {"Av": 107.8682, "El": 1.0}
T = """<simulation verbosity='low'>
  <output prefix='sim'>
    <properties stride='{stride}' filename='out'> [ step, time{{femtosecond}}, conserved{{electronvolt}}, potential{{electronvolt}}, kinetic_md{{electronvolt}}, kinetic_cv{{electronvolt}}, temperature{{kelvin}}, temperature(El){{kelvin}} ] </properties>
    <trajectory filename='pos' stride='{tstride}' format='xyz'> positions{{angstrom}} </trajectory>
  </output>
  <total_steps>{steps}</total_steps>
  <prng><seed>{seed}</seed></prng>
  <ffsocket name='eRFF' mode='unix' pbc='false'>
    <address>{address}</address>
    <consolidate_messages> True </consolidate_messages>
    <batch_size> {beads} </batch_size>
  </ffsocket>
  <system>
    <initialize nbeads='{beads}'>
      <cell mode='abc' units='angstrom'> [200, 200, 200] </cell>
      <file mode='xyz' units='angstrom'> init.xyz </file>
      <masses mode='manual' units='dalton'> [{masses}] </masses>
      <velocities mode='thermal' units='kelvin'> {temp} </velocities>
    </initialize>
    <forces><force forcefield='eRFF'/></forces>
    <motion mode='dynamics'>
      <fixcom> False </fixcom>
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
    ap = argparse.ArgumentParser(); ap.add_argument("case"); ap.add_argument("out")
    ap.add_argument("--beads", type=int, default=8); ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--dt-fs", type=float, default=0.25); ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--mode", default="nve"); ap.add_argument("--thermo", default="pile_l"); ap.add_argument("--tau-fs", type=float, default=100.0)
    ap.add_argument("--all-ring", action="store_true"); ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--stride", type=int, default=1); ap.add_argument("--tstride", type=int, default=10)
    ap.add_argument("--workers", type=int, default=2); a = ap.parse_args()
    case, out = pathlib.Path(a.case).resolve(), pathlib.Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    sym, xyz = [], []
    for l in (case / "geo").read_text().splitlines():
        if l.startswith("HETATM"):
            sym.append(l[16:21].strip() or l.split()[2]); xyz.append((float(l[30:40]), float(l[40:50]), float(l[50:60])))
    sym = [l.split()[2] for l in (case / "geo").read_text().splitlines() if l.startswith("HETATM")]
    (out / "init.xyz").write_text("\n".join([str(len(sym)), "# CELL(abcABC): 200 200 200 90 90 90 positions{angstrom} cell{angstrom}"] + ["%s %.8f %.8f %.8f" % ("Ag" if s == "Av" else s, *p) for s, p in zip(sym, xyz)]) + "\n")
    ag = [i for i, s in enumerate(sym) if s == "Av"]
    cl = "" if a.all_ring else "<classical_atoms> [%s] </classical_atoms>" % ", ".join(map(str, ag))
    addr = "ag%d" % (os.getpid() % 100000)
    (out / "input.xml").write_text(T.format(stride=a.stride, tstride=a.tstride, steps=a.steps, seed=a.seed, beads=a.beads, address=addr,
        masses=", ".join(str(MASS[s]) for s in sym), temp=a.temp, mode=a.mode, dt=a.dt_fs, thermo=a.thermo, tau=a.tau_fs, classical=cl))
    env = dict(os.environ, PYTHONPATH=str(IPI), PYTHONUNBUFFERED="1", OMP_NUM_THREADS="1", EREAC_BIN=str(EREAC), EREAC_CASE=str(case),
               EREAC_WORKERS=str(a.workers), EREAC_TMP=str(out / "drv_tmp"))
    srv = subprocess.Popen([sys.executable, str(IPI / "bin/i-pi"), "input.xml"], cwd=out, env=env, stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
    time.sleep(3)
    cli = subprocess.Popen([sys.executable, str(IPI / "bin/i-pi-py_driver"), "-m", "custom", "-P", str(PROJ / "tools/ipi_ereac_driver.py"), "-u", "-a", addr],
                           cwd=out, env=env, stdout=open(out / "driver.log", "w"), stderr=subprocess.STDOUT)
    t0 = time.time()
    try: srv.wait()
    finally:
        cli.terminate()
        if srv.poll() is None: srv.kill()
    print("i-PI exit", srv.returncode, "%.0f s" % (time.time() - t0)); sys.exit(srv.returncode)
main()
