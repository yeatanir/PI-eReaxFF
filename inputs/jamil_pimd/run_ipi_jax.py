#!/usr/bin/env python3
"""Build an i-PI input for the Jamil graphene/electron system and run it with the JAX-kernel driver (plan phase 2).

Carbons are classical (<classical_atoms>) unless --all-ring; every El is a P-bead ring polymer (mass 1 amu).
All P beads go to the driver as one batched request (ffsocket batch_size = P).

Example:  python3 run_ipi_jax.py --geo <geo> --ffield <ffield> --out <dir> --beads 4 --steps 300 --dt-fs 0.1 --mode nve
"""
import argparse
import os
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve()
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "tools"))
sys.path.insert(0, str(PROJ.parent / "forks" / "i-pi"))
IPI = PROJ.parent / "forks" / "i-pi"
VENV_PY = "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_pimd_modular/.venv-fast/bin/python"
MASS = {"C": 12.011, "El": 1.0}

TEMPLATE = """<simulation verbosity='low'>
  <output prefix='sim'>
    <properties stride='{stride}' filename='out'> [ step, time{{femtosecond}}, conserved{{electronvolt}}, potential{{electronvolt}}, kinetic_md{{electronvolt}}, kinetic_cv{{electronvolt}}, temperature{{kelvin}} ] </properties>
    <trajectory filename='cen' stride='{tstride}' format='xyz'> x_centroid{{angstrom}} </trajectory>
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
        <thermostat mode='{thermo}'>
          <tau units='femtosecond'> {tau} </tau>
        </thermostat>
        {classical}
      </dynamics>
    </motion>
    <ensemble><temperature units='kelvin'> {temp} </temperature></ensemble>
  </system>
</simulation>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geo", required=True)
    ap.add_argument("--ffield", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--beads", type=int, default=4)
    ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--dt-fs", type=float, default=0.1)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--mode", default="nve", choices=["nve", "nvt"])
    ap.add_argument("--thermo", default="pile_l")
    ap.add_argument("--tau-fs", type=float, default=100.0)
    ap.add_argument("--all-ring", action="store_true", help="make carbons ring polymers too (no classical_atoms)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=86400)
    a = ap.parse_args()

    import ipi_ereaxff_jax_driver as drv
    symbols, pos, q, _ = drv.read_geo(a.geo)
    out = pathlib.Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    lines = [str(len(symbols)), "# CELL(abcABC): 200 200 200 90 90 90 positions{angstrom} cell{angstrom}"]
    lines += ["%s %.8f %.8f %.8f" % (s, *p) for s, p in zip(symbols, pos)]
    (out / "init.xyz").write_text("\n".join(lines) + "\n")
    carbons = [i for i, s in enumerate(symbols) if s == "C"]
    classical = "" if a.all_ring else "<classical_atoms> [%s] </classical_atoms>" % ", ".join(map(str, carbons))
    address = "erff%d" % (os.getpid() % 100000)
    (out / "input.xml").write_text(TEMPLATE.format(
        stride=a.stride, tstride=max(a.stride, 10), steps=a.steps, seed=a.seed, beads=a.beads, address=address,
        masses=", ".join(str(MASS[s]) for s in symbols), temp=a.temp, mode=a.mode, dt=a.dt_fs,
        thermo=a.thermo, tau=a.tau_fs, classical=classical))
    env = dict(os.environ, PYTHONPATH=str(IPI), PYTHONUNBUFFERED="1", OMP_NUM_THREADS="4", XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=4", EREAXFF_FFIELD=str(a.ffield), EREAXFF_GEO=str(a.geo))
    srv = subprocess.Popen([sys.executable, str(IPI / "bin" / "i-pi"), "input.xml"], cwd=out, env=env,
                           stdout=open(out / "ipi.log", "w"), stderr=subprocess.STDOUT)
    time.sleep(3)
    cl = subprocess.Popen([VENV_PY, str(IPI / "bin" / "i-pi-py_driver"), "-m", "custom", "-P",
                           str(PROJ / "tools" / "ipi_ereaxff_jax_driver.py"), "-u", "-a", address],
                          cwd=out, env=env, stdout=open(out / "driver.log", "w"), stderr=subprocess.STDOUT)
    t0 = time.time()
    try:
        srv.wait(timeout=a.timeout)
    finally:
        cl.terminate()
        if srv.poll() is None:
            srv.kill()
    print("i-PI exit %s after %.0f s" % (srv.returncode, time.time() - t0))
    sys.exit(srv.returncode)


if __name__ == "__main__":
    main()
