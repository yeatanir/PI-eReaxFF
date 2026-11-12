import importlib.util, pathlib, sys
IPI = pathlib.Path("/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_build/forks/i-pi")
spec = importlib.util.spec_from_file_location("tca", IPI / "ipi_tests/classical_atoms/test_classical_atoms.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
# physical parameters of the in-house test, converted to a.u.: 1 eV/A^2 = 0.010290 Ha/bohr^2
EV_A2 = 0.529177210903**2 / 27.211386245988
m.K_TRAP = 1.0 * EV_A2
m.K_SPR = 0.26 * EV_A2
seed = int(sys.argv[1])
wd = pathlib.Path("seed%d" % seed)
# atom 0 = C (12.011 amu), atom 1 = H label (mass ~1.008 amu); dt = 0.5 fs = 20.67 a.u.; centroid tau = 20 fs = 826.8 a.u.
m.run_ipi(wd, nbeads=8, steps=int(sys.argv[2]), dt=20.67, tau=826.8, thermo="pile_l", classical="[0]", mode="nvt",
          stride=10, seed=seed, timeout=7200)
