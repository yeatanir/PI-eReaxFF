"""Compare the JAX-kernel driver with the exact rebuilt ereac (full-precision export) on the V2 reference geometry.

Run with the JAX environment:  <venv>/bin/python validate_vs_ereac.py
(needs i-PI importable: PYTHONPATH=<forks/i-pi>)
"""
import os
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve()
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "tools"))
REF = PROJ / "reference_data/phase0_jamil_225atom/np_export"
os.environ["EREAXFF_FFIELD"] = str(REF / "ffield")
os.environ["EREAXFF_GEO"] = str(REF / "geo")
import ipi_ereaxff_jax_driver as d  # noqa: E402

KCAL_EV = 0.043364  # the value the old tools used; exact ereac gradient is in kcal/mol/A

drv = d.EReaxFF_JAX_driver()
symbols, pos_A, q, _ = d.read_geo(REF / "geo")
lines = (REF / "ereac_export.out").read_text().splitlines()
grad, estrc = [], None
step = None
for l in lines:
    w = l.split()
    if w[0] == "STEP":
        step = int(w[1])
    elif step == 0 and w[0] == "ESTRC":
        estrc = float(w[1])
    elif step == 0 and w[0] == "GRAD":
        grad.append([float(x) for x in w[2:5]])
grad = np.array(grad)
e_ref_eV = estrc * KCAL_EV
f_ref = -grad * KCAL_EV  # eV/A
e, f = drv.evaluate_positions([pos_A.flatten() / d.BOHR_A])
e_eV, f_eV = e[0] * d.HARTREE_EV, f[0] * d.HARTREE_EV / d.BOHR_A
print("energy  ereac %.10f eV   JAX %.10f eV   diff %.3e eV" % (e_ref_eV, e_eV, e_eV - e_ref_eV))
err = np.abs(f_eV - f_ref)
print("forces  max |dF| = %.3e eV/A (max |F| = %.4f)   rms dF = %.3e   worst atom %d" % (
    err.max(), np.abs(f_ref).max(), np.sqrt((err ** 2).mean()), int(np.unravel_index(err.argmax(), err.shape)[0]) + 1))
