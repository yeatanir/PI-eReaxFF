# Phase 2 — JAX kernel as an i-PI driver (interim engine)

Files:
- `../../tools/ipi_ereaxff_jax_driver.py` — i-PI Python PES driver (`-m custom -P`), wraps the old prototype's JAX kernel (read-only import from `ereaxff_pimd_modular/.../fast_backend`). Atom types and fixed charges (MOLCHARGE) come from the geo file; units converted a.u. <-> A/eV. Handles i-PI batched requests (all P beads in one call).
- `../../inputs/jamil_pimd/run_ipi_jax.py` — builds the i-PI input (112 C classical via `<classical_atoms>`, 113 El as P-bead rings, mass 1 amu) and runs server + driver.
- `validate_vs_ereac.py` — compares the driver with the exact rebuilt `ereac` export.

Physics: central-cell (non-periodic) fixed-charge eReaxFF of the Jamil target, i.e. a **finite flake**; the rebuilt `ereac` (non-periodic, electron taper 10 A) is the exact reference for it (see phase 0). Stability work only, no science claims. The ffield is the repo's (`np_export/ffield`, parameter 38 = 0 as written!).

Important: the repo ffield has `0.0` in general parameter 38; the JAX kernel ignores that slot and effectively uses a 10 A electron taper (it matches V2). For `ereac` the taper must be set explicitly to 10 (done in the scratch copy). Which value is right physically is an open owner question.

## Results

| Check | Result |
|---|---|
| Driver vs exact `ereac` export (225 atoms, step 0) | max force difference **6.3e-7 eV/A** (max force 12.34 eV/A, rms difference 9.8e-8); energy difference 8.95e-6 eV (1.7e-8 relative) |
| Chain smoke test (i-PI + classical atoms + batched driver, P = 2, 2 steps) | runs; conserved quantity changes by 7e-4 eV in the first step (see probe below) |

Run the validation: `PYTHONPATH=<forks/i-pi> <venv-with-jax>/bin/python validate_vs_ereac.py`.

## F2 probe (is eReaxFF PIMD stable with a correct integrator?)

Probe: P = 4, NVE, dt = 0.1 fs, 400 steps planned, thermal start at 300 K, **unrelaxed** geometry (the raw Jamil flake with electrons at their generated positions). About 20 s per step on this machine (the JAX kernel uses all cores).

Partial result at step 174 (run still going when written):
- Conserved quantity (printed per bead): -519.2895 -> -519.4960 eV, span 0.475 eV (0.09 % of its value) while the potential energy fell from -541.6 to -658.7 eV (117 eV) and the kinetic energy rose from 22 to 127 eV (T 304 -> 1205 K).
- Reading: the integrator conserves the extended energy to about 0.4 % of the energy released, so the strong heating is **not** an integrator defect: the starting geometry is far from equilibrium and releases ~115 eV. This suggests the heating seen in the old engine's runs (carbon temperature 1156-2139 K) may also have been the unrelaxed start, not only integrator defects. To be tested properly: relax the geometry first (or thermostat from the start), then repeat NVE.
- Caveats: one run, one P, 17 fs; the conserved-quantity span is only meaningful relative to the energy scale of the motion, so a relaxed start is needed for a clean drift number.
