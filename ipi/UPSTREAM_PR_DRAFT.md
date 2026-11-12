# DRAFT — not submitted. Upstream pull-request description for `classical_atoms`

**Title:** Add `<classical_atoms>` to `<dynamics>`: classical particles inside a ring-polymer simulation

## Motivation
All atoms of an i-PI system share one bead number, so a system with a few quantum particles (e.g. explicit electrons or protons) in a classical environment must either ring-polymerise everything or use an ad-hoc setup. This adds an option that keeps selected atoms classical while the rest are ring polymers, in the same run, with the usual thermostats and estimators.

## What it does
`<classical_atoms> [0, 5, ...] </classical_atoms>` inside `<dynamics>` (atom indices, 0-based). For these atoms the internal normal modes of the ring polymer are held at zero (q_k = p_k = 0, k >= 1), so all their beads coincide.
In i-PI's convention (ω_P = P k_B T/ħ, thermostat temperature P·T, unscaled bead potential) such an atom is exactly a classical particle with its physical mass, at temperature T, moving under the bead-averaged force: the centroid mode obeys dP/dt = (1/P) Σ_j f_j.

## Changes
- `ipi/engine/motion/dynamics.py`: new argument, `_setup_classical_atoms` (validation), `DummyIntegrator.project_classical`, projection in `pconstraints`, booking of removed thermostat noise in `ensemble.eens` (`NVTIntegrator.tstep`). A guard raises if the projected positions are ever non-zero.
- `ipi/inputs/motion/dynamics.py`: the input field (stored in checkpoints).
- `ipi/engine/properties.py`: `temperature` accounts for the projected modes; `kinetic_td` refuses classical atoms (use `kinetic_cv`, which is already correct).
- `ipi/pes/coupled_harmonic.py`: test potential (classical trap + quadratic/quartic coupling to a ring).
- Tests: `ipi_tests/classical_atoms/` (statistical, conservation, equivalence) and `ipi_tests/regression_tests/tests/CLASSICAL_ATOMS/`.

## Energy bookkeeping
Momentum removed right after a force kick is **not** booked (the constrained dynamics never applies that kick); thermostat noise that is added and then removed **is** booked in `eens`. A negative control that also books the kick momentum makes the conserved quantity drift (0.031 kT per 8000 steps vs 7e-8 kT).

## Supported / refused (explicit errors)
Supported: `nve` and `nvt` ensembles; thermostats none, Langevin, PILE-L, PILE-G, SVR; `obabo` and `baoab` splittings; multiple time stepping; exact, Cayley and `bab` free-ring propagators; fft and matrix transforms; restart.
Refused: other ensembles (npt, nst, nvt-cc, sc, driven), GLE-type thermostats, normal-mode mode other than `rpmd`, overlap with `bosons`/`open_paths`/`fixatoms`, ring-polymer contraction, out-of-range indices.

## Evidence (P = 8 unless stated, atomic units, 300 K)
See `validation/phase1_ipi_mixed/README.md` for numbers. Summary: patch inactive = byte-identical to upstream; all atoms classical at P = 8 reproduces P = 1; coupled harmonic and quartic systems match exact finite-P references within statistical error; extended energy conserved to ~1e-3 kT in NVE/NVT for both splittings, MTS, all propagators; restart matches an uninterrupted run to checkpoint float precision (5e-9 relative, same as unpatched).

## Known limits / reviewer questions
- Only the quantum-classical *separation* is implemented; no estimator for classical-atom kinetic energy beyond `kinetic_md`/`kinetic_cv`.
- Fixed `fixatoms` code books both kick and thermostat momentum; this patch deliberately does not change it.
- Is a per-atom option under `<dynamics>` the preferred place, or should it live under `<normal_modes>`?
- Open: exchange/bosons with classical atoms, NPT, GLE thermostats.
