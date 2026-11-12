# Phase 1 — i-PI mixed classical / ring-polymer dynamics

Code: `../../../forks/i-pi`, branch `classical-atoms`, commits `1b8f05e` and `1a47729`.

Changed files:
- `ipi/engine/motion/dynamics.py` — `classical_atoms` option, validation, `project_classical`, projection in `pconstraints`, energy booking after the thermostat step.
- `ipi/inputs/motion/dynamics.py` — `<classical_atoms>` input field (stored in checkpoints).
- `ipi/engine/properties.py` — `temperature` corrected for the projected modes; `kinetic_td` refuses classical atoms.
- `ipi/pes/coupled_harmonic.py` — test potential (new).
- `ipi_tests/classical_atoms/test_classical_atoms.py` — tests (new).

Algorithm: Appendix C of `PI_eReaxFF_BUILD_PLAN.md`. A classical atom keeps q_k = p_k = 0 for k >= 1.
Energy bookkeeping: momentum removed right after a force kick is **not** booked (it was never part of the
constrained dynamics); noise added by the thermostat and then removed **is** booked in `ensemble.eens`.
(The existing `fixatoms` code books both, which would make `econs` drift for forces acting on frozen modes; not changed.)

## Results (2026-11-11, P = 8, 300 K, atomic units)

| Test | Plan rung | Result |
|---|---|---|
| Beads of the classical atom identical | L-P3 | < 1e-9 |
| Ring bead variance about the classical atom, 4 seeds × 25 000 steps | L-P3 | 0.08586 ± 0.00072 vs exact 0.08591 (−0.1 σ) |
| Ring centroid about the classical atom | L-P3 | 0.04671 ± 0.00044 vs exact 0.04750 (−1.8 σ, −1.7 %) |
| Classical atom ⟨X²⟩ per component (stiffer trap, 4 seeds) | L-P3 | 0.00990 ± 0.00042 vs exact 0.00950 (+1.0 σ). The earlier single-seed 2.2 σ deficit was a test-design artefact. |
| Reported temperature | L-P5 | 301.2 ± 0.7 K vs 300 K |
| Conserved quantity, NVE, 8000 steps | L-P4 | span 7.3e-4 kT, drift 7e-8 kT |
| Conserved quantity, NVT PILE-L / PILE-G | L-P4 | span 1.2e-3 kT, drift 1.4e-4 / 3.6e-5 kT |
| Patch inactive vs unpatched i-PI 93162d2, 400 steps | L-P1 | **byte-identical** (energies, bead positions, centroid, velocities) |
| All atoms classical at P = 8 vs P = 1 NVE, 1500 steps | L-P2 | positions equal to output resolution (1e-6 a.u.) over 2 a.u. of motion; potential energy equal to 7e-10 Ha |
| Negative control: also booking the kick momentum | — | conserved quantity drifts 0.031 kT per 8000 steps (vs 7e-8 kT); the drift test would catch it |

Guard: `project_classical` raises if the internal-mode positions of a classical atom are not already zero.
The centroid check (−1.8 σ) is the one result still to watch; it passes the 3 σ / 3 % gate.

**Not done:** restart round-trip of `<classical_atoms>`, MTS (`nmts`), `baoab` vs `obabo`, Cayley propagator, an anharmonic coupling test, i-PI's own regression-test layout, upstream pull request.

Run: `cd forks/i-pi && IPI_UPSTREAM_ROOT=<unpatched checkout> python3 -m pytest -q -s ipi_tests/classical_atoms/test_classical_atoms.py` (about 10 minutes; `IPI_CA_STEPS` sets the statistical run length).

## Loose ends (workstream D, 2026-11-11/12) — all run with OMP_NUM_THREADS=1

i-PI branch `classical-atoms`, commits `2db34b8`, (baoab/MTS), `1c13466`, `588ce4e`. All tests in `ipi_tests/classical_atoms/test_classical_atoms.py`.

| Task | Test | Result |
|---|---|---|
| 1 Restart | checkpoint at step 100, restart to 200, compare with uninterrupted run | `<classical_atoms>` is in the checkpoint; outputs agree to **5e-9 relative** (positions to the 1e-6 output resolution). Not bitwise: the checkpoint stores floats with limited precision. **Unpatched i-PI shows the same 4e-9**, measured as control. Classical beads stay identical. |
| 2 Splittings | BAOAB NVT PILE-L / PILE-G, 8000 steps; BAOAB coupled-harmonic statistics (3 seeds) | span 1.2e-3 / 1.3e-3 kT; ⟨X²⟩ +0.4σ, bead variance −0.2σ, centroid −0.3σ |
| 3 MTS | two force components (slow trap, fast spring), `nmts=[1,4]`, NVE and NVT, obabo and baoab | span ≤ 1.2e-3 kT, drift ≤ 1.2e-4 kT. Works; no refusal needed. |
| 4 Propagators | exact / cayley / bab × fft / matrix, NVT 4000 steps | span 5e-3 kT, drift ≤ 1.8e-4 kT; classical beads identical. All six pass; no refusal needed. |
| 5 Anharmonic | quartic coupling g=0.1 (PES forces FD-checked to 1e-12); reference = exact finite-P transfer-matrix quadrature (validated against the harmonic analytic result to 1e-4); harmonic would be 2× larger, so the test is sensitive | dt=6, 3 seeds×25000: ⟨X²⟩ +0.3σ, centroid +0.4σ, bead variance **+2.9σ (+1.1 %)**. Rerun dt=3, 3 seeds×50000: bead +1.3σ (+1.2 %), centroid +0.6σ (+1.4 %). Residual ≈ +1 % is consistent with time-step/sampling error; **not resolved beyond that**. |
| 6 Regression layout | `ipi_tests/regression_tests/tests/CLASSICAL_ATOMS/nvt_pilel_{obabo,baoab}/coupled_harmonic_python` with references, driver registered in `test_tools.py` | Reference outputs regenerate exactly (all 5 files match, both splittings) when rerun by hand. **i-PI's own `test_run.py` harness did not complete** (hung with no output for >5 min, both through pytest and as a script; not diagnosed), so the layout has not been accepted by the official runner. |

Harness note: the test file now takes `splitting`, `nmts`, extra drivers, normal-mode options; a restart file records the socket address, which `launch` reads back.
Draft PR text (not submitted): `UPSTREAM_PR_DRAFT.md`.

Not done: official regression runner, fixing/diagnosing its hang, upstream submission, NPT/GLE/bosons with classical atoms (refused by design).
