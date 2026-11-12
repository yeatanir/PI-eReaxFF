# PI-eReaxFF build — progress report, 2026-11-11

Plan followed: `PI_eReaxFF_BUILD_PLAN.md` (Phases 0 and 1 started in parallel, as the plan's critical path prescribes).
Evidence labels: **[ran]** I ran it and saw the output; **[read]** read in source only; **[not done]**.

## 1. Summary

| Phase | State | One-line result |
|---|---|---|
| Layout | done | New folder `ereaxff_build/` with `forks/`, `project/`, `scratch/`; old `ereaxff_pimd_modular` untouched |
| 1.1 i-PI `classical_atoms` patch | implemented, tested, **not committed** | Classical atoms with ring-polymer electrons work in i-PI; conservation and analytic checks pass |
| 1.2 test potential | done | `coupled_harmonic` PES with exact references |
| 1.3 i-PI regression layout, L-P1, L-P2 | **not done** | |
| 1.4 fix F1 in the in-house engine | **not done** | |
| 0.1 rebuild original eReaxFF | partly done | Builds after a one-line fix; runs; gives NaN on the Jamil geometry |
| 0.2, 0.3, 0.4 | **not done** | 0.3 is where the NaN must be resolved |
| Phases 2–7 | not started | |

## 2. Folder structure created

```
ereaxff_build/
  README.md                      layout and rules
  forks/i-pi/                    i-PI clone (shallow, upstream 93162d2), branch `classical-atoms`
  forks/metal-reax/              reserved, not cloned yet
  project/                       git repo (2 commits)
    docs/STATUS.md, docs/PROGRESS_REPORT_2026-11-11.md
    inputs/{mixed_harmonic,jamil_classical,jamil_pimd}/   empty
    validation/phase0_reference/README.md
    validation/phase1_ipi_mixed/README.md
    validation/{phase3_engine,phase4_classical,phase5_pimd}/   empty
    reference_data/, tools/                                  empty
  scratch/phase0/                ereac_build/, co2_run/, jamil_vpar38_0/, jamil_vpar38_10/ (throwaway, git-ignored)
```

Outside this tree, one addition, made at your request: `Literature/.../Margaret/eReaxFF_src/build_2026-11-11/` (new subfolder; nothing existing overwritten).

## 3. Phase 1 — classical atoms in i-PI

### 3.1 What was changed (`forks/i-pi`, branch `classical-atoms`, 126 lines added to 3 files, 2 new files)

| File | Change |
|---|---|
| `ipi/engine/motion/dynamics.py` | New `classical_atoms` argument of `Dynamics`. `_setup_classical_atoms()` validates and builds the degree-of-freedom list. `DummyIntegrator.project_classical()` holds internal normal modes at zero. `pconstraints()` now also zeroes the internal-mode momenta of classical atoms. `NVTIntegrator.tstep()` books the removed thermostat noise into `ensemble.eens`. Initial projection at bind time. |
| `ipi/inputs/motion/dynamics.py` | New `<classical_atoms>` field inside `<dynamics>`, stored and fetched (so it survives checkpoints). |
| `ipi/engine/properties.py` | `temperature` corrected for the projected modes (otherwise a classical atom would report T/P). `kinetic_td` (primitive estimator) raises a clear error when classical atoms exist; use `kinetic_cv`. |
| `ipi/pes/coupled_harmonic.py` (new) | Test potential V = ½K\|X\|² + ½k\|x−X\|². |
| `ipi_tests/classical_atoms/test_classical_atoms.py` (new) | Runs i-PI and its Python driver as subprocesses and checks the results. |

Refused with an explicit error: ensembles other than nve/nvt, thermostats other than none/Langevin/PILE-L/PILE-G/SVR (GLE-type excluded), normal-mode mode other than `rpmd`, overlap with `bosons`, `open_paths` or `fixatoms`, ring-polymer contraction, out-of-range indices.

### 3.2 The scheme (why it is correct)

In i-PI's convention (ring frequency ω_P = P·k_B·T/ħ, thermostat temperature P·T, unscaled bead potential), a classical atom is a ring polymer whose internal modes k ≥ 1 are held at zero. Then the centroid equation reduces to dP/dt = (1/P) Σ_j f_j: physical mass, bead-averaged force, temperature T. Derivation in plan Appendix C.

Energy bookkeeping decision: momentum removed right after a force kick is **not** booked, because the constrained dynamics never applies that kick; thermostat noise that is added and then removed **is** booked in `eens`. (i-PI's existing `fixatoms` code books both. I did not change it.)

### 3.3 Test results [ran] (P = 8, 300 K, atomic units; classical atom 12 amu, quantum atom 1 amu)

| Check | Result | Reference | Verdict |
|---|---|---|---|
| Classical atom's beads identical | max difference < 1e-9 | 0 | pass |
| Ring bead variance about the classical atom (50 000 steps) | 0.0847 ± 0.0017 | 0.0859 (finite-P harmonic, exact) | within 1 σ |
| Ring centroid about the classical atom | 0.0464 ± 0.0016 | 0.0475 (k_B T/k) | within 1 σ |
| Classical atom ⟨X²⟩ per component | 0.0799 ± 0.0070 | 0.0950 (k_B T/K) | **2.2 σ low**; passes the 3 σ gate; weakest result |
| Reported temperature | 299.3 ± 1.0 K | 300 K | pass |
| Conserved quantity, NVE, 8000 steps, dt 4 | span 7e-4 k_B T, drift 1e-6 k_B T | 0 | pass |
| Conserved quantity, NVT PILE-L | span 1.2e-3 k_B T, drift 1.4e-4 k_B T | 0 | pass |
| Conserved quantity, NVT PILE-G | span 1.1e-3 k_B T, drift 2e-6 k_B T | 0 | pass |

Notes on strength: the ⟨X²⟩ check is limited by sampling (the classical oscillator has a ~9000 a.u. period; 400 000 a.u. gives ~40 periods). Longer runs (`IPI_CA_STEPS=400000`) are the remedy. The tests use one seed and one P. The ring-polymer checks are meaningful because the system is exactly solvable and the classical atom is coupled to the ring.

### 3.4 Not yet done in Phase 1

- **L-P1** (patch inactive equals unpatched i-PI bitwise) and **L-P2** (all atoms classical at P = 8 equals a P = 1 run): not run.
- Negative control: show that booking the kick energy would break conservation, to confirm the bookkeeping decision.
- Restart round-trip of `<classical_atoms>`, MTS (`nmts`), the `baoab` vs `obabo` split, and the Cayley propagator.
- Lennard-Jones or other anharmonic coupling test; thermostat-off equipartition per mode.
- Layout in i-PI's own regression-test format; upstream pull request.
- Nothing is committed in `forks/i-pi`.
- Plan task 1.4 (fix F1 in the in-house Fortran engine and add the production-path harmonic test).

## 4. Phase 0 — ground truth for eReaxFF

### 4.1 Rebuild of the original source [ran]

Source: `Margaret/eReaxFF_src` (`reac.f` hash 17ab34a5, `poten.f` 1f0f9274, `cbka.blk` 533eeff2).

- **Finding P0-1:** the current `reac.f` does not compile with the shipped full-size `cbka.blk`. The code that fills `znucel` (around `reac.f:4585`) needs `znucel` declared; only `cbka_pimd.blk` declares it. gfortran then reports "Syntax error in IF-clause" because `znucel(i1)=...` parses as a statement function. Consequence: the shipped full-size binaries (dated 2026-04-10) were built from an earlier source state, and the plan's finding V-14 (the source carries local April-2026 edits) is confirmed.
- **Fix used (scratch only):** add `     $       znucel(nsort),` after the `vkap(nvatym),` line of a copy of `cbka.blk`. Stock `makefile` then builds with 0 errors (gfortran 15.2, `-O3 -march=native`, arm64). Binary hash c9963253bb299d74.
- **Copied to the original folder:** `Margaret/eReaxFF_src/build_2026-11-11/` containing the patched `cbka.blk`, the `makefile`, the `ereac` binary and a README. Sources were not duplicated.

### 4.2 Runs [ran]

1. `test_co2` inputs (shipped): the driver stops silently after 0.01 s when only `geo` is present. The shipped folder is itself a failed run (its `run.log` shows a dyld error) and its `fort.73` holds an old BH3 result, so it cannot serve as a reference.
2. Discovered input conventions of this driver: **force field is read as `fort.4`, geometry as `fort.3`**.
3. Jamil 112 C / 113 El geometry with the V2 `control`, `ffield`, `geo` (`icharg 8`, fixed charges): the run **stops with "Energy not a number"**.
4. Same case with general parameter 38 (electron taper radius per the plan's finding V-13) changed from 0 to 10: **identical output and the same NaN**. So the NaN is not caused by the zero taper radius. The printed energy vector before the stop has a second entry of −336.372 and NaN in the vdW and Coulomb positions.

One transient `dyld ... syscall to map cache into shared region failed` occurred on a run launched through a Python subprocess; the same binary then ran normally from the shell. I did not investigate further.

### 4.3 What the NaN means for the plan

- Plan task 0.3 (2×2 table of periodic/non-periodic × parameter-38 values compared against the saved V2 energy terms) cannot be completed until the base code produces finite energies for this geometry.
- Hypotheses, none tested: (a) the V2 geometry or `control` needs entries the standalone driver handles differently from the V2 oracle wrapper (V2 used its own geometry reader); (b) the `El` radius or `ielec` handling; (c) a limit or array-size issue in the full-size `cbka.blk` (neem = 4000 vs ACKS2 sizing); (d) a charge issue despite `icharg 8`.
- Planned bisection: run a carbon-only subset; run with `El` atoms and electron terms removed; print charges and bond orders; compare with the saved `electron_trace.out` and `single_point_bond_orders.out` from the V2 run directory.

## 5. Evidence and limits

- Everything marked [ran] was run by me on this machine (macOS, arm64, Python 3.14.3, numpy 2.4.4, gfortran 15.2). Test runs use one random seed.
- i-PI tests need about 6 minutes in total at the default length.
- I did not run i-PI's own regression suite, so I have not shown that the patch leaves unrelated behaviour unchanged (that is L-P1).
- No change was made to `ereaxff_pimd_modular`, Metal-Reax, or any legacy source file. The only change outside `ereaxff_build/` is the new `build_2026-11-11/` subfolder.

## 6. Recommended next steps (in order)

1. Commit the i-PI work on the `classical-atoms` branch, then finish L-P1 and L-P2 and the negative control (about 1 day).
2. Run the longer ⟨X²⟩ check (400 000 steps, about 25 minutes) to confirm or remove the 2.2 σ deviation.
3. Resolve the NaN in the rebuilt `ereac` by bisection (task 0.3); decide with you whether the saved V2 outputs remain a reference.
4. Fix F1 in the in-house engine (task 1.4) and add the production-path harmonic test, as an independent check of the i-PI result.
5. Phase 2: wrap the JAX kernel as an i-PI driver and run the first eReaxFF path-integral NVE test; this is the earliest scientific signal about the instability seen earlier.
6. Clone Metal-Reax into `forks/metal-reax` and start Phase 3 once Phase 0 gives a trustworthy reference.

## 7. Questions for you

1. Commit the i-PI branch now, or after L-P1/L-P2?
2. If the standalone `ereac` cannot be made to reproduce the saved V2 energy terms, are you willing to treat the JAX port and the rebuilt `ereac` as the two references, as plan risk R1 proposes?
3. The plan's open owner questions (electron mass model, authoritative ffield for the Jamil benchmark, periodic vs central-cell target, public release of legacy-derived code) are still unanswered; the first two affect Phase 0 directly.
