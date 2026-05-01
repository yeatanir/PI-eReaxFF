# Project Status

Generated: 2026-04-30

## Inputs Provided

- eReaxFF V1: `/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1`
- eReaxFF V2: `/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2`
- i-PI local source: `/Users/yeatanir/PhD/PIMD/i-pi`
- Architecture report:
  `/Users/yeatanir/PhD/eReaxFF PIMD/i_pi_vs_ereaxff_pimd_architecture_report.txt`

## Initial Findings

- V2 contains a useful force-only boundary in `src/ereaxff_force_oracle.f`.
- V2 phase4/phase5/phase6 and harmonic drivers duplicate normal-mode,
  thermostat, propagation, diagnostics, and output logic in monolithic drivers.
- V1 debug PIMD paths contain comparable duplicated normal-mode routines.
- The report and direct search both identify the legacy `sqrt(P)` ring-polymer
  prefactor in V1/V2 PIMD paths.
- The new code will default to the standard i-PI closed-ring convention.

## Checklist

- [x] Create modular directory tree.
- [x] Add README, DESIGN, VALIDATION, and PROJECT_STATUS.
- [x] Add docs/LEGACY_AUDIT.md.
- [x] Implement `core/kinds.f90`.
- [x] Implement `core/constants.f90`.
- [x] Implement `core/units.f90`.
- [x] Implement `core/errors.f90`.
- [x] Implement `core/rng.f90`.
- [x] Implement `pimd/normal_modes.f90`.
- [x] Implement `pimd/ring_polymer.f90`.
- [x] Implement `pimd/spring_energy.f90`.
- [x] Implement `pimd/free_ring_propagator.f90`.
- [x] Add constants/unit/RNG regression test.
- [x] Add normal-mode regression test.
- [x] Add free ring-polymer energy conservation regression test.
- [ ] Add system state and particle-group modules.
- [x] Add thermostat base helpers.
- [x] Add `none` thermostat.
- [x] Add Cartesian Langevin thermostat with heat accounting.
- [x] Add PILE-L thermostat with mode-dependent damping.
- [x] Add PILE-G stub that fails clearly.
- [x] Add thermostat registry name checks.
- [x] Add Langevin OU variance/heat regression test.
- [x] Add PILE-L mode-temperature/heat regression test.
- [x] Add PILE-G negative regression test.
- [x] Add force provider result/interface scaffolding.
- [x] Add dummy force provider.
- [x] Add harmonic force provider.
- [x] Add bead force manager with `1/P` physical force scaling.
- [x] Add charge coupling policy names and bead-local charge memory helper.
- [x] Add force provider/scaling regression test.
- [x] Add charge policy/memory regression test.
- [x] Add eReaxFF oracle adapter failure stub.
- [x] Add system state module.
- [x] Add selector module for `all`, `element == X`, `id == N`, `id != N`, and `id in [...]`.
- [x] Add particle group assignment module.
- [x] Add mixed classical/PIMD/frozen group regression test.
- [x] Add harmonic finite-P quantum variance utility.
- [x] Add minimal harmonic BAOAB PIMD step.
- [x] Add harmonic quantum variance regression for `P = 1, 2, 4, 8, 16, 32`.
- [x] Add stochastic BAOAB+PILE-L harmonic variance smoke test.
- [ ] Connect eReaxFF oracle adapter to legacy V2 oracle.

## First Test Command

```sh
cd /Users/yeatanir/PhD/eReaxFF\ PIMD/ereaxff_pimd_modular
make test
```

## Current Test Summary

Last local run:

```text
test_00_build: PASS
test_01_normal_modes: PASS
test_02_free_ring_polymer: PASS
test_03_harmonic_quantum_variance: PASS
test_04_thermostat_ou_statistics: PASS
test_05_pile_mode_temperatures: PASS
test_05_pile_g_stub: PASS
test_06_force_scaling: PASS
test_07_charge_memory: PASS
test_09_mixed_classical_pimd_config: PASS
All implemented tests passed.
```

The free ring-polymer regression ran 20,000 exact propagation steps and reported
absolute free-energy drift of about `8.9e-13 eV`.

The Langevin OU test sampled velocity variance within about `0.7%` of the
Maxwell target.  The PILE-L internal mode test reported maximum internal mode
temperature relative error of about `2.3%`.

The force-scaling test verifies harmonic `F = -kx`, dummy zero forces,
bead-averaged PIMD potential, and centralized `1/P` physical force scaling.  The
charge-memory test verifies valid policy names and independent bead-local charge
storage.

The mixed group test verifies config-style selection of classical carbons, one
PIMD electron with 16 beads, and frozen background electrons without source-code
special cases.

The harmonic quantum variance test verifies finite-P convergence toward the
analytic quantum oscillator variance, verifies that the legacy `sqrt(P)`
prefactor is much worse at `P=16`, and samples a `P=16` BAOAB+PILE-L trajectory
within about `2.1%` of the finite-P reference.
