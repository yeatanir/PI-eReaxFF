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
- [ ] Add thermostat interface and PILE-L.
- [ ] Add force provider interface.
- [ ] Add eReaxFF oracle adapter.

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
All implemented tests passed.
```

The free ring-polymer regression ran 20,000 exact propagation steps and reported
absolute free-energy drift of about `8.9e-13 eV`.
