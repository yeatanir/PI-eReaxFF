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
- [x] Implement `pimd/normal_modes.f90`.
- [x] Add normal-mode regression test.
- [ ] Add full core units/errors/RNG/state modules.
- [ ] Add free ring-polymer propagator.
- [ ] Add thermostat interface and PILE-L.
- [ ] Add force provider interface.
- [ ] Add eReaxFF oracle adapter.

## First Test Command

```sh
cd /Users/yeatanir/PhD/eReaxFF\ PIMD/ereaxff_pimd_modular
make test
```

