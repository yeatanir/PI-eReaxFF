# Lessons Learnt

This file records implementation problems, scientific risks, and fixes as the
modular PI-eReaxFF engine is built.  It should be updated at each phase before
the phase is considered complete.

## 2026-04-30: Initial Scaffold And Normal Modes

Problem:

- The legacy V1/V2/Python PIMD code duplicated normal-mode and PILE logic in
  several drivers.
- The old paths used `omega_P = sqrt(P) * kB*T/hbar`, while standard closed-ring
  PIMD and i-PI use `omega_P = P * kB*T/hbar`.
- The first `make test` tried to use the platform default `FC=f77`, which was
  not installed.

Resolution:

- Created a new standalone source tree instead of modifying V1/V2 in place.
- Centralized normal-mode construction and `omega_k` calculation in
  `src/pimd/normal_modes.f90`.
- Added a regression test that checks orthogonality, roundtrip, centroid
  consistency, and explicitly fails if the legacy `sqrt(P)` prefactor reappears.
- Updated the `Makefile` so the default compiler resolves to `gfortran` unless
  the caller explicitly overrides `FC`.

Validation:

- `make test` passed with normal-mode orthogonality error around `9e-16` and
  roundtrip error around `3e-14`.

## 2026-04-30: Phase 1/2 Next-Step Plan

Problem:

- Normal modes alone are not enough; before thermostats or eReaxFF coupling, the
  code needs tested unit conversion, spring-energy, and exact free-ring
  propagation machinery.

Resolution:

- Proceed with a narrow Phase 1/2 increment: units, errors, RNG scaffold,
  ring-polymer initialization, spring energy, exact free propagator, and
  conservation tests.

Validation:

- `make test` now runs:
  - constants/unit/RNG smoke test;
  - normal-mode transform and standard-frequency regression;
  - exact free-ring propagation conservation test.
- The free-ring test propagated 20,000 steps with absolute energy drift below
  `1e-12 eV`.

Problem found during implementation:

- The first constants test used a reference value for
  `amu*Angstrom^2/fs^2 -> eV` derived from a slightly different AMU constant
  than the one encoded in `src/core/units.f90`.

Resolution:

- Kept the explicit CODATA AMU value in the unit module and updated the test
  reference to match that exact constant set.  This keeps the test strict while
  avoiding a false failure from mixed reference constants.
