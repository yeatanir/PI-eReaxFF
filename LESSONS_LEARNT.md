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

## 2026-04-30: Thermostat Layer

Problem:

- Thermostats need to be first-class modules, but adding a full polymorphic
  hierarchy before an integrator exists would create unnecessary surface area.
- PILE-L and PILE-G must not be conflated.  PILE-G requires global stochastic
  velocity rescaling on the centroid, which is not implemented yet.
- The first `none` thermostat implementation produced compiler warnings because
  no-op velocity arguments were not referenced.

Resolution:

- Implemented concrete thermostat modules with a minimal shared convention:
  `init`, `step`, and accumulated `heat_eV`.
- Added `thermostat_base.f90` for kinetic-energy and Maxwell velocity-scale
  helpers so heat accounting and OU variance use one unit convention.
- Implemented Cartesian Langevin with exact OU velocity updates.
- Implemented PILE-L with centroid friction from input and internal mode
  friction `gamma_k = 2*lambda*omega_k`.
- Added a PILE-G stub that stops with:
  `PILE-G requested but not implemented; use PILE-L or implement global centroid SVR.`
- Added a negative test that confirms PILE-G fails with that message.
- Removed no-op warnings by explicitly referencing the velocity arrays in the
  `none` thermostat without changing values or heat.

Validation:

- `test_04_thermostat_ou_statistics` checks Langevin velocity variance against
  `kBT/m`, verifies nonzero heat exchange, verifies the `none` thermostat, and
  checks registry behavior.
- `test_05_pile_mode_temperatures` checks PILE-L internal mode temperatures and
  heat accounting.
- `test_05_pile_g_stub` confirms PILE-G does not silently run as PILE-L.
