# PI-eReaxFF

PI-eReaxFF is a new modular path-integral molecular dynamics engine for
eReaxFF-style force and charge models.  The goal is to build a robust standalone
PIMD codebase with an i-PI-like internal architecture, rather than another
hard-coded research driver for a single phase or system.

The core design principle is separation of responsibility:

- the eReaxFF force/charge oracle evaluates energy, forces, charges, and
  diagnostics;
- the dynamics engine owns state, PIMD propagation, thermostats, integration,
  restart, and output;
- particle grouping, bead counts, charge-coupling policy, and thermostat choice
  are input-level decisions, not source edits.

Legacy eReaxFF V1/V2 code is treated as source material and an eventual oracle
dependency.  It is not modified in place.

## Why This Exists

Existing eReaxFF-PIMD prototypes contain useful science and force-oracle work,
but too much algorithmic machinery is embedded in monolithic drivers:

- normal-mode transforms;
- ring-polymer frequencies;
- BAOAB propagation;
- PILE/Langevin thermostat logic;
- bead force orchestration;
- ACKS2 charge-memory handling;
- diagnostics and ReaxFF-compatible output.

This repository centralizes those invariants into reusable modules so that new
systems and approximations are controlled through configuration and tests.

## First Scientific Guardrail

This code uses the standard closed-ring PIMD frequency convention by default:

```text
omega_P = P * kB * T / hbar
omega_k = 2 * omega_P * sin(pi * k / P)
```

The older eReaxFF PIMD prototypes used `sqrt(P) * kB*T/hbar` in several paths.
That convention makes internal ring-polymer frequencies too small by
`sqrt(P)` and spring stiffness too weak by `P` relative to standard i-PI-style
closed-ring PIMD.  For `P = 16`, frequencies are 4x too small and spring
stiffness is 16x too weak.  The normal-mode tests in this tree explicitly fail
if the
`sqrt(P)` prefactor reappears.

## Planned Architecture

```text
src/
  core/          numeric kinds, constants, units, errors, RNG, state
  force/         force-provider interface, harmonic/dummy providers, eReaxFF adapter
  pimd/          normal modes, ring polymer propagation, springs, estimators
  thermostat/    none, Langevin, PILE-L, PILE-G/SVR interfaces
  integrator/    BAOAB/OBABO drivers using thermostat and force interfaces
  output/        PIMD-native diagnostics and ReaxFF compatibility projections
  drivers/       thin executable entry points
  ipi_socket/    future i-PI socket adapter plan and stub
```

The future i-PI socket path will reuse the same eReaxFF oracle adapter as the
standalone engine.  It should not become a second independent force path.

## Current Scope

The current implementation is the first foundation deliverable:

- repository scaffold and architecture documentation;
- core numeric kinds and physical constants;
- real orthonormal closed-ring normal-mode transform;
- standard i-PI ring-polymer frequency calculation;
- bead-to-normal-mode and normal-mode-to-bead transforms;
- regression test for orthogonality, roundtrip, centroid consistency, and the
  standard `omega_P` prefactor.

eReaxFF oracle coupling, thermostats, mixed classical/PIMD grouping, restart,
and compatibility outputs are planned but intentionally not implemented before
the normal-mode foundation is tested.

## Build And Test

Requirements:

- GNU Make
- `gfortran` with Fortran 2018 support

Run:

```sh
make test
```

or:

```sh
tests/run_all.sh
```

Expected current summary:

```text
test_00_build: PASS
test_01_normal_modes: PASS
test_02_free_ring_polymer: PASS
test_04_thermostat_ou_statistics: PASS
test_05_pile_mode_temperatures: PASS
test_05_pile_g_stub: PASS
All implemented tests passed.
```

## Documentation

- [DESIGN.md](DESIGN.md): ownership boundaries and scientific conventions.
- [VALIDATION.md](VALIDATION.md): staged validation ladder.
- [PROJECT_STATUS.md](PROJECT_STATUS.md): current checklist and inspected
  legacy paths.
- [docs/LEGACY_AUDIT.md](docs/LEGACY_AUDIT.md): migration notes from V1/V2 and
  i-PI comparison.
- [docs/INPUT_SPEC.md](docs/INPUT_SPEC.md): target configuration schema.

## Development Roadmap

1. Complete system-state and particle-group modules.
2. Add harmonic/dummy force providers and bead force manager.
3. Validate harmonic quantum variance versus analytic reference.
4. Wrap the legacy V2 eReaxFF force oracle without editing legacy code in place.
5. Add mixed classical/PIMD group selection and charge-coupling policies.
6. Add PIMD-native outputs, restart, and ReaxFF compatibility projections.

## Status

This repository is pre-production research software.  The normal-mode,
free-ring propagation, and initial thermostat foundation is tested; production
eReaxFF coupling is intentionally still behind the validation gates described
above.
