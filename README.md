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

The current implementation now spans the validated modular foundation and the
first working mixed RPMD/PIMD eReaxFF dynamics boundary:

- repository scaffold and architecture documentation;
- core numeric kinds and physical constants;
- real orthonormal closed-ring normal-mode transform;
- standard i-PI ring-polymer frequency calculation;
- bead-to-normal-mode and normal-mode-to-bead transforms;
- regression test for orthogonality, roundtrip, centroid consistency, and the
  standard `omega_P` prefactor;
- eReaxFF oracle adapter boundary with repeatability and NaN rejection;
- force-manager `P=1` reproduction of the validated V2 oracle;
- mixed classical/PIMD grouping with explicit thermostat selection;
- mixed BAOAB restart and compatibility output handling;
- multiple PIMD particles in one mixed run when they share a common bead
  count;
- explicit pseudo-mass and restart state for each quantum particle;
- persistent eReaxFF oracle worker and bounded worker-pool backends;
- Jamil 112C/113El graphene/electron replication inputs;
- P=16 eReaxFF one-step runs through the persistent worker pool;
- physical versus extended ring-polymer temperature diagnostics;
- centroid/internal electron normal-mode kinetic partition diagnostics;
- optional one-time initial PIMD normal-mode velocity rescaling for small
  stochastic systems;
- force-off mixed BAOAB diagnostics for separating thermostat/free-ring
  behavior from eReaxFF force feedback.

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
test_03_harmonic_quantum_variance: PASS
test_04_thermostat_ou_statistics: PASS
test_05_pile_mode_temperatures: PASS
test_05_pile_g_stub: PASS
test_06_force_scaling: PASS
test_07_charge_memory: PASS
test_08_reaxff_oracle_smoke: PASS
test_09_mixed_classical_pimd_config: PASS
test_10_harmonic_driver_outputs: PASS
test_11_structural_stress: PASS
test_11_invalid_group_overlap: PASS
test_12_oracle_force_manager_p1: PASS
test_13_oracle_charge_policy_modes: PASS
test_14_mixed_oracle_manager_fixture: PASS
test_15_mixed_oracle_manager_matrix: PASS
test_16_mixed_oracle_manager_restart_output: PASS
test_17_ereaxff_mixed_driver_smoke: PASS
test_18_ereaxff_mixed_baoab_matrix: PASS
test_19_ereaxff_mixed_baoab_restart: PASS
test_20_ereaxff_mixed_baoab_policy_matrix: PASS
test_21_ereaxff_mixed_baoab_config_smoke: PASS
test_22_ereaxff_mixed_baoab_config_matrix: PASS
test_23_ereaxff_mixed_baoab_geometry_matrix: PASS
test_24_ereaxff_mixed_baoab_user_geo_matrix: PASS
test_25_phase6_like_temperature_partition: PASS
test_26_crystx_periodic_loader: PASS
test_27_multi_quantum_restart: PASS
test_28_paired_graphene_60el_60c_smoke: PASS
test_29_oracle_dispatch_modes: PASS
test_30_persistent_oracle_worker: PASS
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

1. Extend the eReaxFF-backed mixed driver from smoke coverage into broader
   modular dynamics and restart coverage.
2. Add ReaxFF compatibility projections on top of the mixed driver.

## Status

This repository is still pre-production research software, but the current
mixed RPMD/PIMD eReaxFF plumbing is now working for the Jamil graphene/electron
target.

The normal-mode, free-ring propagation, thermostat, dummy/harmonic force,
mixed group assignment, adapter-mediated V2 oracle coverage, P=1 force-manager
integration, centroid charge policy, mixed oracle-manager fixture/matrix,
restart/output, periodic geometry loading, and persistent oracle worker-pool
foundations are tested.

The mixed eReaxFF BAOAB executable now reads the YAML-subset config directly,
loads user-supplied legacy `geo` files with `CRYSTX` cell data, runs P=16
Jamil 112C/113El one-step eReaxFF dynamics through the persistent worker pool,
and writes ReaxFF-compatible projections plus PIMD-native diagnostics.  The
electron temperature diagnostics distinguish the physical aggregate temperature
from the extended ring-polymer temperature, and
`electron_mode_partition.out` splits kinetic energy into centroid and internal
normal-mode contributions.

For small stochastic systems, `simulation.rescale_initial_pimd_velocities:
true` can be used to remove unlucky hot initial Maxwell draws.  On the Jamil
P=16, `dt=0.05 fs` one-step oracle check, this reduced the physical electron
temperature from about `489 K` to about `369 K`; the force-off validation starts
exactly at `300 K` and averages near `300 K` over 10000 steps.

The current Jamil production-style config uses `mass_amu: 1.0` for `El`, which
is a fictitious eReaxFF electron mass convention rather than the physical
electron mass.  Force-off mass checks show that the kinetic temperature
diagnostic is mass-independent for the same random draw, but bead spread is not:
using the physical electron mass produces much larger bead motion and is likely
not stable unless the model is parameterized for it.

The main remaining blocker is speed.  A P=16 one-step Jamil run currently makes
80 eReaxFF oracle calls and takes about 7.5 minutes on the tested machine, even
with four persistent workers.  The next engineering target is to bring the
oracle path much closer to native ReaxFF throughput so this can become usable
for longer MD trajectories.

The fast-backend work has started behind a gated interface:

- the exact legacy oracle remains the default and validated backend;
- `force_provider.type: fast_jax_oracle` is parsed but intentionally stops at
  driver startup until parity validation exists;
- `projects/jamil_replication/fast_backend/build_oracle_dataset.py` can extract
  JSONL reference cases from exact Jamil restart or adapter-work output for
  JAX/native backend development;
- `projects/jamil_replication/fast_backend/jax_backend.py` is a fail-closed
  optional entrypoint for JAX availability, dataset self-tests, and labeled
  oracle readiness checks;
- the first JAX charge kernel reproduces the current Jamil `icharg=8`
  fixed-`MOLCHARGE` charge path in batched form;
- the force-provider base now has a batched result type so future backends can
  evaluate all beads in one compiled/vectorized request.
