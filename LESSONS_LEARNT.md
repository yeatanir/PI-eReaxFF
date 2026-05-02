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

## 2026-04-30: Force Provider And Bead Force Manager

Problem:

- Physical PIMD force scaling must not be scattered through future drivers.
  The physical potential convention is `V_PIMD = (1/P) sum_s V(R_s)`, so
  physical forces from each bead must be scaled by `1/P`.
- eReaxFF/ACKS2 is not ready to connect yet, but the force-manager contract must
  be testable before that coupling.
- The first force-manager compile produced `gfortran` warnings around local
  derived-type temporaries with allocatable components.

Resolution:

- Added a shared `force_result` container and concrete dummy/harmonic force
  providers.
- Added `bead_force_manager.f90` as the only current location that applies
  physical bead-force scaling.
- Added `charge_coupling_policy.f90` with explicit policy names and bead-local
  charge memory helpers.
- Kept `ereaxff_oracle_adapter.f90` as an explicit failure stub until the
  harmonic/dummy path is validated.
- Initialized local `force_result` temporaries before provider calls to avoid
  compiler warnings with allocatable components.

Validation:

- `test_06_force_scaling` checks harmonic `F = -kx`, dummy zero force,
  bead-averaged PIMD potential, and `1/P` scaled bead forces.
- `test_07_charge_memory` checks policy validation, bead-local memory detection,
  per-bead charge sums, max charge difference across beads, and independence of
  bead-local charge storage.

## 2026-04-30: System State And Particle Groups

Problem:

- The engine needs to support mixed classical/PIMD/frozen atoms and electrons
  through input-like group definitions, not source edits.
- A full YAML parser and expression language would be premature before the
  integrator exists.
- The first selector implementation used `.and.` with impure helper calls,
  which triggered a compiler warning because Fortran may short-circuit
  evaluation.

Resolution:

- Added `system_state.f90` for ids, elements, masses, positions, and velocities.
- Added a deliberately small selector subset: `all`, `element == X`, `id == N`,
  `id != N`, `id in [...]`, and one `and`.
- Added `particle_groups.f90` to validate group specs and assign per-particle
  dynamics, bead count, mobility, thermostat, and temperature.
- Rewrote compound selector evaluation with explicit temporaries so compiler
  output stays warning-free.

Validation:

- `test_09_mixed_classical_pimd_config` builds a C/El toy system and assigns
  mobile classical carbons, one 16-bead PIMD electron, and frozen background
  electrons using selector strings.
- The test also checks invalid selector rejection and `id in [...]` behavior.

## 2026-04-30: Harmonic Quantum Variance And BAOAB Path

Problem:

- The harmonic quantum variance test is the first assembled PIMD dynamics gate,
  so it must catch convention errors, not merely run.
- Reinitializing one `ring_polymer_state` for several bead counts failed because
  `init_ring_polymer` did not deallocate already allocated arrays.
- The first finite-P variance formula mixed bead-averaged force scaling with the
  wrong canonical temperature/spring convention, which made the convergence
  trend fail.

Resolution:

- Fixed `init_ring_polymer` so it can be safely reused for scans over bead
  count.
- Added `harmonic_quantum_variance.f90` for analytic quantum and finite-P
  oscillator variance checks.
- Added a minimal harmonic BAOAB step that uses the existing normal-mode free
  propagator and PILE-L thermostat module.
- For the harmonic validation path, used the standard finite-P canonical
  representation: standard `omega_P`, physical harmonic force scale of `1`, and
  effective bead temperature `P*T`.
- Documented that force scaling, spring convention, and thermostat effective
  temperature are coupled and must not be mixed silently in future production
  drivers.

Validation:

- `test_03_harmonic_quantum_variance` checks finite-P convergence for
  `P = 1, 2, 4, 8, 16, 32`.
- The same test verifies that the legacy `sqrt(P)` prefactor is farther from
  the analytic quantum variance than the standard convention at `P=16`.
- A stochastic `P=16` BAOAB+PILE-L harmonic run sampled the finite-P variance
  within about `2.1%`.

## 2026-04-30: Harmonic Driver And PIMD-Native Outputs

Problem:

- The project needed its first executable assembly path, but it should not jump
  to eReaxFF before harmonic validation is exercised end-to-end.
- The first driver build failed because a long Fortran string exceeded the
  compiler line-length limit.
- The first driver smoke run exposed an output-path bug: `mkdir -p` was called
  without quoting, and this repository path contains a space.

Resolution:

- Added `harmonic_pimd_main.f90` as a thin executable that uses the existing
  ring-polymer, BAOAB, PILE-L, and output modules.
- Added PIMD-native outputs: `pimd_energy.out`, `bead_spread.out`,
  `thermostat_heat.out`, and `centroid.xyz`.
- Split the long header string and quoted the output directory in the
  `execute_command_line` call.
- Added a shell smoke test that runs the driver in a path with spaces and
  verifies expected output files and frame counts.

Validation:

- `test_10_harmonic_driver_outputs` builds and runs `harmonic_pimd_main`,
  verifies output files are nonempty, checks the standard ring-polymer
  convention diagnostic, checks 21 energy frames, and verifies thermostat heat
  changes.

## 2026-04-30: Harmonic Restart

Problem:

- Restart must reproduce a stochastic trajectory when a run is split, so
  positions and velocities alone are insufficient.
- Langevin/PILE restart also needs thermostat heat and RNG state to continue
  diagnostics and random streams consistently.

Resolution:

- Added `restart.f90` with a small harmonic normal-mode restart format.
- Extended `harmonic_pimd_main` with positional arguments for total steps,
  restart input, and restart output.
- Stored normal-mode positions, normal-mode velocities, thermostat heat, bead
  count, step, and RNG state.

Validation:

- The harmonic driver smoke test now compares a continuous 2000-step run against
  a split 1000+1000-step run.  The final restart files match for step, bead
  count, RNG state, thermostat heat, `q_nm`, and `v_nm`.

## 2026-04-30: Structural Stress And Edge-Case Gate

Problem:

- Before coupling eReaxFF, we needed a broader vulnerability check for lifecycle
  and scaling problems that the small tests might miss.
- The first stress target produced a compiler warning around a reused
  `ring_polymer_state`, even though runtime behavior was correct.

Resolution:

- Added `test_11_structural_stress` for:
  - normal-mode orthogonality/roundtrip from `P=1` through `P=128`;
  - repeated ring-polymer initialization;
  - 1024-particle/32-bead harmonic force-manager scaling;
  - 2048-particle/32-bead charge-memory diagnostics;
  - 10,000-particle mixed group assignment.
- Added `test_11_invalid_group_overlap`, a negative test that confirms
  overlapping selectors fail instead of silently assigning a particle twice.
- Reworked the repeated-initialization stress test to avoid the compiler warning
  while still exercising reinitialization.

Validation:

- `make test` now includes the structural stress suite and the overlapping-group
  negative test.  Both pass cleanly.

## 2026-04-30: External V2 Oracle Smoke Gate

Problem:

- The next risk is eReaxFF coupling, but linking legacy V2 directly into the
  new engine before a smoke gate would make failures hard to isolate.
- Building V2 from this repo would write artifacts into the legacy tree, which
  conflicts with the rule that V1/V2 should not be modified in place.
- The first parser assumed key/value output, but the V2 single-point oracle
  writes one header line followed by one compact data line.

Resolution:

- Added `test_08_reaxff_oracle_smoke` as an external validation gate using the
  already-built V2 `phase4_single_point_oracle` and
  `single_point_repeatability` executables.
- Copied V1/V2 fixtures into `build/test_v2_oracle_smoke` and kept all
  generated geometry, logs, and summaries inside the new repo.
- Ran both V2 executables under the existing timeout wrapper and parsed the
  real table format from `phase4_single_point_oracle.out`.
- Kept `src/force/ereaxff_oracle_adapter.f90` as a clean failure stub until the
  wrapped force-provider path is implemented deliberately.

Validation:

- `test_08_reaxff_oracle_smoke` checks a 61-particle static fixture for status
  zero, finite energy, finite max force, finite charge sum, and three repeated
  evaluations with zero reported force and charge differences.
- `make test` now includes this oracle smoke gate and passes all implemented
  regression tests.

## 2026-05-02: Adapter-Mediated V2 Oracle Boundary

Problem:

- The adapter needed to keep the legacy oracle force-only while also producing
  the `P=1` equivalence gate and NaN rejection checks.
- The first attempt used inline shell env assignments, which failed because the
  current command execution path did not interpret them reliably.
- A direct geometry writer produced a malformed legacy atom record, which the
  oracle rejected as an unknown atom type.
- The smoke run also needed the repeatability executable, because the force and
  charge trace is produced there rather than by the single-point oracle.

Resolution:

- Switched the adapter to write the phase4 geometry through the validated legacy
  `make_geo.py` script with explicit electron coordinates.
- Used a shell wrapper in the adapter work directory to run the legacy oracle
  and repeatability binaries with the correct cwd and logged outputs.
- Parsed `phase4_single_point_oracle.out` for summary values and
  `single_point_forces.out` for per-atom forces, charges, and solver memory.
- Added the adapter regression binary to check `P=1` energy, forces, charges,
  repeatability across two calls, and a deliberate NaN rejection path.

Validation:

- `make test_08_reaxff_oracle_smoke` now passes directly.
- `make test` passes the full suite with the adapter gate included.

## 2026-05-02: Force-Manager P=1 Integration

Problem:

- The adapter was validated in isolation, but the next risk is the modular
  manager layer accidentally scaling or reordering oracle results.
- The first manager integration needed to preserve bead-local charge memory and
  force scaling while still calling the oracle through the new manager path.
- A direct manager call without a dedicated summary type would have hidden bead
  status, solver-memory, and charge outputs.

Resolution:

- Added `oracle_bead_summary` and `evaluate_oracle_beads` to
  `bead_force_manager.f90`.
- Restricted the first manager path to explicit `bead_local_acks2` so the
  manager cannot silently cross into centroid or frozen charge modes yet.
- Added a `P=1` regression that checks the manager path against the validated
  V2 oracle values for energy, max force, charge sum, solver memory, and charge
  memory update.

Validation:

- `test_12_oracle_force_manager_p1` passes and confirms that the manager path
  preserves the validated oracle outputs at `P=1`.
- `make test` still passes the entire suite after the manager integration.
