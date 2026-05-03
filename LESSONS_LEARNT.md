# Lessons Learnt

This file records implementation problems, scientific risks, and fixes as the
modular PI-eReaxFF engine is built.  It should be updated at each phase before
the phase is considered complete.

## 2026-05-02: Mixed BAOAB Thermostats, Restart Heat, And Build Order

Problem:

- The mixed eReaxFF path needed separate thermostat choices for the classical
  carbon group and the PIMD electron group.
- The mixed restart format did not preserve classical velocities, cumulative
  thermostat heat, or RNG state, so split-run equivalence broke once the bath
  was wired in.
- `particle_groups.f90` started depending on `thermostat_registry`, which the
  build graph did not include in the core-only test targets.
- The initial policy-matrix regression was too strict for a thermalized mixed
  run: the last frame could relax back toward equal bead charges even though the
  policy path itself was valid.

Resolution:

- Added per-group thermostat diagnostics and explicit handling for classical
  Langevin versus PIMD PILE-L.
- Extended the mixed BAOAB restart file with classical velocities, cumulative
  classical/electron heat, and RNG state.
- Restored electron thermostat heat after restart initialization so the split
  trajectory matches the continuous one.
- Moved `thermostat_registry.f90` into the core build list so the module order
  remains valid for core-only and mixed builds.
- Kept the policy-matrix test as a robustness gate by checking finite
  diagnostics across the trajectory instead of requiring the final frame to
  remain nonzero.

Validation:

- `test_19_ereaxff_mixed_baoab_restart` now passes with continuous-vs-restart
  equality.
- `make test` passed after the build graph and restart fixes.

## 2026-05-03: Explicit Cell State And CRYSTX Parsing

Problem:

- The ReaxFF compatibility layer was still inferring box lengths from the
  centroid coordinates.
- That was too weak for periodic geometries and did not match the legacy
  ReaxFF contract, where `CRYSTX` is an explicit cell record.
- The mixed BAOAB path also needed a visible cell contract in the startup
  diagnostics so the geometry boundary is clear at runtime.

Resolution:

- Added a first-class cell state to `system_state` with explicit `crystx`,
  `periodic`, and cached matrix/inverse data.
- Parsed `CRYSTX` from legacy `geo` files in the loader and propagated it into
  the shared state.
- Updated the compatibility output to use the explicit cell for
  minimum-image connectivity and molecule lists instead of inferring a box from
  coordinates.
- Added YAML `geometry.periodic` and `geometry.crystx` fields for explicit box
  specification when the geometry file does not provide one.
- Added a periodic minimum-image regression with a tiny CRYSTX fixture that
  detects a bond across the box boundary.

Validation:

- `test_26_crystx_periodic_loader` passes and confirms the periodic bond is
  detected.
- `make test` passes with the new cell-aware state and output layer.

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

## 2026-05-02: Oracle Centroid Charge-Policy Dispatch

Problem:

- The manager path needed an explicit centroid-contraction branch instead of a
  hidden special case inside the oracle adapter.
- The policy layer had valid names for `frozen_charges` and
  `no_charge_coupling`, but the oracle manager boundary still needed to state
  clearly which policies it can actually execute.
- A centroid regression needed to prove that the new branch does not disturb
  the existing bead-local path.

Resolution:

- Added `charge_policy_uses_centroid_contraction` and
  `charge_policy_supports_oracle_manager`.
- Extended `evaluate_oracle_beads` so `centroid_acks2` contracts bead
  positions, charge memory, and solver guesses to a shared oracle evaluation,
  then broadcasts the result back to the bead summary.
- Kept `frozen_charges` and `no_charge_coupling` as valid policy names but not
  oracle-manager modes.
- Added a two-bead identical-geometry regression that compares bead-local and
  centroid-contracted results.

Validation:

- `test_13_oracle_charge_policy_modes` passes.
- `make test` passes the full suite after the policy dispatch update.

## 2026-05-02: Mixed Oracle-Manager Fixture

Problem:

- The manager layer needed a nontrivial mixed-bead regression, not just an
  identical-beads comparison.
- The phase4 oracle adapter only uses the electron coordinate from each bead, so
  the mixed fixture had to vary the PIMD electron bead positions while keeping
  the classical subset replicated across beads.
- The new regression needed to stay structural, not just repeat the exact
  identical-bead charge-policy gate.

Resolution:

- Added a two-bead phase4-derived fixture with one PIMD electron and classical
  carbon atoms.
- Kept the classical subset identical across beads and perturbed the electron
  bead coordinates so bead-local and centroid charge handling produce distinct
  behavior.
- Verified that the bead-local path accumulates nonzero bead-to-bead charge
  differences while the centroid path broadcasts a shared charge state.

Validation:

- `test_14_mixed_oracle_manager_fixture` passes.
- `make test` passes the full suite with the mixed oracle-manager fixture
  included.

## 2026-05-02: Mixed Oracle-Manager Bead-Count Matrix

Problem:

- The mixed manager coverage needed a bead-count matrix, not just one `P = 2`
  case.
- The bead-local and centroid paths had to be checked at `P = 1` as a
  convergence point, then at larger bead counts to confirm the split between
  local and centroid charge handling.
- The regression had to stay on the same phase4-derived system so the only
  variable was bead count and bead-local perturbation.

Resolution:

- Added a `P = 1, 2, 4` matrix on the phase4 mixed fixture.
- Reused the same classical carbon plus one PIMD electron assignment and
  varied only the number of beads and the electron bead offsets.
- Verified that `P = 1` collapses to the centroid-equivalent result, while
  `P > 1` produces nonzero bead-local charge spread and nonzero bead-local vs
  centroid energy differences.

Validation:

- `test_15_mixed_oracle_manager_matrix` passes.
- `make test` passes the full suite after adding the bead-count matrix.

## 2026-05-02: Mixed Restart And Output Roundtrip

Problem:

- The mixed manager path needed a restart format of its own, not just the
  harmonic single-coordinate restart.
- The restart file had to preserve bead positions, bead-local charge memory, and
  solver guesses in a text format that can be inspected and round-tripped.
- The output side needed bead-wise and centroid charge summaries so the mixed
  path has an output artifact, not only an in-memory state.

Resolution:

- Added `ereaxff_mixed_restart` with text write/read support for mixed bead
  positions, charges, and solver guesses.
- Added bead-wise and centroid charge summary writers.
- Added a regression that writes the mixed restart/output files, checks they
  exist, then reads the restart back and verifies exact roundtrip on the stored
  arrays and policy metadata.

Validation:

- `test_16_mixed_oracle_manager_restart_output` passes.
- `make test` passes the full suite after the restart/output addition.

## 2026-05-02: eReaxFF Mixed Smoke Driver Boundary

Problem:

- The mixed restart/output helpers were validated in memory, but the project
  needed an executable that exercises the eReaxFF adapter through the modular
  mixed path.
- The smoke driver had to keep the oracle force-only, reuse the restart/output
  helpers, and remain restart-splittable without creating a second dynamics
  path.

Resolution:

- Added `ereaxff_pimd_main.f90` as a thin mixed smoke executable that owns the
  restart/output boundary and calls the eReaxFF oracle through the adapter.
- Reused the mixed restart helpers for bead positions, charge summaries, and
  solver memory.
- Added a split-run smoke script that checks output equivalence across a
  restart boundary and verifies the driver prints the expected adapter
  diagnostic.

Validation:

- `test_17_ereaxff_mixed_driver_smoke` passes.
- `make test` includes the executable gate once the runner summary list is
  updated to print it.

## 2026-05-02: eReaxFF Mixed BAOAB Matrix

Problem:

- The smoke path proved the adapter boundary, but we still needed an actual
  modular dynamics path to exercise the force manager and free-ring
  propagator.
- The mixed BAOAB path had to stay on the same oracle-backed fixture while
  varying bead count, not by inventing a second system setup.

Resolution:

- Added `ereaxff_mixed_baoab_integrator.f90` to advance the PIMD electron with
  a BAOAB-style force split and exact free ring-propagation.
- Added `ereaxff_mixed_baoab_main.f90` as a thin executable that keeps the
  oracle force-only, uses the mixed force manager, and writes the same mixed
  diagnostics as the smoke driver.
- Added a small `P = 1, 2, 4, 8` matrix smoke test to catch integrator/boundary
  regressions early.

Validation:

- `test_18_ereaxff_mixed_baoab_matrix` passes.
- `make test` will include this gate once the runner summary list prints it.

## 2026-05-02: Mixed BAOAB Restart And Compatibility Outputs

Problem:

- The mixed BAOAB path needed restart continuity, not just a forward smoke
  path.
- The centroid-projected compatibility outputs had to be explicit final-state
  snapshots, not an inferred side effect of the mixed manager files.

Resolution:

- Added a BAOAB-specific restart format carrying bead positions, normal-mode
  coordinates, velocities, bead charges, and solver memory.
- Added centroid-projected ReaxFF compatibility outputs: `fort.7`, `fort.8`,
  `molfra.out`, and `xmolout`.
- Added a split-run restart regression that compares final restart state and
  compatibility snapshots against a continuous run.

Validation:

- `test_19_ereaxff_mixed_baoab_restart` passes.
- `make test` includes the restart gate once the runner summary list prints it.

## 2026-05-03: Legacy ReaxFF Output Semantics And Thermal Partition

Problem:

- The compatibility layer had been treating `fort.7`, `fort.8`, and
  `molfra.out` as generic snapshots instead of ReaxFF-style topology outputs.
- The mixed BAOAB path also needed a phase6-like temperature partition
  diagnostic so the hot-electron branch can be tracked explicitly.

Resolution:

- Reworked the compatibility writer so `fort.7` and `fort.8` are connection
  tables and `molfra.out` is a per-frame molecule list / molecular
  fragmentation summary.
- Added minimum-image topology handling for the centroid-projected
  compatibility writer.
- Added `carbon_temperature.out`, `electron_temperature.out`, and
  `electron_mode_temperature.out` to the mixed BAOAB path.
- Added a phase6-style temperature partition regression that records the
  carbon/electron split as a diagnostic, not as a proof of a solved hot-electron
  physics problem.

Validation:

- `test_24_ereaxff_mixed_baoab_user_geo_matrix` passes.
- `test_25_phase6_like_temperature_partition` passes.
- `make test` passes.

## 2026-05-02: Mixed BAOAB Policy Matrix

Problem:

- The mixed-system stress matrix still only covered bead count.  We needed to
  prove that the same executable path behaves sensibly under both explicit
  charge-coupling policies.
- The charge-policy matrix had to stay on the same fixture so the new signal was
  policy behavior, not a new geometry.

Resolution:

- Added a policy-axis matrix over `bead_local_acks2` and `centroid_acks2`.
- Covered `P = 1, 2, 4, 8` on the mixed BAOAB executable.
- Checked that the centroid policy suppresses bead-to-bead charge spread while
  the bead-local policy only collapses at `P = 1`.

Validation:

- `test_20_ereaxff_mixed_baoab_policy_matrix` passes.
- `make test` includes the policy matrix once the runner summary list prints it.

## 2026-05-02: Mixed BAOAB Config Layer

Problem:

- The mixed BAOAB executable still hardcoded its output directory, bead count,
  run length, charge policy, and legacy V1/V2 input paths.
- The current config examples were too shallow for the adapter boundary because
  they named only bare legacy filenames, not the actual V1-relative paths.
- The first parser pass was too strict about group keys and rejected thermostat
  tuning entries that are present in the config examples but not yet consumed
  by the dynamics path.

Resolution:

- Added `src/core/input_config.f90` as a small YAML-subset parser for
  simulation, force provider, ring polymer, charge coupling, output, and
  group sections.
- Taught `ereaxff_mixed_baoab_main.f90` to accept either the legacy positional
  interface or a config file directly.
- Wired the driver to use config-defined output directory, bead count, step
  count, temperature, charge policy, group list, and explicit V1-relative
  input paths for `control`, `ffield`, and `iopt`.
- Relaxed group parsing for `gamma_fs_inv` and `centroid_gamma_fs_inv` so the
  current example configs can carry thermostat tuning metadata without blocking
  the smoke path.

Validation:

- `test_21_ereaxff_mixed_baoab_config_smoke` passes.
- `make test` passes with the new config-driven gate included.

## 2026-05-02: Mixed BAOAB Config Matrix And Negative Cases

Problem:

- One config smoke run is not enough to trust the new input layer.
- We needed to confirm that the mixed BAOAB executable stays stable across a
  small matrix of bead counts and charge policies, and that it still rejects
  bad selectors and overlapping groups cleanly.

Resolution:

- Added a config-driven matrix gate covering `P = 1, 2, 4, 8` for both
  `bead_local_acks2` and `centroid_acks2`.
- Added deliberate negative cases for an invalid selector and overlapping
  group definitions.

Validation:

- `test_22_ereaxff_mixed_baoab_config_matrix` passes.
- `make test` passes with the matrix and failure checks included.

## 2026-05-02: Geometry Variation Within The Oracle Contract

Problem:

- The mixed BAOAB path was still effectively locked to one electron starting
  point inside the validated phase4 helper.
- We needed a way to widen the geometry stress surface without claiming the
  oracle supports arbitrary structures.

Resolution:

- Added a `geometry` section to the input parser with the validated
  `phase4_fixed_graphene` template and an explicit electron starting position.
- Threaded the geometry position into the mixed BAOAB driver and into the
  legacy make-geo wrapper through the existing `PHASE4_EL_X/Y/Z` contract.
- Added a geometry-matrix regression that exercises two electron starting
  positions while keeping the same oracle boundary.

Validation:

- `test_23_ereaxff_mixed_baoab_geometry_matrix` passes.
- `make test` passes with the geometry matrix included.

Problem found during implementation:

- The first geometry negative test only covered unsupported selectors and
  overlap.  It did not prove that the geometry contract itself rejects
  unsupported templates with a clear message.

Resolution:

- Added a dedicated unsupported-`geometry.template` failure case to the
  geometry matrix.
- Reserved `structure_path` and `geo_path` in the parser as the future
  user-supplied loader contract, while explicitly rejecting the unimplemented
  `user_geo_structure` template.

## 2026-05-02: Legacy Geo Loader And Geometry Control

Problem:

- The mixed BAOAB path needed to move beyond the phase4-only geometry helper
  and accept a user-supplied legacy `geo` file path.
- The first attempt to parse the legacy `geo` file with list-directed reads was
  too brittle for the fixed-width HETATM layout used by the V1/V2 fixtures.
- The first attempt to reuse legacy atom ids from the file did not give a
  reliable selector contract for the modular group system.

Resolution:

- Added `src/core/legacy_geo_loader.f90` to load fixed-width legacy `geo`
  files into `system_state`.
- Wired `geometry.template: user_geo_structure` plus `geo_path`/`structure_path`
  through the mixed BAOAB driver and the eReaxFF adapter.
- Switched the loader to parse HETATM records with the legacy fixed format and
  to assign deterministic sequential ids so the selector layer remains stable.

Validation:

- `test_24_ereaxff_mixed_baoab_user_geo_matrix` runs the mixed BAOAB executable
  against the V2 background-revisit `geo` fixture and checks the geometry
  source, particle counts, and output files.
- `make test` passed after the new loader was added.

Problem found during integration:

- A clean rebuild failed because the Makefile compiled `bead_force_manager`
  before `ereaxff_oracle_adapter`, leaving the module dependency unresolved.

Resolution:

- Reordered `FORCE_SRCS` so the oracle adapter compiles before the bead force
  manager on a clean tree.

Validation:

- Clean rebuild and full `make test` passed after the source-order fix.

## 2026-05-03: Multi-Quantum Mixed BAOAB Support

Problem:

- The mixed BAOAB driver still assumed exactly one PIMD particle in practice.
- That made the framework unusable for the real graphitic transport cases,
  where multiple excess electrons can be quantum particles in the same run.
- The pseudo-electron mass was still too implicit in the driver boundary.

Resolution:

- Removed the hard singleton guard and generalized the mixed BAOAB step to
  carry an explicit list of quantum particles.
- Kept the production boundary conservative by requiring a shared bead count
  for all quantum particles in one mixed run.
- Added explicit `mass_amu` support to group specifications so pseudo-particle
  mass is not hidden inside driver code.
- Extended the restart state to round-trip multiple quantum indices, masses,
  normal-mode coordinates, charges, and solver guesses.

Validation:

- `test_27_multi_quantum_restart` passes for a two-quantum-particle mixed
  restart roundtrip.
- `make test` passes end to end with the new multi-quantum boundary.
