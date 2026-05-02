# Project Status

Generated: 2026-05-02

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
- [x] Add system state and particle-group modules.
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
- [x] Add eReaxFF oracle adapter and adapter-mediated V2 smoke/repeatability gate.
- [x] Add oracle force-manager P=1 regression gate.
- [x] Add oracle centroid charge-policy regression gate.
- [x] Add mixed oracle-manager fixture regression gate.
- [x] Add mixed oracle-manager bead-count matrix regression gate.
- [x] Add mixed oracle-manager restart/output regression gate.
- [x] Add eReaxFF mixed smoke driver regression gate.
- [x] Add eReaxFF mixed BAOAB dynamics matrix regression gate.
- [x] Add eReaxFF mixed BAOAB restart-continuity regression gate.
- [x] Add eReaxFF mixed BAOAB charge-policy matrix regression gate.
- [x] Add YAML-subset input parser for the mixed BAOAB driver.
- [x] Add config-driven mixed BAOAB smoke regression gate.
- [x] Add config-driven mixed BAOAB matrix and negative-case regression gate.
- [x] Add geometry-variation regression gate for the mixed BAOAB oracle path.
- [x] Add system state module.
- [x] Add selector module for `all`, `element == X`, `id == N`, `id != N`, and `id in [...]`.
- [x] Add particle group assignment module.
- [x] Add mixed classical/PIMD/frozen group regression test.
- [x] Add harmonic finite-P quantum variance utility.
- [x] Add minimal harmonic BAOAB PIMD step.
- [x] Add harmonic quantum variance regression for `P = 1, 2, 4, 8, 16, 32`.
- [x] Add stochastic BAOAB+PILE-L harmonic variance smoke test.
- [x] Add PIMD-native output helpers.
- [x] Add thin `harmonic_pimd_main` driver.
- [x] Add harmonic driver output smoke test.
- [x] Add harmonic restart read/write scaffold.
- [x] Add split-run restart equivalence regression including RNG state.
- [x] Add structural stress/edge-case test suite.
- [x] Add negative overlapping-group failure test.
- [x] Connect eReaxFF oracle adapter to legacy V2 oracle.

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
test_08_reaxff_oracle_smoke: PASS
test_09_mixed_classical_pimd_config: PASS
test_10_harmonic_driver_outputs: PASS
test_11_structural_stress: PASS
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
test_11_invalid_group_overlap: PASS
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

The eReaxFF oracle adapter now owns the legacy V2 subprocess boundary.  It
generates the phase4 geometry through the validated legacy script, runs
`phase4_single_point_oracle` and `single_point_repeatability` in a repo-local
work directory, returns eV energy with eV/Angstrom forces and charges, and
rejects NaN input.  The adapter gate verifies `P=1` reproduction against the
validated V2 oracle, repeatability across two adapter calls, and a negative
NaN-rejection path.

The force-manager gate then exercises the same adapter through
`evaluate_oracle_beads` at `P=1`, confirming that bead-local charge memory,
scaled forces, and PIMD averaging remain consistent when the oracle is reached
through the modular manager layer.

The centroid-policy gate exercises the same manager path with explicit
`centroid_acks2` contraction on a two-bead identical-geometry fixture and
confirms that the manager classifies `bead_local_acks2` and
`centroid_acks2` as oracle-supported while leaving `frozen_charges` and
`no_charge_coupling` outside the oracle manager boundary.

The mixed oracle-manager fixture then exercises a nontrivial group selection on
the phase4 template: one PIMD electron, classical carbon atoms, and a bead-
local versus centroid charge-policy comparison on a mixed bead geometry where
the electron bead positions differ.  It confirms that the classical subset is
replicated across beads, bead-local charge memory diverges across beads, and
centroid contraction broadcasts a shared charge state.

The mixed oracle-manager matrix extends that check across `P = 1, 2, 4`, which
confirms that the bead-local path collapses to the centroid path at `P = 1`
and splits away at larger bead counts while keeping the classical subset
replicated across all beads.

The mixed restart/output gate writes a text restart containing bead positions,
bead-local charges, and solver guesses, then reads it back and verifies exact
roundtrip.  It also writes bead-wise and centroid charge summaries so the mixed
manager path has a tested output/restart boundary, not just an in-memory
fixture.

The eReaxFF mixed smoke driver exercises that same mixed restart/output
boundary through the executable path.  It calls the oracle adapter through the
modular mixed driver, writes the same restart and charge summaries, and
verifies split-run equivalence across a restart boundary.

The eReaxFF mixed BAOAB driver then advances the PIMD electron against the
eReaxFF-backed force manager with the modular integrator layer and verifies a
small bead-count matrix at `P = 1, 2, 4, 8`.

The mixed BAOAB restart gate then confirms that a split run reproduces the
same final restart state and centroid-projected compatibility snapshots as the
continuous run.

The mixed BAOAB policy matrix then exercises both `bead_local_acks2` and
`centroid_acks2` across `P = 1, 2, 4, 8`, which widens the stress matrix to
cover charge-coupling policy as well as bead count.

The harmonic quantum variance test verifies finite-P convergence toward the
analytic quantum oscillator variance, verifies that the legacy `sqrt(P)`
prefactor is much worse at `P=16`, and samples a `P=16` BAOAB+PILE-L trajectory
within about `2.1%` of the finite-P reference.

The harmonic driver smoke test builds `harmonic_pimd_main`, runs it in a path
containing spaces, and verifies PIMD-native outputs:
`pimd_energy.out`, `bead_spread.out`, `thermostat_heat.out`, and `centroid.xyz`.
It also verifies continuous-vs-split restart equivalence.

The structural stress suite covers normal-mode transforms through `P=128`,
repeated ring-polymer reinitialization, a 1024-particle/32-bead force-manager
case, 2048-particle/32-bead charge memory, 10,000-particle group assignment, and
overlapping-group rejection.
