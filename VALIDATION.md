# Validation Ladder

Validation proceeds from algorithmic invariants to coupled eReaxFF behavior.

## Phase 1: Constants And Units

- print core constants;
- verify selected unit conversions;
- keep constants in one module.

Current status: implemented for core constants, unit conversion, and
deterministic RNG smoke testing.

## Phase 2: Normal Modes And Ring Polymer

Required gates:

- `max |C*C^T - I| < 1e-12`;
- bead-to-normal-mode-to-bead roundtrip `< 1e-12`;
- centroid mode equals bead average;
- `omega_P = P*kB*T/hbar`;
- `omega_k = 2*omega_P*sin(pi*k/P)`;
- explicit regression fails if `sqrt(P)*kB*T/hbar` is restored.

## Phase 3: Free Ring Polymer

- exact free ring-polymer propagation conserves free ring-polymer energy in
  NVE;
- centroid drifts freely and internal modes rotate exactly.

Current status: implemented with a 20,000-step exact propagation conservation
test.

## Phase 4: Thermostats

- Langevin OU velocity variance approaches `kBT/m`;
- PILE-L internal mode temperatures are reasonable;
- thermostat heat exchange is accumulated;
- invalid names fail cleanly;
- PILE-G fails with a clear message until implemented.

Current status: implemented for `none`, Cartesian Langevin, PILE-L, registry
name checks, heat accounting, and a PILE-G negative test.  PILE-G remains a
stub by design.

## Phase 5: Harmonic Quantum Benchmarks

For harmonic oscillator tests, compare against:

```text
<x^2> = hbar/(2*m*omega) * coth(beta*hbar*omega/2)
```

Run `P = 1, 2, 4, 8, 16, 32` and report sampled centroid variance, bead
variance, total variance, analytic quantum variance, and percent error.

Current status: implemented finite-P bead variance convergence checks for
`P = 1, 2, 4, 8, 16, 32` and a stochastic `P=16` BAOAB+PILE-L harmonic
variance check.  The test also verifies that the legacy `sqrt(P)` prefactor is
farther from the analytic quantum variance than the standard convention.

## Phase 6: Force And Charge Coupling

- harmonic and dummy providers return known forces;
- bead force manager scales physical PIMD forces by `1/P`;
- bead-local charge memory is independent per bead;
- centroid charge contraction is explicit in diagnostics.

Current status: implemented for dummy/harmonic force providers, bead-averaged
potential, centralized `1/P` force scaling, policy-name validation, and
bead-local charge memory diagnostics.  eReaxFF/ACKS2 coupling is not connected
yet.

## Phase 7: eReaxFF Oracle

- one static geometry smoke test;
- A/B/A/B repeatability test for energies, forces, and charges;
- NaN/Inf rejection;
- timeout wrapper for loader or solver hangs;
- charge sum/status diagnostics per bead.

Current status: implemented the adapter-mediated V2 oracle gate.  The new
`ereaxff_oracle_adapter.f90` writes the legacy phase4 geometry, runs
`phase4_single_point_oracle` and `single_point_repeatability` through a
repo-local wrapper, returns eV energy with eV/Angstrom forces and charges, and
rejects NaN input.  The regression checks `P=1` adapter equivalence, adapter
repeatability, and the negative NaN path before any mixed-classical/PIMD work
builds on top of it.

Current status for the manager layer: `evaluate_oracle_beads` now exercises the
same adapter at `P=1`, updates bead-local charge memory, and verifies that the
force manager preserves the validated oracle outputs when the oracle is reached
through the modular bead loop.

Current status for charge-policy dispatch: `bead_local_acks2` and
`centroid_acks2` are both accepted by the oracle manager.  The centroid branch
contracts bead positions, charge memory, and solver guesses to a shared
centroid evaluation, then broadcasts the result back to the bead summary.  The
policy layer still treats `frozen_charges` and `no_charge_coupling` as valid
names, but not as oracle-manager modes.

Current status for mixed oracle-manager fixtures: a nontrivial phase4-derived
test now exercises one PIMD electron plus classical carbon atoms on a mixed
two-bead geometry.  It verifies that the classical subset is replicated across
beads, that bead-local charge memory diverges when the PIMD electron bead
positions differ, and that centroid contraction collapses the bead-local charge
differences back to a shared state.

Current status for mixed oracle-manager bead-count coverage: the same fixture
has now been run at `P = 1, 2, 4`.  That matrix confirms that the manager path
reduces to the centroid-equivalent result at `P = 1`, while the bead-local path
splits away at larger bead counts in a controlled, reproducible way.

Current status for mixed restart/output coverage: the mixed manager fixture now
round-trips through a text restart file containing positions, bead-local
charges, and solver guesses.  The regression also writes bead-wise and
centroid-charge summary files and verifies that the restart reload reproduces
the stored mixed state exactly.

Current status for the executable mixed eReaxFF path: `ereaxff_pimd_main`
wraps the same mixed restart/output helpers and the V2 oracle adapter in a
standalone smoke driver.  The smoke test runs the executable twice with a
restart boundary, verifies the output files are identical across the split
run, and confirms the driver prints the expected eReaxFF adapter diagnostic.

Current status for the modular mixed dynamics path: `ereaxff_mixed_baoab_main`
advances the PIMD electron with the force-manager and free-ring propagator
layer, then checks a small bead-count matrix at `P = 1, 2, 4, 8` with the same
oracle-backed mixed fixture.

Current status for restart continuity: the mixed BAOAB driver now writes a
BAOAB-specific restart containing bead positions, normal-mode coordinates, and
velocities.  The restart test confirms that a split run matches the continuous
run in its final restart state and centroid-projected compatibility outputs.

## Phase 8: Mixed Classical/PIMD Configuration

- one PIMD electron plus fixed classical graphene;
- one PIMD electron plus mobile classical graphene;
- all-classical eReaxFF smoke;
- harmonic-only PIMD;
- selection changes require only config edits.

Current status: implemented the system-state, selector, and particle-group
assignment layer needed for mixed classical/PIMD/frozen selection.  The current
test covers mobile classical carbons, one PIMD electron, and frozen background
electrons.

## Phase 9: Outputs And Restart

- PIMD-native files are written at the configured stride;
- ReaxFF compatibility outputs exist as centroid projections;
- restart can continue a run;
- diagnostics include spring energy, thermostat heat, mode temperatures, bead
  spread, charge by bead, charge centroid, and localization metrics.

Current status: implemented first PIMD-native harmonic outputs for
`pimd_energy.out`, `bead_spread.out`, `thermostat_heat.out`, and `centroid.xyz`.
The harmonic driver writes and reads restart state, including RNG state, and the
test suite verifies continuous-vs-split restart equivalence.  ReaxFF
compatibility outputs remain future work.

The mixed eReaxFF smoke driver now extends that restart/output boundary to the
oracle-backed mixed path.  It writes `pimd_energy.out`,
`charge_by_bead.out`, `charge_centroid.out`, `electron_localization.out`, and
`restart.dat`, then checks split-run equivalence through the same executable
path.

The mixed BAOAB driver extends that boundary one layer further by integrating
the PIMD electron against the eReaxFF force manager and checking a small
bead-count matrix with the modular integrator boundary in place.

The restart gate then compares the continuous and split mixed BAOAB runs at
the end of the trajectory and verifies the final restart and compatibility
snapshots match exactly.

## Structural Stress Gate

Before eReaxFF coupling, the suite includes stress and negative tests for:

- normal-mode orthogonality/roundtrip across odd, even, small, and larger bead
  counts through `P=128`;
- repeated ring-polymer reinitialization;
- large bead force-manager arrays;
- large bead-local charge memory;
- 10,000-particle group assignment;
- overlapping group rejection.
