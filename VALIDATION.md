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

Current status: implemented an external V2 oracle smoke gate that runs the
existing `phase4_single_point_oracle` and `single_point_repeatability`
executables from V2 in a repo-local work directory.  The test copies legacy
fixtures into `build/`, runs both commands under timeout, and verifies finite
single-point diagnostics plus repeatability of force and charge results.  The
new modular `ereaxff_oracle_adapter.f90` is still a failure stub; the next step
is linking this validated oracle boundary through the new force-provider
interface without moving dynamics into the oracle.

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

## Structural Stress Gate

Before eReaxFF coupling, the suite includes stress and negative tests for:

- normal-mode orthogonality/roundtrip across odd, even, small, and larger bead
  counts through `P=128`;
- repeated ring-polymer reinitialization;
- large bead force-manager arrays;
- large bead-local charge memory;
- 10,000-particle group assignment;
- overlapping group rejection.
