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

## Phase 5: Harmonic Quantum Benchmarks

For harmonic oscillator tests, compare against:

```text
<x^2> = hbar/(2*m*omega) * coth(beta*hbar*omega/2)
```

Run `P = 1, 2, 4, 8, 16, 32` and report sampled centroid variance, bead
variance, total variance, analytic quantum variance, and percent error.

## Phase 6: Force And Charge Coupling

- harmonic and dummy providers return known forces;
- bead force manager scales physical PIMD forces by `1/P`;
- bead-local charge memory is independent per bead;
- centroid charge contraction is explicit in diagnostics.

## Phase 7: eReaxFF Oracle

- one static geometry smoke test;
- A/B/A/B repeatability test for energies, forces, and charges;
- NaN/Inf rejection;
- timeout wrapper for loader or solver hangs;
- charge sum/status diagnostics per bead.

## Phase 8: Mixed Classical/PIMD Configuration

- one PIMD electron plus fixed classical graphene;
- one PIMD electron plus mobile classical graphene;
- all-classical eReaxFF smoke;
- harmonic-only PIMD;
- selection changes require only config edits.

## Phase 9: Outputs And Restart

- PIMD-native files are written at the configured stride;
- ReaxFF compatibility outputs exist as centroid projections;
- restart can continue a run;
- diagnostics include spring energy, thermostat heat, mode temperatures, bead
  spread, charge by bead, charge centroid, and localization metrics.
