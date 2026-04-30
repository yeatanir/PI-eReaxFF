# Design

## Architecture

The engine is organized around a small set of ownership boundaries:

- `core`: numeric kinds, units, constants, input/state containers, selectors,
  restart, errors, and RNG.
- `force`: force provider interfaces, harmonic/dummy providers, the eReaxFF
  oracle adapter, bead force orchestration, and charge coupling policy.
- `pimd`: ring-polymer normal modes, free propagation, spring energy, and PIMD
  estimators.
- `thermostat`: first-class thermostat modules selected by name and bound to
  particle groups or normal modes.
- `integrator`: BAOAB/OBABO-style split integrators that call thermostat and
  force interfaces without owning their implementation details.
- `output`: PIMD-native diagnostics and ReaxFF compatibility projections.
- `drivers`: thin programs that parse input, build components, run loops, and
  delegate all physics operations.
- `ipi_socket`: future socket interface notes and stubs that reuse the same
  eReaxFF force oracle adapter.

## Force Oracle Boundary

The eReaxFF oracle adapter must only initialize force-field state and evaluate
energy, forces, charges, virial/status, and diagnostics for a supplied geometry.
It must not propagate positions, update velocities, thermostat, own MD state, or
write trajectories except through an output adapter.

## PIMD Convention

The central convention is `standard_ipi`:

```text
omega_P = P * kB * T / hbar
omega_k = 2 * omega_P * sin(pi*k/P), k = 0..P-1
omega_0 = 0
```

The implementation stores frequencies in `fs^-1` using `kB` in `eV/K` and
`hbar` in `eV*fs`.  All modules and tests must call the shared PIMD modules
instead of duplicating transform or frequency logic in drivers.

The legacy `sqrt(P)` prefactor is not supported as a production convention.  If
a future nonstandard Hamiltonian is added, it must have an explicit name,
derivation, unit tests, and diagnostics.

## Normal-Mode Transform

The initial transform is the real orthonormal closed-ring DFT matrix:

- row 1 is the centroid mode, `1/sqrt(P)`;
- paired cosine/sine rows span internal modes;
- for even `P`, the final Nyquist row is `(-1)^j/sqrt(P)`.

Transforms use:

```text
q_nm    = C * q_bead
q_bead  = transpose(C) * q_nm
```

The physical centroid is `q_nm(1) / sqrt(P)`, equal to the bead average.

## Charge Coupling Policy

The planned policy names are:

- `bead_local_acks2`: independent force/charge solve and charge memory per
  bead;
- `centroid_acks2`: centroid charge solve distributed as a documented
  approximation;
- `frozen_charges`: fixed charges for validation/debugging;
- `no_charge_coupling`: harmonic/dummy tests.

The default production-like full PIMD path will be `bead_local_acks2`.

## Outputs

PIMD-native outputs will be separate from ReaxFF compatibility outputs.
Compatibility files such as `fort.7`, `fort.8`, `molfra.out`, and `xmolout`
are centroid-projected classical views unless explicitly configured otherwise.

