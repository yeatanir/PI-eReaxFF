# Thermostats

Planned thermostat names:

- `none`
- `langevin`
- `pile_l`
- `pile_g`
- `svr`

`pile_l` means local Langevin on centroid and internal normal modes.  `pile_g`
means internal normal-mode Langevin plus global stochastic velocity rescaling on
the centroid.  Until implemented, `pile_g` must fail clearly and must not be
silently treated as `pile_l`.

Every thermostat step must accumulate:

```text
delta_E_bath = KE_before_O - KE_after_O
```

## Current Implementation

Implemented modules:

- `thermostat_base.f90`: kinetic-energy and Maxwell velocity-scale helpers.
- `thermostat_none.f90`: no-op thermostat with zero heat exchange.
- `thermostat_langevin.f90`: exact OU velocity update in Cartesian coordinates.
- `thermostat_pile_l.f90`: local Langevin thermostat in normal-mode space.
- `thermostat_pile_g.f90`: explicit failure stub.
- `thermostat_registry.f90`: supported-name and stub-status checks.

For PILE-L, the centroid uses the configured `centroid_gamma_fs_inv`; internal
modes use:

```text
gamma_k = 2 * lambda * omega_k
```

The current tests validate OU variance, heat accounting, internal mode
temperatures, and the PILE-G failure path.
