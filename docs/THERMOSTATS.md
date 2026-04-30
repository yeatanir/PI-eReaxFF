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

