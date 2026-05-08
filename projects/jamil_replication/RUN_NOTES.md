# Jamil Replication Run Notes

## 2026-05-07 Startup

- Generated the first Fig. 3-like target geometries:
  - `geometries/jamil_graphene_112c_112el.geo`
  - `geometries/jamil_graphene_112c_113el.geo`
- The 112C/113El target has a 17.22871 x 17.05200 x 20.00000 A orthorhombic cell.
- Initial strict 0-step run failed because `geo` had the Jamil cell but the copied legacy `control` file still had the old TrackB2 `axis1/axis2/axis3` values. The oracle adapter now rewrites those control axes from configured `CRYSTX` for each call.
- After the cell fix, the 225-particle single-point oracle passed:
  - `energy_kcal = -1.220699157847738E+04`
  - `max_force = 1.002694565160808E-01`
  - `charge_sum = -1.416644579421700E-13`
  - `status = 0`
- The strict run still failed because `single_point_repeatability 3` hit the current 90 s wrapper timeout. For 112C/113El, production/stress configs disable repeatability and rely on the finite single-point oracle path.
- Production/stress configs use `build/oracle_single_point_force_export` as the single-point binary so forces are exported without relying on the repeatability executable.

## Timing Ledger

| Case | Config | Steps | P | Repeatability | Result | Wall time |
| --- | --- | ---: | ---: | --- | --- | --- |
| Strict single-point | `config_jamil_112c_113el_p1_300K.yaml` | 0 | 1 | on | single point passed; repeatability timed out | 143.76 s |
| Production smoke | `config_jamil_112c_113el_p1_300K_prod.yaml` | 1 | 1 | off | passed; 5 oracle calls | 258.13 s |

The production smoke gives about 51.6 s per serial force-export oracle call on this 112C/113El target. A P=16 one-step run would require many serial bead force calls on the current adapter path, so the first P=16 launch should be a short bounded smoke or should wait for parallel force-call dispatch.

## 2026-05-07 P=1 Production Smoke

- Command:
  `./build/ereaxff_mixed_baoab_main projects/jamil_replication/config_jamil_112c_113el_p1_300K_prod.yaml build/jamil_replication/run_jamil_112c_113el_p1_prod_steps1/out`
- Final carbon temperature: 292.972 K.
- Final electron temperature: 130.170 K.
- Final `energy_eV`: -497.510355976870.
- Final `max_force_eV_per_A`: 3.58940646785550E-02.
- Thermostat heat after one step:
  - classical Berendsen: -1.01771946427576E-04 eV
  - electron PILE: -2.13726349951010E-03 eV
- Host assignment check for the input geometry: 113 electrons assigned across 112 unique carbon hosts, as expected for one excess electron.
