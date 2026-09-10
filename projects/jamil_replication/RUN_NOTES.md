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
| Parallel PIMD smoke | `config_jamil_112c_113el_p16_300K_prod.yaml` | 1 | 16 | off | passed with `OMP_NUM_THREADS=16`; 80 oracle calls | 4066.99 s |

The production smoke gives about 51.6 s per serial force-export oracle call on this 112C/113El target.
The adapter now supports three dispatch modes:

- `serial`
- `openmp_all`
- `bounded_pool`

For the legacy external oracle backend, the recommended baseline is:

```bash
export OMP_NUM_THREADS=1
```

with config:

```yaml
oracle_dispatch_mode: bounded_pool
oracle_max_workers: 2   # or 4
oracle_timing_log: true
oracle_stagger_launch_ms: 0
oracle_fail_fast: true
```

The P=16 smoke confirmed that 16 concurrent legacy force-export binaries scale poorly on this machine. The next benchmark sweep should start from bounded batches, not full fan-out.

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

## 2026-05-07 P=16 Parallel Smoke

- Command:
  `env OMP_NUM_THREADS=16 ./build/ereaxff_mixed_baoab_main projects/jamil_replication/config_jamil_112c_113el_p16_300K_prod.yaml build/jamil_replication/run_jamil_112c_113el_p16_prod_steps1_parallel/out`
- Result: passed one 0.1 fs step.
- Oracle calls: 80 complete `adapter_call_*` directories.
- Adapter work size: 715 MB.
- Final carbon temperature: 292.973 K.
- Final electron temperature: 505.433 K.
- Final `energy_eV`: -511.163605892981.
- Final `max_force_eV_per_A`: 2.24834187306845E-03.
- Thermostat heat after one step:
  - classical Berendsen: -1.01771946427576E-04 eV
  - electron PILE: 6.48222086563056E-01 eV

## 2026-05-08 Dispatch Sweep Template

Use the production P=16 geometry with `OMP_NUM_THREADS=1` and sweep `oracle_max_workers` over `1, 2, 4, 6, 8`.

Command template:

```bash
OMP_NUM_THREADS=1 ./build/ereaxff_mixed_baoab_main \
  projects/jamil_replication/config_jamil_112c_113el_p16_300K_prod.yaml \
  build/jamil_replication/run_jamil_112c_113el_p16_w4/out
```

Benchmark table:

| Workers | Dispatch mode | Wall time | Oracle calls | Batch count | Min call | Mean call | Max call | Notes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | bounded_pool |  |  |  |  |  |  | serial baseline |
| 2 | bounded_pool |  |  |  |  |  |  |  |
| 4 | bounded_pool |  |  |  |  |  |  |  |
| 6 | bounded_pool |  |  |  |  |  |  |  |
| 8 | bounded_pool |  |  |  |  |  |  |  |

Interpretation:
- pick the smallest wall time that keeps all forces finite and does not trigger long-tail stalls.
- if worker counts above 4 slow down, keep the backend bounded and do not fan out all 16 beads.

## Oracle Call Count Audit

The observed `80` oracle calls in the P=16 one-step smoke are not extra validation traffic. They come from the current mixed BAOAB schedule:

- `initial_force`
- `group_pre_force`
- `quantum_pre_force`
- `quantum_post_force`
- `group_post_force`

That is `5` oracle evaluations per step, multiplied by `16` bead-local jobs per evaluation.

Production-only validation traffic is already controlled by `repeatability_check`; the current production configs set it to `false`. No additional validation-only call path was identified in the P=16 smoke. So the next gain should come from bounded dispatch tuning, not from hiding extra oracle calls.

The timing log now records `stage_kind` and `stage_label` for each call, so the 80-call trace is visible in `oracle_timing.log` as a stage-by-stage series instead of a single undifferentiated bead-local batch.

## 2026-09-09 Persistent Worker Pool

Added an exact `persistent_worker_pool` backend. Each pool slot is an independent persistent worker process with its own initialized legacy ReaxFF COMMON state, so bead-local PIMD force calls can run concurrently without sharing thread-unsafe ReaxFF state.

Use:

```yaml
oracle_backend_mode: persistent_worker_pool
persistent_worker_binary: "/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_pimd_modular/build/persistent_oracle_worker"
oracle_dispatch_mode: bounded_pool
oracle_max_workers: 4
oracle_timing_log: true
repeatability_check: false
```

P=4/one-step verification:

- Command:
  `env OMP_NUM_THREADS=4 ./build/ereaxff_mixed_baoab_main projects/jamil_replication/config_jamil_112c_113el_p4_pool_verify.yaml build/jamil_replication/verify_p4_pool_steps1/out`
- Result: passed one 0.1 fs step with 20 oracle requests and four workers shut down cleanly.
- Parallel worker initialization took about `34.3 s` before the first batch.
- Exact four-bead batches took `19.5-20.5 s` each on the 225-particle Jamil system.
- Final carbon temperature: `292.973 K`.
- Final electron physical temperature: `316.700 K`; extended ring-polymer kinetic temperature: `1266.801 K`, matching the expected P=4 diagnostic scaling.
- Thermostat heat after one step:
  - classical Berendsen: `-1.01771946427576E-04 eV`
  - electron PILE: `1.43071183013831 eV`

This is the first real throughput win: the previous direct single-worker P=4 timing was about four serial bead calls per force stage, while the pool now makes each P=4 force stage cost close to one 225-particle direct call.

## 2026-05-08 Partial Dispatch Sweep

I started the strict bounded-pool sweep with `JAMIL_DISPATCH_SKIP_W1=1` and `OMP_NUM_THREADS=1`.

Observed so far for `oracle_max_workers=2`:

- step 0 `initial_force` completed as 5 bead batches of 2 jobs each
- each 2-job batch was about `96-97 s`
- each bead call was still about `48-49 s`
- after 5 batches the run had already spent about `484 s` just on `initial_force`

That is the same negative shape as before: two-way concurrency did not improve throughput. The full sweep was stopped before `workers=4,6,8` because the `workers=2` result already showed that bounded process-level concurrency is not fixing the legacy backend cost.
