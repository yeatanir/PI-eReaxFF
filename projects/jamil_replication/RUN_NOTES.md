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

## 2026-09-10 P Sweep At 4 Workers

Command:

```bash
projects/jamil_replication/run_jamil_p_sweep.sh
```

Environment defaults:

- `JAMIL_P_SWEEP_WORKERS=4`
- `JAMIL_P_SWEEP_OMP_THREADS=4`
- `JAMIL_P_SWEEP_STEPS=1`
- `JAMIL_P_SWEEP_DT_FS=0.1`

Results:

| P | Wall time | Oracle calls | Batches | Mean batch | Electron T | Extended T | Expected extended T | Ratio | Carbon T | Quantum heat | Bead RMS spread | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 4 | 137.893 s | 20 | 5 | 20.254 s | 316.700 K | 1266.801 K | 1200.000 K | 4.000 | 292.973 K | 1.430712 eV | 0.005341 A | PASS |
| 8 | 237.698 s | 40 | 10 | 19.954 s | 368.316 K | 2946.529 K | 2400.000 K | 8.000 | 292.973 K | -4.582652 eV | 0.008395 A | PASS |
| 16 | 432.735 s | 80 | 20 | 19.833 s | 505.433 K | 8086.934 K | 4800.000 K | 16.000 | 292.973 K | 10.371555 eV | 0.013097 A | PASS |

Interpretation:

- The pool speed scales as expected: the per-batch force cost stays about `20 s`, and wall time grows with `ceil(P / 4)`.
- The extended-to-physical diagnostic ratio is exactly P, so that part is just the scaled-Hamiltonian bookkeeping.
- The physical electron temperature still increases with P at fixed `dt_fs=0.1`, so the hot-electron issue is real and is not only a diagnostic-label problem.
- Carbon temperature remains stable near `293 K`, so the P-sensitive heating is isolated to the electron PIMD path.

## 2026-09-10 P Sweep At 4 Workers, 0.05 fs

Command:

```bash
JAMIL_P_SWEEP_DT_FS=0.05 projects/jamil_replication/run_jamil_p_sweep.sh
```

Results:

| P | Wall time | Oracle calls | Batches | Mean batch | Electron T | Extended T | Expected extended T | Ratio | Carbon T | Quantum heat | Bead RMS spread | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 4 | 133.612 s | 20 | 5 | 19.957 s | 322.719 K | 1290.876 K | 1200.000 K | 4.000 | 292.973 K | 1.068152 eV | 0.002693 A | PASS |
| 8 | 238.493 s | 40 | 10 | 20.073 s | 379.417 K | 3035.338 K | 2400.000 K | 8.000 | 292.973 K | -3.048897 eV | 0.004251 A | PASS |
| 16 | 432.502 s | 80 | 20 | 19.666 s | 489.111 K | 7825.775 K | 4800.000 K | 16.000 | 292.972 K | 5.482156 eV | 0.006594 A | PASS |

Interpretation:

- Halving the timestep reduces bead spread and P=16 quantum heat, but it does not remove the P-dependent physical electron heating.
- P=16 physical electron temperature only drops from `505.433 K` at `0.1 fs` to `489.111 K` at `0.05 fs`.
- The next check should vary electron thermostat coupling/initialization, not just timestep.

## 2026-09-10 P=16 Centroid Thermostat Sweep

Command:

```bash
JAMIL_P_SWEEP_CASES=16 JAMIL_P_SWEEP_GAMMAS="0.001 0.005 0.01" \
  JAMIL_P_SWEEP_DT_FS=0.05 projects/jamil_replication/run_jamil_p_sweep.sh
```

Results:

| Centroid gamma | Wall time | Oracle calls | Batches | Mean batch | Electron T | Extended T | Expected extended T | Ratio | Carbon T | Quantum heat | Bead RMS spread | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.001 fs^-1 | 442.092 s | 80 | 20 | 20.327 s | 489.111 K | 7825.775 K | 4800.000 K | 16.000 | 292.972 K | 5.482156 eV | 0.006594 A | PASS |
| 0.005 fs^-1 | 441.559 s | 80 | 20 | 20.256 s | 489.024 K | 7824.388 K | 4800.000 K | 16.000 | 292.972 K | 5.709646 eV | 0.006594 A | PASS |
| 0.010 fs^-1 | 447.885 s | 80 | 20 | 20.406 s | 488.959 K | 7823.349 K | 4800.000 K | 16.000 | 292.972 K | 5.879248 eV | 0.006594 A | PASS |

Interpretation:

- Increasing centroid friction across this range does not fix the P=16 physical electron heating.
- The physical electron temperature stays near `489 K`, while the carbon temperature remains near `293 K`.
- The heat term increases with gamma even though the measured physical electron temperature barely moves, so this is not a simple underdamped-centroid problem.
- The next diagnostic should split electron kinetic energy by normal mode and run a force-off/harmonic-only P=16 check. That will separate incorrect initialization or thermostat sampling from force-coupled heating.

## 2026-09-10 P=16 Force-Off Mode Diagnostics

Added `electron_mode_partition.out` to split the primary PIMD electron kinetic energy into:

- centroid extended mode temperature
- internal extended mode temperature averaged over modes `2:P`
- all-mode extended temperature
- physical aggregate electron temperature
- centroid/internal kinetic-energy fraction

Also added a config-driven `force_provider.type: dummy` mode for the mixed BAOAB driver. This keeps the Jamil geometry, masses, velocity initialization, free ring propagation, and PILE-L thermostat, but replaces all force evaluations with zero energy and zero force.

One-step force-off check:

```bash
JAMIL_FORCE_OFF_STEPS=1 JAMIL_FORCE_OFF_DT_FS=0.05 JAMIL_FORCE_OFF_P=16 \
  JAMIL_FORCE_OFF_MASS_AMU=1.0 projects/jamil_replication/run_jamil_force_off_check.sh
```

Result:

- physical electron temperature: `489.111 K`
- all-mode extended temperature: `7825.777 K`
- expected extended mode temperature: `4800.000 K`
- centroid extended mode temperature: `2103.860 K`
- internal extended mode temperature: `8207.238 K`
- centroid kinetic-energy fraction: `0.0168`
- bead RMS spread: `0.006594 A`

This matches the one-step oracle P=16/dt=0.05/gamma=0.001 temperature, so the immediate one-step hot-electron diagnostic is reproducible without eReaxFF forces. It is not force-feedback heating on the first step.

Mass check:

- With `mass_amu=1.0`, the one-step force-off bead RMS spread is `0.006594 A`.
- With physical electron mass `mass_amu=0.00054858`, the same one-step force-off kinetic temperatures are unchanged, but bead RMS spread grows to `0.281545 A`.
- After 100 force-off steps with physical electron mass, bead RMS spread reaches `6.743825 A`.

Interpretation:

- The sampled kinetic temperature is mass-independent for the same RNG stream because the Maxwell velocity width scales as `1/sqrt(m)` and the diagnostic uses `m v^2`.
- The mass is still crucial for dynamics: physical electron mass produces much larger bead motion and is likely too delocalized/unstable for this eReaxFF electron coordinate unless the model was parameterized for it.
- The current `mass_amu=1.0` Jamil setting is therefore a fictitious electron mass convention, not the physical electron mass.
- The immediate `489 K` reading is mostly a finite-sample/initialization/internal-mode issue. With one PIMD electron and `P=16`, the kinetic diagnostic has only `48` velocity degrees of freedom, so a single-seed one-step temperature can fluctuate substantially.

## 2026-05-08 Partial Dispatch Sweep

I started the strict bounded-pool sweep with `JAMIL_DISPATCH_SKIP_W1=1` and `OMP_NUM_THREADS=1`.

Observed so far for `oracle_max_workers=2`:

- step 0 `initial_force` completed as 5 bead batches of 2 jobs each
- each 2-job batch was about `96-97 s`
- each bead call was still about `48-49 s`
- after 5 batches the run had already spent about `484 s` just on `initial_force`

That is the same negative shape as before: two-way concurrency did not improve throughput. The full sweep was stopped before `workers=4,6,8` because the `workers=2` result already showed that bounded process-level concurrency is not fixing the legacy backend cost.
