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

## 2026-09-10 Initial PIMD Velocity Rescale

Added `simulation.rescale_initial_pimd_velocities`. When enabled for a fresh run, the mixed BAOAB driver rescales the initial PIMD normal-mode velocities only once:

- centroid mode kinetic energy is set to its target extended-mode temperature
- internal modes `2:P` are collectively set to their target extended-mode temperature
- the thermostat is not clamped after initialization, so subsequent stochastic PILE-L fluctuations remain canonical

Force-off validation with `P=16`, `dt_fs=0.05`, `mass_amu=1.0`:

- step 0 physical electron temperature: `300.000 K`
- step 0 all-mode extended temperature: `4800.000 K`
- step 0 centroid extended temperature: `4800.000 K`
- step 0 internal extended temperature: `4800.000 K`
- 10000-step physical electron temperature average: `299.207 K`
- last-5000-step physical electron temperature average: `301.249 K`

Jamil oracle one-step check:

```bash
JAMIL_P_SWEEP_CASES=16 JAMIL_P_SWEEP_GAMMAS=0.001 JAMIL_P_SWEEP_DT_FS=0.05 \
  projects/jamil_replication/run_jamil_p_sweep.sh
```

Result:

- status: `PASS`
- physical electron temperature: `368.562 K`
- all-mode extended temperature: `5896.991 K`
- expected extended temperature: `4800.000 K`
- centroid extended temperature: `4772.871 K`
- internal extended temperature: `5971.932 K`
- carbon temperature: `292.972 K`
- bead RMS spread: `0.005462 A`
- oracle calls: `80`
- wall time: `449.718 s`

This fixes the hot initial draw at step 0 and improves the one-step P=16 oracle diagnostic by about `120.5 K` versus the previous `489.111 K`. The remaining one-step overshoot comes from the stochastic PILE-L update, not from initial kinetic-energy partitioning or eReaxFF force feedback.

## 2026-05-08 Partial Dispatch Sweep

I started the strict bounded-pool sweep with `JAMIL_DISPATCH_SKIP_W1=1` and `OMP_NUM_THREADS=1`.

Observed so far for `oracle_max_workers=2`:

- step 0 `initial_force` completed as 5 bead batches of 2 jobs each
- each 2-job batch was about `96-97 s`
- each bead call was still about `48-49 s`
- after 5 batches the run had already spent about `484 s` just on `initial_force`

That is the same negative shape as before: two-way concurrency did not improve throughput. The full sweep was stopped before `workers=4,6,8` because the `workers=2` result already showed that bounded process-level concurrency is not fixing the legacy backend cost.

## 2026-09-10 Fast JAX Mixed-Driver Timing

Added `projects/jamil_replication/run_fast_vs_oracle_timing.sh` to compare the connected fast JAX driver path against the exact persistent oracle pool on the same Jamil 112C/113El one-step BAOAB schedule.

Command:

```bash
projects/jamil_replication/run_fast_vs_oracle_timing.sh
```

Default benchmark set:

- fast JAX: `P=4`, `P=8`
- exact oracle pool: `P=4`
- `steps=1`, `dt_fs=0.1`, `OMP_NUM_THREADS=4`, `oracle_max_workers=4`

Results:

| backend | P | steps | wall_s | force stages | bead force cases | s/stage | s/bead case | electron_T_K | extended_T_K | carbon_T_K | energy_eV | status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| fast JAX | 4 | 1 | 12.139 | 5 | 20 | 2.428 | 0.607 | 295.036 | 1180.146 | 293.429 | -510.332554 | PASS |
| fast JAX | 8 | 1 | 12.003 | 5 | 40 | 2.401 | 0.300 | 297.356 | 2378.850 | 293.475 | -517.060030 | PASS |
| oracle pool | 4 | 1 | 146.497 | 5 | 20 | 29.299 | 7.325 | 295.036 | 1180.146 | 293.429 | -510.332544 | PASS |

Fast P=16 was run separately with the oracle leg skipped:

```bash
JAMIL_FAST_TIMING_FAST_CASES=16 JAMIL_FAST_TIMING_ORACLE_CASES= \
  projects/jamil_replication/run_fast_vs_oracle_timing.sh
```

Result:

| backend | P | steps | wall_s | force stages | bead force cases | s/stage | s/bead case | electron_T_K | extended_T_K | carbon_T_K | energy_eV | status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| fast JAX | 16 | 1 | 12.745 | 5 | 80 | 2.549 | 0.159 | 400.393 | 6406.292 | 293.394 | -511.895829 | PASS |

Interpretation:

- The connected fast path is already about `12.1x` faster than the exact P=4 oracle pool for the same one-step schedule.
- Fast JAX wall time is nearly flat from `P=4` to `P=16`, which means batching all beads through one JAX call is doing the right thing.
- The P=4 fast/oracle physical diagnostics match closely; the scalar potential differs by about `9.9e-6 eV`, consistent with the current force-kernel parity gate.
- P=16 fast clears the near-term target from `FAST_BACKEND_PLAN.md`: one-step Jamil wall time is below `60 s`.
- The next speed target is not process-level parallelism. It is profiling and simplifying the JAX kernel itself, then removing redundant force stages only if trajectory-equivalence tests allow it.

## 2026-09-10 Fast JAX Stability Checks

Added `projects/jamil_replication/run_fast_stability_check.sh` and the short regression target `make test_33_fast_jax_stability`.  The script generates guarded `fast_jax_oracle` configs, runs the mixed BAOAB driver, checks finite diagnostics and `.xyz` outputs, verifies worker shutdown, and records per-request worker timing from `fast_jax_work/worker_timing.csv`.

Short P=16 checks:

| P | steps | dt_fs | wall_s | s/step | s/force stage | worker force sum_s | electron_T_K | extended_T_K | carbon_T_K | bead spread_A | conserved span_eV | status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 16 | 10 | 0.05 | 83.530 | 8.353 | 2.037 | 80.513 | 283.339 | 4533.430 | 317.117 | 0.050581 | 88.720 | PASS |
| 16 | 10 | 0.10 | 83.158 | 8.316 | 2.028 | 80.173 | 305.601 | 4889.609 | 395.333 | 0.090160 | 276.077 | PASS |

Medium P=16 check:

| P | steps | dt_fs | wall_s | s/step | s/force stage | worker force sum_s | electron_T_K | extended_T_K | carbon_T_K | bead spread_A | conserved span_eV | status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 16 | 100 | 0.05 | 848.873 | 8.489 | 2.117 | 820.014 | 251.321 | 4021.139 | 2138.742 | 0.148771 | 1038.979 | THERMAL FAIL |
| 16 | 100 | 0.01 | 816.429 | 8.164 | 2.036 | 789.258 | 205.324 | 3285.179 | 399.733 | 0.080271 | 303.617 | PASS |

Interpretation:

- The fast driver path is numerically robust enough for short P=16 MD: no NaN/Inf diagnostics, `.xyz` trajectories are written, and the worker exits cleanly.
- Worker timing shows the bottleneck clearly: the 10-step run spends about `80.5 s` of `83.5 s` in JAX force evaluation. Request parsing/writing is below `0.4 s` total, so protocol overhead is not the main limiter.
- The 100-step `dt=0.05` trajectory is not acceptable production MD yet. Carbon temperature rises monotonically to about `2139 K` over 5 fs while electron temperature remains near target. This originally pointed to force/integrator drift or incomplete off-equilibrium force parity, not the earlier hot-electron initialization issue.

2026-09-11 corrected-force follow-up:

- After fixing off-equilibrium force parity, the `dt=0.05 fs` 100-step run still fails the thermal guard: final carbon temperature is `2136.095 K`, final physical electron temperature is `251.373 K`, and conserved-like span is `1038.982 eV`.
- The same corrected-force 100-step run at `dt=0.01 fs` passes: final carbon temperature is `399.733 K`, physical electron temperature is `205.324 K`, bead spread is `0.080271 A`, and conserved-like span is `303.617 eV`.
- The carbon heating is therefore timestep/integrator stability, not the sampled fast-vs-oracle force parity issue.  Next check should map the stable window with intermediate timesteps such as `0.02 fs` and `0.025 fs`.

## 2026-09-10 Fast Off-Equilibrium Oracle Parity

Added `projects/jamil_replication/run_fast_offequilibrium_parity.sh` and `make test_34_fast_offequilibrium_parity`.  The script samples selected bead frames from a fast `.xyz` trajectory, builds adapter-call-like exact oracle inputs, labels those frames with `oracle_single_point_force_export`, converts them to JSONL, and runs the JAX force-term comparison.

Regression subset from the P=16 10-step `dt=0.05` trajectory:

- sampled frames: `step=0,10`, `bead=1`
- cases: `2`
- max force error: `3.73e-7 eV/A`
- status: `PASS` with threshold `1e-5 eV/A`

Broader diagnostic subset from the thermally bad P=16 100-step `dt=0.05` trajectory:

- sampled frames: `step=0,50,100`, `beads=1,8,16`
- cases: `9`
- max force error before electron-sharing fix: `3.81 eV/A`
- worst frames: `step=100`, atom `114` (`C`), mostly `x/y` components
- step 50 is still close: worst force error about `5.36e-3 eV/A`

Interpretation:

- The fast force kernel is accurate near the initial and short-trajectory manifold.
- The 100-step carbon heating correlates with a real off-equilibrium force mismatch, not just thermostat noise.

Follow-up correction:

- `electron_corrected_valencies` now matches the Jamil legacy central-cell electron-sharing convention used by the electron/nonbonded buckets.  This removes the periodic-image overcoordination error for `step=100`, `bead=8` at atom `114`.
- After that fix, the broader diagnostic max force error was reduced to `3.61 eV/A`, with no NaNs in the force output.
- `finite_difference_oracle_force.py` showed the smooth negative-side oracle-energy derivative matched the exported legacy force, while JAX was missing the remaining residual in the `lone_pair` bucket.  The missing source term was the legacy C-C C2 correction accumulated into `elp`.
- After adding the C2 correction, the same broad 9-frame diagnostic reaches max force error `6.59e-7 eV/A`; the previous atom-114 `step=100` mismatch is gone.
