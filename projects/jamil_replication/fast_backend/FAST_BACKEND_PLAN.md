# Fast eReaxFF/RPMD Implementation Plan

## Target

Reach native-ReaxFF-like throughput for the Jamil 112C/113El RPMD case without
removing the exact legacy oracle.

## Near-Term Milestones

1. Build reference cases from existing exact Jamil runs.
2. Add a batched backend protocol that can submit all beads for one force stage
   in one request.
3. Implement a JAX backend for the Jamil force-field subset and compare every
   output to the oracle dataset.
4. Add force-stage caching and remove duplicate force calls only where the
   trajectory remains equivalent under the validation suite.

## Acceptance

- Exact oracle remains the default backend.
- No fast backend can run production MD until its parity test passes.
- First useful target: P=16 one-step Jamil wall time below 60 seconds.
  Achieved for the connected JAX path on 2026-09-10: `12.745 s` for
  `P=16`, one 0.1 fs step, 5 force stages, and 80 bead-force cases.
- Long target: P=16 Jamil dynamics within roughly 2-5x the same-system
  classical ReaxFF cost.

## Current Implementation Slice

- Adapter-work JSONL extraction now produces labeled oracle cases from
  `geo`, `charge_memory.in`, `phase4_single_point_oracle.out`, and
  `single_point_forces.out`.
- The JAX entrypoint can run a labeled dataset readiness check using x64 JAX
  arrays.
- The first JAX charge kernel supports the actual Jamil `icharg=8`
  fixed-`MOLCHARGE` path and reports charge parity against oracle labels.
- A strict oracle-replay force kernel now validates batched response shape and
  parity plumbing, while failing closed unless replay is explicitly requested.
- An experimental fixed-charge Coulomb-only JAX kernel now exercises real
  energy/force autodiff and batching behind `--allow-unvalidated`.
- An experimental source-derived shielded/tapered eReaxFF Coulomb slice now
  parses the oracle call `ffield` and uses the Jamil `poten.f`/`reac.f`
  formula for the Coulomb term.
- An experimental source-derived no-inner-wall, shielded-vdW nonbonded kernel
  now combines the Jamil `nonbonnoinS` vdW and Coulomb terms in a batched JAX
  path.
- An experimental source-derived raw bond-order and uncorrected covalent
  bond-energy kernel now parses bond/off-diagonal `ffield` sections and runs
  in a batched JAX path.
- An experimental source-derived corrected bond-order and covalent bond-energy
  kernel now applies the Jamil `boncor` overcoordination and 1-3 BO correction
  formulas before the `covbon` energy expression.
- An experimental aggregate kernel now sums the implemented corrected-bond,
  shielded-vdW, and shielded/tapered-Coulomb slices in one batched JAX call.
- The mixed BAOAB driver can now launch the persistent JAX worker through
  `force_provider.type: fast_jax_oracle`, guarded by
  `fast_backend_allow_unvalidated: true`, and batch all beads in one force
  request per force stage.
- The connected fast path is about `12.1x` faster than the exact P=4
  persistent oracle pool for the same Jamil one-step schedule, with near-flat
  wall time from `P=4` to `P=16`.
- Short P=16 fast MD checks pass finite/output/worker-shutdown guards for
  10 steps at `dt_fs=0.05` and `0.1`, but the 100-step `dt_fs=0.05` run is not
  thermally acceptable even after force-parity fixes: carbon temperature rises
  to about `2136 K` over 5 fs.  The matching 100-step `dt_fs=0.01` run passes
  the thermal guard with final carbon temperature about `400 K`.
- Worker timing shows request parsing/writing is negligible; the current
  bottleneck is the JAX force kernel itself at roughly `1.9-2.1 s` per P=16
  force-stage request after startup.
- Off-equilibrium oracle parity now labels selected fast trajectory frames.
  The near-equilibrium 10-step subset stays below `1e-5 eV/A`, and the
  thermally bad 100-step sampled trajectory now stays below `1e-6 eV/A` after
  fixing the electron-sharing convention and the legacy C-C C2 correction in
  the `lone_pair` bucket.
- Full dynamic ACKS2 matrix construction/solve remains fail-closed until it has
  oracle labels that exercise nonzero solver state.  The connected Jamil force
  path is still explicitly unvalidated and should remain behind the opt-in gate
  while off-equilibrium frame parity, drift diagnostics, and profiling are
  added.

## Next Correctness Target

Map the stable timestep window before optimizing speed.  The corrected force
kernel passes broad off-equilibrium parity and the `dt_fs=0.01` stability check,
but `dt_fs=0.05` still overheats carbon.  Run intermediate `dt_fs` values
(`0.02`, then `0.025` or `0.03`) and inspect conserved-energy drift to choose a
usable production timestep or identify the integrator/scaling term that makes
`0.05 fs` unstable.
