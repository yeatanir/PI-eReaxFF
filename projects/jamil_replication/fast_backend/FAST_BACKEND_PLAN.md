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
- Full dynamic ACKS2 matrix solving and physics energy/force evaluation remain
  fail-closed until their parity kernels are implemented.  The next force work
  should add `boncor` bond-order corrections and then valence/overcoordination
  terms against oracle labels.
