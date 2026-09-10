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

