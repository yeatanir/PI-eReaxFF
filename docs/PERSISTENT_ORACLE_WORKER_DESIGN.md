# Persistent Oracle Worker Design

## Problem Statement

The current oracle backend launches a fresh legacy eReaxFF binary for each force evaluation. That is correct but expensive. On the Jamil 112C/113El P=16 case, the dominant cost is not the PIMD integrator; it is repeated external process startup, file staging, parsing, and teardown.

The bounded dispatch refactor makes this visible:

- `serial` is stable but slow.
- `openmp_all` is functionally correct but performs poorly under concurrent legacy launches.
- `bounded_pool` gives controlled concurrency, but it still pays the full process-launch cost for every bead-level oracle call.

The next backend step is to keep the legacy oracle process alive across calls.

## Goal

Replace repeated process launch with persistent workers that:

- initialize the legacy backend once
- keep a private scratch/work directory
- service repeated force requests
- return energy, forces, charges, and solver state
- preserve bead-local charge-memory behavior
- remain compatible with the current driver and timing/logging model

## Non-Goals

- Do not change the physics.
- Do not change geometry setup or charge assignment.
- Do not remove the current external-process backend.
- Do not introduce asynchronous complexity before the synchronous persistent path is working.

## Proposed Backend Shape

Add a new backend implementation alongside the current external-process adapter:

- `external_process_oracle_backend`
- `persistent_worker_oracle_backend`

The driver and bead-force manager should talk to a small backend interface, not to the legacy binary directly.

## Worker Lifecycle

Each worker should:

1. start once
2. copy or initialize its own legacy inputs
3. wait for requests
4. evaluate one request at a time
5. write outputs to its private scratch tree
6. return results
7. remain alive until shutdown

For the first implementation, one worker can simply be a long-lived subprocess with a line-oriented stdin/stdout protocol or a small control file protocol. A socket protocol can come later.

## Scratch Directory Ownership

Each worker must own a distinct scratch root:

- worker 1: `.../worker_0001`
- worker 2: `.../worker_0002`
- etc.

Each request should get its own call subdirectory inside that worker root.

This avoids collisions from:

- geometry files
- `control` copies
- restart/charge files
- force export outputs
- timing logs

## Request / Response

Keep the first protocol simple and synchronous.

Request fields:

- request id
- bead index
- step index
- stage label
- positions
- current charges
- solver guess
- charge policy
- optional metadata for logging

Response fields:

- request id
- status
- energy
- forces
- charges
- solver guess
- diagnostic text

The request id should map cleanly to the current `call_slot` concept.

## Charge-Memory Handling

Charge-memory updates should remain on the driver side, after all bead-local worker calls complete.

That means:

- workers may compute bead-local charges
- workers must not mutate shared ring-polymer charge memory
- the bead-force manager commits charge-memory updates in deterministic bead order after the batch finishes

This preserves the current race-free behavior.

## Timing and Logging

The current timing model should stay:

- per-call timing records
- per-batch summaries
- stage labels
- batch summaries in `oracle_timing.log`

Persistent workers should record:

- worker id
- request id
- stage label
- request start/end
- elapsed time
- status
- worker scratch path

The first implementation can keep the same log format and just substitute the worker backend.

## Failure and Timeout Handling

The persistent backend needs bounded failure behavior:

- if a worker fails to start, fail early
- if a request returns nonzero status, preserve current fail-fast behavior
- if a request hangs, enforce a timeout and restart or fail the worker
- if a worker crashes, mark its pending request failed

For now, synchronous request handling is enough. No speculative retry logic is needed in the first pass.

## Legacy Fortran State Risk

The main risk is hidden global state in the legacy backend:

- static arrays
- fixed filenames
- implicit working-directory assumptions
- solver state carried across calls
- outputs that depend on file system side effects

This is why the worker must keep a private scratch tree and why the first implementation should be conservative:

- one worker process per scratch root
- one request at a time per worker
- explicit reset or reinitialization hooks if the legacy code needs them

## Minimal Implementation Plan

1. Extract a small backend interface from the current adapter.
2. Implement `persistent_worker_oracle_backend` as a synchronous worker loop.
3. Keep the old external-process path available for regression and fallback.
4. Run the current regression test suite unchanged.
5. Benchmark Jamil P=16 with `serial`, `bounded_pool`, and persistent worker mode.
6. Only after that, consider a socket protocol or fully asynchronous scheduling.

## Expected Benefit

If process startup and file staging dominate, persistent workers should reduce wall time more than any worker-count tweak.

In the best case, the expensive loop becomes:

- initialize once
- evaluate many requests
- shut down once

instead of:

- launch binary
- stage files
- parse control
- compute
- tear down
- repeat 80 times per step
