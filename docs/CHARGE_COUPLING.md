# Charge Coupling

Initial policy names:

- `bead_local_acks2`
- `centroid_acks2`
- `frozen_charges`
- `no_charge_coupling`

`bead_local_acks2` owns independent charge memory per bead.  `centroid_acks2`
is an explicit approximation and must be printed in diagnostics.

## Current Implementation

`charge_coupling_policy.f90` currently provides:

- policy-name validation;
- detection of policies that require bead-local memory;
- a small `charge_memory_state` with `charges_by_bead(nparticles, beads)`;
- per-bead charge updates;
- charge sum by bead;
- maximum charge difference across beads.

The current implementation is a validation scaffold for later ACKS2 coupling.
It deliberately does not call eReaxFF yet.
