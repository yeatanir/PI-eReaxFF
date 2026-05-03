# Multi-Particle PIMD Architecture Plan

## Why this exists

The current mixed eReaxFF BAOAB path was built as a validation scaffold around
one selected PIMD electron plus a classical substrate. That was enough to prove
the modular boundary, thermostat plumbing, restart/output handling, and the
`omega_P = P k_B T / \hbar` ring-polymer convention.

It is not enough for the real problem class.

The target workload includes:

- one or more excess electrons as PIMD particles;
- classical carbon nuclei with `nbeads = 1`;
- frozen particles;
- multiple quantum groups in the same system;
- explicit cell and periodicity handling;
- explicit charge-coupling policy per group or per subset.

The architecture must stop assuming:

- exactly one PIMD particle;
- exactly one PIMD group;
- one electron per system;
- one ring-polymer bead count for every particle;
- one thermostat mode for all particles.

## Scientific model

The intended model is mixed quantization:

- carbon nuclei stay classical unless a future test explicitly promotes them;
- explicit electrons can be classical or PIMD, depending on input;
- PIMD is a property of a particle group, not a hard-coded special case;
- the force oracle remains a force/charge evaluator only;
- PIMD beads live entirely outside the oracle boundary.

For the current eReaxFF graphitic systems, the realistic near-term model is:

- carbon atoms: `nbeads = 1`;
- selected excess electrons: `nbeads = P`;
- background electrons: classical, frozen, or quantum only when explicitly
  requested.

## Main architectural change

Replace the current single-quantum-particle assumption with a particle-group
layout that can carry multiple quantum subsets at once.

The core idea is:

- particle groups define dynamics, bead count, thermostat, temperature, and
  mobility;
- the system state holds per-particle metadata and explicit cell state;
- the integrator operates on groups;
- the force manager evaluates per group and per bead according to policy;
- the output layer writes centroid-projected compatibility views plus PIMD-
  native bead diagnostics.

## Required data-model changes

### 1. Per-particle bead count must be first-class

The current assignment layer already stores `beads(i)`, but the driver path
still uses a single selected PIMD particle in practice.

Next state:

- every particle carries a group assignment;
- every group carries its own bead count;
- the driver builds one ring polymer per PIMD particle or per PIMD group
  member;
- classical particles remain single-copy states;
- frozen particles remain excluded from propagation and thermostatting.

### 2. Multiple PIMD particles must be represented explicitly

The current `pimd_count /= 1` guard must be removed.

The new layout should support:

- zero, one, or many PIMD particles;
- different PIMD bead counts in the same system;
- multiple electrons with independent ring polymers;
- mixed classical/PIMD/frozen subsets in one run.

### 3. Charge state must be bead-aware and particle-aware

The current `charge_memory_state` is bead-local, which is correct for
bead-local ACKS2, but it must extend naturally to multiple quantum particles.

Required contract:

- charge memory is indexed by particle and bead;
- centroid contraction remains explicit and opt-in;
- each quantum group can choose bead-local, centroid, frozen, or no coupling;
- charge diagnostics must report per-particle and per-bead spread.

### 4. Cell handling remains shared

Cell state must stay in `system_state`, not in the oracle or thermostat layer.

Required:

- explicit `periodic` flag;
- explicit `crystx` or cell matrix;
- shared minimum-image helper for all connectivity and analysis code.

## Force-manager changes

The force manager needs to evolve from a single-electron bead loop into a
general particle-group evaluator.

Target behavior:

- evaluate classical particles once per step;
- evaluate each PIMD particle or PIMD group member across its bead set;
- keep physical force scaling centralized;
- keep charge policy explicit per subset;
- preserve bead-local charge memory for ACKS2;
- average only where the physics requires averaging;
- do not smear all quantum particles into one shared ring polymer.

The manager must eventually expose a structure like:

```text
evaluate_particles(state, groups, positions, charge_state, policy) -> summary
```

where `summary` contains:

- energy;
- per-particle forces;
- per-particle charges;
- per-bead diagnostics;
- charge convergence status;
- thermostat heat bookkeeping;
- error/status information.

## Integrator changes

The integrator must stop assuming one quantum electron and one thermostat.

Required shift:

- the integrator loops over particle groups, not hard-coded species;
- each group chooses its own thermostat object;
- classical and PIMD update paths are both supported in one step;
- multiple quantum particles do not share bead coordinates by accident;
- BAOAB/OBABO ordering remains defined in one place.

The first implementation target should be:

- classical carbon nuclei with `nbeads = 1`;
- multiple excess electrons each with `nbeads = P`;
- a single mixed BAOAB step that advances all groups consistently.

## Input contract changes

The YAML contract must be expanded beyond one electron.

Required additions:

- group-level bead count remains mandatory;
- group-level dynamics must support multiple PIMD groups;
- group-level electron mass or pseudo-mass must be configurable;
- group-level thermostat and temperature remain separate;
- group-level charge policy must be allowed where needed.

For the electron model specifically:

- the pseudo-particle mass must be explicit and configurable;
- it must not be hard-coded as a special value in the driver;
- the test harness must be able to switch between Jamil-style pseudo-electron
  mass and the current validation mass where needed.

## Output changes

The output layer must become multi-particle aware.

Needed outputs:

- centroid positions per PIMD particle or per quantum group;
- bead positions for each quantum particle;
- gyration radius per quantum particle;
- layer assignment / defect proximity / localization metrics per particle;
- per-particle and per-group charge summaries;
- per-group thermostat heat;
- compatibility outputs projected from the centroid view only when requested.

## Testing strategy

The test suite needs to shift from single-electron regression to a matrix of
particle counts and group patterns.

### Phase A: configuration and model checks

- parse multiple PIMD groups;
- allow more than one quantum electron;
- verify group-specific bead counts;
- verify explicit pseudo-electron mass selection;
- reject accidental `pimd_count == 1` assumptions.

### Phase B: force-manager checks

- multiple quantum particles in one force evaluation;
- bead-local charge memory for each quantum particle;
- centroid contraction only where explicitly configured;
- classical carbon forces remain `P=1`.

### Phase C: integrator checks

- one-step mixed update with multiple quantum particles;
- thermostat heat accounting per group;
- restart continuity with multiple quantum particles;
- deterministic split-run equivalence.

### Phase D: physics-regression checks

- pristine graphene with one excess electron;
- graphite with one excess electron in a selected layer;
- bilayer graphene with linker distortion and two excess electrons;
- vacancy localization case;
- per-particle gyration radius and layer occupancy.

## Implementation order

1. Remove the hard `pimd_count == 1` assumption from the driver path.
2. Extend the group assignment and state plumbing to support many PIMD
   particles at once.
3. Make pseudo-electron mass a per-group configuration value.
4. Generalize the force-manager loops to operate on multiple quantum particles.
5. Generalize the mixed integrator to propagate all PIMD groups.
6. Extend restart/output to carry multiple quantum particles cleanly.
7. Add the new replication/physics regression matrix.

## What not to do

- Do not add another single-purpose `phaseX` driver.
- Do not share one ring polymer across unrelated electrons.
- Do not let the force oracle own bead propagation or thermostatting.
- Do not hide centroid contraction as if it were full PIMD.
- Do not leave the pseudo-electron mass implicit.
- Do not preserve `pimd_count == 1` in the final production path.

## Current boundary

The current repository state is still valid as a validation scaffold, but the
multi-particle architecture above is the next necessary step before the code
can represent the full family of Jamil-style graphitic electron transport
problems.
