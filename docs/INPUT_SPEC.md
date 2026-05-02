# Input Specification

The target production input format is YAML.  A simpler structured text parser
may be used temporarily during early Fortran development, but the user-facing
schema should retain these sections:

- `simulation`: timestep, steps, seed, integrator, output stride, restart
  stride;
- `force_provider`: provider type and eReaxFF input paths;
- `ring_polymer`: convention, bead count, propagator;
- `charge_coupling`: mode and tolerances;
- `groups`: selectors, dynamics type, mobility, thermostat, temperature, and
  bead count;
- `output`: PIMD-native and compatibility output choices.

Current implementation note: the eReaxFF adapter still uses legacy V1/V2
fixture paths for `control`, `ffield`, and `iopt`, and it generates the phase4
geometry through the validated legacy helper in its work directory.  That is
intentional for the current validation stage; the user-facing YAML path is the
next layer to standardize once the mixed dynamics boundary is stable.

## Current Selector Subset

The current Fortran selector scaffold supports:

- `all`
- `element == C`
- `element == El`
- `id == 61`
- `id != 61`
- `id in [61,64]`
- one `and` combining two supported clauses, for example
  `element == El and id == 61`

This is intentionally small but already covers the first target scenarios:
classical carbons, one selected PIMD electron, and frozen background electrons.
Unsupported selector syntax fails during group validation.
