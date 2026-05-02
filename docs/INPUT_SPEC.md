# Input Specification

The target production input format is YAML.  A simpler structured text parser
may be used temporarily during early Fortran development, but the user-facing
schema should retain these sections:

- `simulation`: timestep, steps, seed, integrator, output stride, restart
  stride;
- `force_provider`: provider type and eReaxFF input paths;
- `ring_polymer`: convention, bead count, propagator;
- `charge_coupling`: mode and tolerances;
- `geometry`: currently the validated phase4 graphene template plus the
  electron starting position in Angstrom;
- `groups`: selectors, dynamics type, mobility, thermostat, temperature, and
  bead count;
- `output`: PIMD-native and compatibility output choices.

Current implementation note: the eReaxFF adapter still uses legacy V1/V2
fixture paths for `control`, `ffield`, and `iopt`, and it generates the phase4
geometry through the validated legacy helper in its work directory.  The mixed
BAOAB driver now reads a YAML-subset config file directly, but the legacy
oracle inputs remain explicit paths relative to the V1 root.  That is
intentional for the current validation stage; the user-facing YAML path is the
next layer to standardize once the mixed dynamics boundary is stable.

## Current YAML-Subset Parser

The Fortran input parser currently accepts:

- `simulation`: `dt_fs`, `steps`, `seed`, `integrator`, `output_stride`,
  `restart_stride`, `output_dir`;
- `force_provider`: `type`, `v1_path`, `v2_path`, `work_root`, `ffield`,
  `geo`, `control`, `iopt`, `smoke_binary`, `repeat_binary`,
  `make_geo_script`, `timeout_runner`;
- `ring_polymer`: `convention`, `default_beads`, `propagator`;
- `charge_coupling`: `mode`, `charge_tolerance`;
- `groups`: `name`, `selector`, `dynamics`, `mobile`, `beads`, `thermostat`,
  `temperature_K`, plus currently ignored thermostat tuning keys such as
  `gamma_fs_inv` and `centroid_gamma_fs_inv`;
- `output`: `centroid_xyz`, `beads_xyz`, `reaxff_compat`, `diagnostics`.

The parser is intentionally strict about structure and selector validity, but
it does not yet implement full YAML features such as anchors, nested maps, or
arbitrary sequences.

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
