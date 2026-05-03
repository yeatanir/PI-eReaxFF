# Input Specification

The target production input format is YAML.  A simpler structured text parser
may be used temporarily during early Fortran development, but the user-facing
schema should retain these sections:

- `simulation`: timestep, steps, seed, integrator, output stride, restart
  stride;
- `force_provider`: provider type and eReaxFF input paths;
- `ring_polymer`: convention, bead count, propagator;
- `charge_coupling`: mode and tolerances;
- `geometry`: the validated phase4 graphene template, the user-supplied
  legacy `geo` loader, the electron starting position in Angstrom, and the
  explicit cell / periodicity contract;
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
- `geometry`: `template`, `structure_path`, `geo_path`, `periodic`,
  `crystx`, and `electron_position_A`;
- `groups`: `name`, `selector`, `dynamics`, `mobile`, `beads`, `thermostat`,
  `temperature_K`, `gamma_fs_inv`, `centroid_gamma_fs_inv`, and `mass_amu`;
- `output`: `centroid_xyz`, `beads_xyz`, `reaxff_compat`, `diagnostics`.

The parser is intentionally strict about structure and selector validity, but
it does not yet implement full YAML features such as anchors, nested maps, or
arbitrary sequences.  The `user_geo_structure` geometry contract is now
implemented for legacy `geo` files containing `HETATM` records and supports
both `structure_path` and `geo_path`.  `CRYSTX` records are parsed into the
shared cell state, and the YAML can also provide explicit `periodic` and
`crystx` values when the geometry file does not.

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

Current mixed-driver note: the classical carbon group now commonly uses a
Cartesian Langevin bath while the selected PIMD electron uses PILE-L.  The
driver prints the per-group thermostat choice at startup so mixed runs do not
silently inherit a single global bath.

Current multi-quantum note: the mixed BAOAB path now accepts multiple PIMD
particles in one run, provided they share a common bead count.  Their
pseudo-masses are specified per group with `mass_amu`, and the restart format
round-trips the full set of quantum indices, masses, normal-mode coordinates,
and cumulative heat.
