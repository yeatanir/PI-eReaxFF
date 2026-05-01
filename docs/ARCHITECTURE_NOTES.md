# Architecture Notes

The new engine follows an i-PI-like internal separation:

```text
configuration -> system/groups -> PIMD/thermostat/integrator -> force manager -> force provider
```

The eReaxFF adapter is one force provider.  Future i-PI socket support must
reuse that adapter rather than creating an independent force path.

## Current Force Layer

Implemented force-layer pieces:

- `force_provider_base.f90`: shared `force_result` container.
- `dummy_force_provider.f90`: zero energy and zero forces.
- `harmonic_force_provider.f90`: harmonic validation provider with
  `E = 1/2 k |x-x0|^2` and `F = -k(x-x0)`.
- `bead_force_manager.f90`: loops over beads, bead-averages the physical
  potential, and applies the central PIMD physical-force scaling `F/P`.
- `ereaxff_oracle_adapter.f90`: explicit not-implemented stub.

The force manager is where PIMD force scaling lives.  Drivers should not
duplicate `1/P` scaling.

## Current System/Group Layer

Implemented core pieces:

- `system_state.f90`: particle ids, elements, masses, positions, and velocities.
- `selectors.f90`: first config-style selector subset.
- `particle_groups.f90`: group specs and per-particle assignment.

Group assignment rejects unsupported selectors, invalid dynamics modes, invalid
bead counts, particles that match multiple groups, and particles that match no
group.

## Current Driver/Output Layer

Implemented pieces:

- `harmonic_pimd_main.f90`: thin harmonic PIMD executable using the validated
  normal-mode, BAOAB, and PILE-L modules.
- `pimd_output.f90`: centroid XYZ and bead-spread helpers.
- `diagnostics.f90`: conserved-like NVT diagnostic helper.
- `output_manager.f90`: path joining.
- `reaxff_compat_output.f90`: explicit not-implemented stub.

The harmonic driver is intentionally not a production eReaxFF driver.  It is the
first executable assembly test for the internal architecture.
