# Legacy Audit

## Scope

Inspected inputs:

- V1: `/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1`
- V2: `/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2`
- i-PI: `/Users/yeatanir/PhD/PIMD/i-pi`
- architecture report in the parent project directory.

## Reusable Code

- `V2/src/ereaxff_force_oracle.f`: the most important reusable boundary.  It is
  already framed as force/energy/charge evaluation without owning dynamics.
- `V2/src/phase4_single_point_oracle*.f90`: useful smoke-test patterns for
  static oracle evaluation.
- `V2/src/harmonic_pimd.f90`: useful benchmark intent, but the reusable math
  should be migrated into shared modules instead of copied.
- `V2/src/phase5_one_el_mobile_graphene_pimd.f90` and
  `V2/src/phase6_background_el_mobile_graphene_pimd.f90`: useful diagnostics,
  charge-memory patterns, and force scaling examples.
- `V1/src/debug/pimd_el_fixed_c.f`, `pimd_level2.f`, and `pimd_multiel.f`:
  useful historical references for normal-mode transform structure and electron
  diagnostics.

## Unsafe Duplicated Code

The following responsibilities are currently duplicated or embedded in
experiment-specific drivers:

- normal-mode matrix construction and transforms;
- ring-polymer frequency calculation;
- free ring-polymer propagation;
- PILE/Langevin coefficients;
- BAOAB ordering;
- bead force orchestration and `1/P` force scaling;
- charge-memory management;
- output/trajectory/diagnostic writing.

These must migrate into shared modules in the new tree.

## Ring-Polymer Convention Bug

Direct search and the architecture report identify:

- V2 `phase4_one_el_fixed_graphene.f90`: `omega_p = sqrt(P) * kB*T/hbar`;
- V2 `phase5_one_el_mobile_graphene_pimd.f90`: same;
- V2 `phase6_background_el_mobile_graphene_pimd.f90`: same;
- V2 `harmonic_pimd.f90`: same;
- V1 `pimd_el_fixed_c.f` and `pimd_multiel.f`: same.

The new code uses:

```text
omega_P = P * kB * T / hbar
omega_k = 2 * omega_P * sin(pi*k/P)
```

This matches the i-PI closed-ring convention described in
`ipi/engine/normalmodes.py` and the architecture report.

## Recommended Migration Path

1. Stabilize the shared normal-mode and ring-polymer modules with tests.
2. Add free-ring propagation and spring-energy tests before thermostats.
3. Add thermostat modules and heat accounting.
4. Add harmonic/dummy force providers and bead force manager.
5. Wrap V2 `ereaxff_force_oracle.f` through an adapter without editing legacy
   source in place.
6. Recreate legacy phase4/phase5/phase6 scenarios as configuration files and
   regression tests, not as new hard-coded drivers.

