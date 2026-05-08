# Jamil eReaxFF Graphitic Anode Replication Plan

Goal: reproduce the classical eReaxFF graphitic-anode results from Hossain,
Pawar, and van Duin (JES 2022), then use this codebase to add controlled PIMD
electron extensions.

Paper reference: `/Users/yeatanir/PhD/PIMD/Jamil_paper.pdf`

## Working Assumptions

- The paper simulations are classical eReaxFF, not PIMD.
- Classical paper setup uses C+/electron pairs, pseudo-electron mass `1.0 amu`,
  velocity Verlet, Berendsen thermostat with `tau = 100 fs`, `dt = 0.1 fs`,
  and mostly `300 K` NVT trajectories.
- PIMD is an extension we will add after reproducing the classical behavior.
- Guard rails stay enabled during validation: repeatability checks, trace output,
  restart output, and short stress tests remain available until the classical
  baseline is trusted.

## Phase 0 - Paper/Artifact Inventory

1. Extract the paper methods and target observables.
2. Locate or reconstruct the exact structures:
   - single-layer graphene with one electron per C;
   - single-layer graphene with one extra electron;
   - four-layer graphite, 32 C per layer, one or more excess electrons;
   - graphene vacancy system;
   - Li/Li+ graphene systems, later.
3. Locate the force-field/control files used by the current V2 Track B2 setup
   and compare them with the paper supporting-information force field.
4. Decide which target is first:
   - recommended: Fig. 3-style in-plane electron mobility on graphene.

Acceptance:
- A table exists listing each target figure/result, required structure, particle
  count, temperature, duration, and observable.

## Phase 1 - Classical Jamil Baseline in This Driver

1. Add a classical-like Jamil config:
   - `P = 1`;
   - pseudo-electron mass `1.0 amu`;
   - `dt_fs = 0.1`;
   - `T = 300 K`;
   - carbon thermostat `berendsen`, `tau_fs = 100`;
   - electron thermostat either Berendsen-compatible classical mode or a
     documented P=1 PIMD stand-in until true classical electron groups are added.
2. Start with the existing 60C/60El paired graphene geometry.
3. Add one extra electron to make a charged system equivalent to the paper's
   "one excess electron" experiment.
4. Run short guard-rail trajectories:
   - 0-step initialization;
   - 1 step;
   - 10 steps;
   - 100 steps.
5. Record timing:
   - wall time;
   - number of oracle calls;
   - seconds per oracle call;
   - output directory size.

Acceptance:
- The run is stable through 100 steps.
- Total charge and charge localization diagnostics are finite and consistent.
- PBC remains present in generated oracle geos.
- Timing is recorded in a small benchmark summary file.

## Phase 2 - Electron Mobility Observable

1. Implement an electron-host assignment analysis:
   - for every output frame, assign each electron to nearest carbon under PBC;
   - identify the excess electron by charge/localization criterion;
   - track visited host carbon IDs.
2. Add path-density output:
   - host visit counts;
   - path length or hop count;
   - per-frame nearest C for the excess electron.
3. Reproduce a Fig. 3-style observable:
   - excess electron visits many C atoms in-plane over time;
   - no systematic out-of-plane drift for single-layer graphene.

Acceptance:
- A script produces a CSV summary and an XYZ/visualization artifact for the
  excess electron path.

## Phase 3 - Match Paper System Sizes

1. Generate a single-layer graphene cell closer to the paper:
   - neutral: one electron per C;
   - charged: one extra electron per 112 C.
2. Validate the generated geometry:
   - C-C distances;
   - electron initial radius around C, about `0.35 A`;
   - CRYSTX/PBC;
   - no duplicate or near-zero distances.
3. Run short tests at `300 K`.
4. Extend to a modest trajectory once stable.

Acceptance:
- The generated geometry passes a geometry-check script.
- The classical baseline reproduces qualitative in-plane mobility.

## Phase 4 - Graphite Out-Of-Plane Conductivity

1. Build four-layer graphite:
   - 32 C per layer;
   - one electron per C;
   - one excess electron in layer 2.
2. Run 300 K, 400 K, and 500 K short trajectories.
3. Track layer assignment of the excess electron.
4. Count out-of-plane electron events.

Acceptance:
- At 300 K, excess electron remains mostly in-plane.
- Higher temperature runs show increased out-of-plane transitions, at least
  qualitatively.

## Phase 5 - Vacancy Localization

1. Generate graphene with one C vacancy.
2. Place one excess electron far from the vacancy.
3. Track whether the excess electron localizes near the vacancy and nearest
   neighbors.
4. Compare against paper's reported localization metric.

Acceptance:
- The excess electron spends elevated time near vacancy/neighbor atoms.
- A localization percentage is reported.

## Phase 6 - PIMD Extension

1. Convert the validated classical setup to electron PIMD:
   - begin with `P = 4`;
   - use explicit electron mass choices:
     - paper-like `1.0 amu` control;
     - smaller pseudo mass;
     - real electron mass only after stability checks.
2. Compare against classical:
   - electron path density;
   - localization near vacancy;
   - ring-polymer spread;
   - thermostat heat and temperature behavior.
3. Scale to `P = 8`, then `P = 16`.

Acceptance:
- PIMD runs remain stable for short trajectories.
- PIMD observables can be compared directly to classical Jamil-style runs.

## Phase 7 - Performance/Guard-Rail Strategy

1. Keep strict mode for stress tests:
   - repeatability enabled;
   - trace output retained;
   - restart output frequent.
2. Add measured timing summaries for strict and faster modes.
3. Only after the physics baseline is trusted, consider production cleanup:
   - disable V2 `electron_trace.out`;
   - avoid keeping every adapter work directory;
   - move toward a persistent oracle interface.

Acceptance:
- Every reported physics run includes its guard-rail mode and timing metadata.

## Immediate Next Tasks

1. Build `projects/jamil_replication/targets.md` from the paper.
2. Add a geometry generator for the first single-layer graphene target.
3. Add a P=1 classical-style config using pseudo-electron mass `1.0 amu`.
4. Add electron-host assignment analysis.
5. Run 0-step and 1-step strict tests and record timing.
