# Adopting the "smooth" control file for runs with mobile atoms — 2026-11-11

Decision (coordinator): legacy thresholds stay the default for all reference comparisons; a second control file, `smooth`, changes BOTH the bond-list cutoff and the 3-/4-body BO threshold; no new code flag. Smooth is for runs with mobile atoms, legacy for validation.
Labels: **[measured]**, **[derived]**.

## The two files (`project/inputs/control/`)
| file | nbrhood_cutoff (bond list, Å) | thb_cutoff (BO threshold, 3-/4-body) | thb_cutoff_sq |
|---|---|---|---|
| `legacy.ctrl` (identical to `pair_style ... NULL`, verified: 0.0 difference in E, F and all terms on 150 frames) | 5.0 | 0.001 | 1e-5 |
| `smooth.ctrl` | 5.2 | 1e-9 | 1e-18 |
The control parser accepts no comments, so the explanation is in `inputs/control/README.md`. `run_pimd_lammps.py --control none|legacy|smooth` selects one and writes `control_used.txt` and a `CONTROL_USED` line; `tools/discontinuity/{run_nve,time_ctrl,frames_compare}.py` and `tools/run_all_checks.py` do the same. **Not yet done:** `inputs/tunneling_ag/run_barrier_lammps.py` (used by the running tunneling batches, so not touched) still writes `NULL` (= legacy) into its LAMMPS input and does not log it; add the option after the batches end.

## 1. Smooth vs legacy on the delivered Ag/El frames [measured]
150 frames (not 130: Ag3 20, Ag4 20, Ag5_2 20, Ag6_2 20, Ag7_1 50, Ag7_2 20), charges ±1, taper 5.7, adapter, same positions.
| case | max \|ΔE\| (kcal/mol) | rms ΔE | max \|ΔF\| per atom | rms ΔF component |
|---|---|---|---|---|
| Ag3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Ag4 | 0.0034 | 0.0012 | 0.117 | 0.012 |
| Ag5_2 | 0.0055 | 0.0023 | 0.111 | 0.011 |
| Ag6_2 | 0.234 | 0.061 | 1.39 | 0.080 |
| Ag7_1 | 0.297 | 0.159 | 1.62 | 0.141 |
| Ag7_2 | 0.442 | 0.181 | 1.61 | 0.171 |
| **all 150** | **0.442** | **0.115** | 1.62 (max component) | **0.118** kcal/mol/Å |
Per term (all 150): bond max 0.074, rms 0.0149; over/under max 0.057, rms 0.0115; **torsion max 0.442, rms 0.117**. Every other term is identical. So the change is almost entirely torsion from bonds with BO between 1e-9 and 0.001 (and the bond-list tail), and it grows with cluster size (Ag3: none; Ag7: 0.1–0.4 kcal/mol). It is a change of the energy surface, not a numerical refinement: relative energies of Ag7 isomers can move by a few tenths of a kcal/mol.

## 2. Does `ereac` accept the same thresholds? [measured]
Yes. In the legacy source the control-file keys `vlbora` (bond-order distance range, default 5.00, reac.f:8444/8614) and `cutof2` (BO threshold for angles/torsions, default 0.001, reac.f:8389/8540) are the same two quantities (`thb_cutoff_sq` has no `ereac` counterpart; it had no effect in the adapter either). No patch was needed; the existing rebuilt `ereac` was used with `    5.2 vlbora` and `1.0E-09 cutof2` added to its control file. Caveat: `ereac` applies `vlbora` when it builds its Verlet pair list (every `nvlist` steps, default 25), so in `ereac` MD the cutoff acts on stale distances; irrelevant for single points.
Validation (30 frames, every 5th, adapter vs `ereac` with the same thresholds, adapter evaluated on the unquantised frame):
| mode | max \|ΔE\| | max \|ΔF\| |
|---|---|---|
| legacy | 6.47e-4 kcal/mol | 1.11e-2 kcal/mol/Å |
| smooth | 6.53e-4 kcal/mol | 1.11e-2 kcal/mol/Å |
These numbers are the 1e-5 Å quantisation of the geo file (the earlier exact-position harness gives 1e-8). Smooth is reproduced at exactly the legacy level, so smooth mode is validated against `ereac` to the quantisation floor, but on frames only; a finer smooth-vs-ereac comparison on identically quantised positions was not run.

## 3. Cost [measured, shared machine, 1 core]
| case | legacy | smooth |
|---|---|---|
| 225 atoms C/El, non-periodic, taper 10 | 3.47 ms / force call | 3.10 ms |
| 120 atoms C/El, periodic (graphene), taper 10 | 11.2 ms | 11.9 ms |
No measurable cost (±10 % is noise at this load). **But on carbon the smooth surface is not usable as is:** at fixed geometry the energy changes by −83.4 kcal/mol (225 atoms) and −1.05 kcal/mol (120 atoms), and a graphene NVE run with `thb_cutoff` 1e-9 was unstable earlier (energy span 38 kcal/mol in 100 fs). Smooth has been validated for Ag/El only; do not use it for C/El without a separate check.

## 4. NVE with mobile Ag under smooth: span versus dt [measured]
Ag7_1 frame 0, 300 K, total-energy span per bead over the run (kcal/mol).
| setup | dt 0.25 fs | dt 0.125 | dt 0.0625 | ratios |
|---|---|---|---|---|
| P = 1 smooth, 1 ps (LAMMPS NVE) | 0.0235 | 0.00596 | 0.00149 | 3.95, 4.00 |
| P = 1 legacy, 1 ps | 0.082 | 0.072 | 0.070 | flat |
| P = 8 smooth, 0.5 ps (i-PI, Ag classical, 8 beads) | 0.0370 | 0.0097 | 0.0024 | 3.8, 4.0 |
| P = 32 smooth, 1 ps (i-PI, thermal start, no equilibration) | 0.0389 | 0.00939 | 0.00282 | 4.14, 3.3 |
Smooth falls as dt² at P = 1, 8 and 32; legacy does not. At P = 8 and 32 the legacy run at dt 0.25 has the same span as smooth in these windows (0.0370 and 0.0389): no threshold crossing occurred in 0.5–1 ps, so P = 8/32 over these short windows cannot discriminate legacy from smooth (P = 1 shows 29 jumps per ps). The longer legacy-versus-smooth P = 32 comparison (4 ps) is below.
The P = 32 runs started from the delivered frame with thermal velocities, not from an equilibrated state; an equilibration run at P = 32 (NVT PILE-L) was stopped because the machine load (load average 25–37 on 18 cores) made it impractically slow.

**Not done:** legacy-versus-smooth over 4 ps at P = 32 (and P = 8). I started both and stopped them: the machine was running at load average 25–37 on 18 cores (the tunneling batches use ~6 cores at 120 %), which made P = 32 about 70× slower than its unloaded 100 steps/s (2200 steps in 10 min). The legacy-versus-smooth separation is therefore shown at P = 1 only (29 jumps per ps there); at P = 8 and 32 the claim is limited to "smooth falls as dt²".

## Recommendation
Adopt `smooth` for Ag/El runs with mobile Ag (silicon-type work), keep `legacy` for all validation; do not use `smooth` for C/El. Report the control file in every run (done in the scripts above). The energy surface differs by up to 0.44 kcal/mol at fixed geometry for Ag7, so smooth results should not be compared with legacy numbers at the 0.1 kcal/mol level.
