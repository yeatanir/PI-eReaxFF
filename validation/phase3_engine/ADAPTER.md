# LAMMPS adapter for eReaxFF (task 3.5) — worktree `forks/metal-reax-lammps`, branch `ereaxff-lammps`

Setup: pinned LAMMPS stable_30Sep2026 @ 8de817dd fetched/verified with `tools/fetch_lammps.sh --full` and built with `tools/build_lammps_reference.sh` into `external/{lammps,build_lammps,install_lammps}` (outside any repo). Plugin built CPU-64 (Metal OFF). 4 build/test threads, OMP_NUM_THREADS=1.
Practical: the build tree path must not contain a space (`eReaxFF PIMD`): the plugin tests failed from that path (14/42), so I used a space-free symlink to the same tree (`scratchpad/eb`).

## 1. Baseline (before any change)
ctest 42 tests: 38 pass, 4 fail. All 4 were a build-system bug, not physics: `$<$<BOOL:METAL>:--metal>` expands to an EMPTY argument when Metal is OFF, so the Python harnesses (decision_census, threshold_scan, fp32_envelope, lammps_eem) got "unrecognized arguments". Fixed (variable `RM_METAL_ARG`); then 42/42 pass, including all in-LAMMPS tests (INT-2, charge models, hybrid, analysis, ...).

## 2. Adapter (commit d57a1cd and later)
`pair_style reaxff/metal NULL ereaxff yes [electron_taper <A>]`
- eReaxFF terms (electron pair, sharing, lone pair, over/under) run inside `compute()`; legacy single-precision constants (c1c, 4/3 exponent) used.
- Charges come from the data file (legacy `icharg 8`); no charge fix allowed (error if present); no polarization energy (legacy total excludes `ech`, reac.f:5383).
- **15th pvector slot `eelect`** (EnergyTerm::Electron); `compute pair reaxff/metal` now returns 15 values. Contract test updated deliberately (14 stock slots + 1). eng_vdwl includes the electron energy.
- Refused with explicit errors: backend metal, lgvdw, per-atom energy/virial, Ho, every element outside the allow-list, electron taper <= 0 or > nonbonded cutoff.
- **Element allow-list** `kEreaxffElements = {C, EL, AV}` in `plugin/adapter/pair_reaxff_metal.cpp` (comment says how to extend). Ag is deliberately NOT accepted (see §5).
- Ghost shell: required width = max(nonb, hbond, 2*bond_cut) + sqrt(27.63/vincr_min) (sharing range where exp(-vincr r^2) < 1e-12); checked at run time with a clear error (`shellcheck no` downgrades to a warning as before). Note: checked in compute() path (host_atom_set), as for the stock plugin, not in init_style. Tests set `comm_modify cutoff 11.5`.
- Negative test: mapping a type to H gives "element 'H' is not in the eReaxFF allow-list".

## 3. Parsed ffield == legacy fixed-column reader (test `ffield_legacy_columns_*`)
`tests/python/legacy_ffield_slice.py` re-implements the legacy formats independently ((f10.4), (1x,a2,10f9.4), (3x,10f9.4), (2i3,8f9.4), (6x,8f9.4)) incl. the post-read rule vval3 := valf for mass < 21 (reac.f:4622-4629); `test_ffield_legacy_columns.cpp` compares every general (39), atom (32 per atom) and bond (16 per bond) value with Metal-Reax's parse. 8 force fields registered (Jamil 225 np and vpar38=0, p60 np/periodic, phase4_frames, sharing A/E/J): all pass (0 mismatches). Negative control: widening `181.0000` to `1810.0000` is caught (exit 1); a benign variant passes. The first run found one real difference class, the vval3 post-read rule above (now in the reference slicer).

## 4. In-LAMMPS validation (`tests/lammps/run_ereaxff.py`; data file from the geo, charges from MOLCHARGE, units real, newton on)
Energy = LAMMPS pe vs reference ESTRC; terms = pvector slots vs the 14 reference terms (relative to max(1,|ref|)); force = -gradient, max |dF| relative to the largest reference force.
| Case | |dE_total| kcal/mol | max term rel (no torsion) | torsion rel | force rel |
|---|---|---|---|---|
| Jamil 225 atoms, non-periodic, taper 10 | 2.0e-7 | 3.2e-10 | 1.0e-10 | 7.5e-10 |
| 60 C/60 El, non-periodic | 2.4e-8 | 2.5e-10 | 1.0e-10 | 6.0e-11 |
| 60 C/60 El, periodic | 9.6e-6 | 1.1e-9 | 6.5e-7 | 9.1e-8 |
| Jamil 225 atoms, periodic | 5.8e-6 | 2.7e-9 | 6.0e-8 | 1.8e-7 |
| **51 ereac MD frames** (phase4_frames, 60 C/60 El periodic, steps 0..3000, electrons 0.11 A from carbons) | worst 1.3e-5 | worst 7.0e-8 (lone_pair, late frames); 1.6e-9 for frames up to step 2000 | worst 6.6e-7 (early frames; falls to 1e-9..1e-8 later) | worst 1.1e-6 (1e-4 kcal/mol/A absolute on forces up to 920) |
Per-frame lines: `run_ereaxff.py ... --frames frames.json --periodic` (gates used: terms 1e-7, torsion 1e-6, force 2e-6; all 51 pass). Torsion carries the known stock-vs-legacy residual (torsion fork); it is largest in the early frames.
Not covered by the frames: two electrons near one carbon, an electron midway between carbons (only the small cases A-J in the standalone harness; not yet run through LAMMPS).

## 5. Short MD through LAMMPS integrators (`tests/lammps/run_ereaxff_md.py`, fix nve, periodic 60 C/60 El, Maxwell 298 K, masses C 12.0, El 0.5)
Forces on El drive fix nve correctly: total energy is constant to 0.05 kcal/mol (dt 0.1) / 0.003 (dt 0.025) for the first ~110 fs, then the system collapses (T -> thousands of K, etotal falls by ~90 kcal/mol) at ~110-120 fs, independent of dt (0.1, 0.05, 0.025 fs give the same collapse time and size). The same buckling/collapse onset is in the ereac reference run (phase4_frames README: "collapse from 200 fs", validation/phase0_reference conclusion 5), so it is the force field/electron model, not the adapter; the script's 0.05 kcal/mol NVE gate therefore FAILS on a 200 fs run (reported honestly: run shorter, <= 100 fs, or accept). fix temp/berendsen not yet run.

## 6. Av in the legacy tables (for the Ag/El target; published ffield, parameter 38 = 5.7)
- elecpos energy side (reac.f:8844): `Av` -> avalc += eshare, stlpc += eshare (same branch as Li, Lv, P, B, Al). Metal-Reax `share_sign_of("AV") = +1`: matches.
- znucel (reac.f:4614): Av = 3.0 (also Ag 4612 and Lv 4613 = 3.0). Metal-Reax core-charge table has AV = 3: matches.
- Derivative side (poten.f sgn tables): Av is in the sgn=1 (i.e. davalc sign +) list at poten.f:493, 517 (list differs: also 'Ag'), 561/589, 1152 and the later copies; one copy at poten.f:517 puts 'Av' in the sgn=2 list too (line 517 block is the hole branch: C,H,Li,B,Lv,P,Av,Ag,Al -> sgn=2 vs O,F,N,Os -> sgn=1), i.e. the opposite sign convention in the hole branch. Hole branches are irrelevant (Ho refused).
- Quirks to keep in mind: (a) `sgn` is not reset per atom, so an element in neither list inherits the previous atom's sign (the audit's R10 stale-variable); (b) `'Ag'` has derivative-table entries but NO elecpos energy branch, so Ag would get no valence shift in the energy but a derivative term: inconsistent, hence Ag is not on the allow-list; (c) poten.f:1152 tests `qa(i1).eq.'Os'` where i7 is meant (typo), which only matters for Os; (d) some copies write 'P' without the trailing blank (character*2 comparison pads, so harmless).
- Not verified by me: that the full Metal-Reax over/under + lone-pair derivatives for an Av/El system reproduce `ereac` (no Av reference data exists yet). Needed: an Av/El export from the rebuilt ereac with the published ffield (taper 5.7) and the same 51-frame protocol; the ffield columns are fixed-format and the legacy post-read rule vval3:=valf applies to mass < 21 only (Av mass is above that, so vval3 is read as written).

## Not done
- fix temp/berendsen run; two-electron / midway-electron cases through LAMMPS; ghost-shell check in init_style (currently at compute); multi-rank; Av/El validation (no reference yet); upstream of nothing.

---
# Addendum — Ag/El target (published ffield, taper 5.7), MISC/fix ipi

## A1. All Ag/El reference frames through LAMMPS (Av on the allow-list; `tests/lammps/run_ereaxff.py --taper 5.7`; masses read from the ffield: Av 107.8682, El 1.0)
**Reference defect found first:** in every delivered `phase4_frames_ag/*/geo` the line `END` comes BEFORE the `MOLCHARGE` lines (the Jamil geos have them before `END`). The legacy reader stops at `END`, so `nmcharge = 0` and all charges are 0: the reference "coulomb" term is 0.00 in every frame and `ereac` is effectively running neutral, uncharged Ag/El. (Reproduced by running the rebuilt export `ereac` on Ag3 frame 0: Coulomb 0.) The references are exact for q = 0; the intended Ag+(+1)/El(-1) system was never exercised. So two comparisons:

| Case (frames) | charges | max \|dE_total\| kcal/mol | worst term rel (no torsion) | worst torsion rel | worst force rel |
|---|---|---|---|---|---|
| Ag3_md (20) | q=0 (as delivered) | 1.7e-10 | 7.2e-11 | - | 1.2e-10 |
| Ag4_md (20) | q=0 | 2.9e-10 | 5.3e-11 | 4.3e-11 | 6.8e-11 |
| Ag5_2_md (20) | q=0 | 4.0e-10 | 2.0e-11 | 5.8e-11 | 6.5e-11 |
| Ag6_2_md (20) | q=0 | 7.3e-9 | 2.1e-11 | 4.2e-10 | 6.0e-10 |
| Ag7_2_md (20) | q=0 | 1.3e-8 | 2.2e-11 | 3.2e-10 | 5.9e-10 |
| Ag7_1_md (50) | q=0 | 8.7e-9 | 1.2e-11 | 2.1e-10 | 4.0e-10 |
| Ag3_md frames 0,5,10,19 | **+1/-1, new exact refs** | 6.2e-11 | 5.0e-11 | - | 8.7e-12 |
| Ag7_1_md frames 0,5,10,19 | **+1/-1, new exact refs** | 8.7e-9 | 1.2e-11 | 2.0e-10 | 1.3e-10 |
The charged references (`frames_charged.json` in Ag3_md and Ag7_1_md, made by `project/tools/make_ag_charged_refs.py`: MOLCHARGE moved before END, one ereac step per frame, positions quantised to the geo format) have Coulomb = -365.58 (Ag3) and -872.91 (Ag7_1) kcal/mol; LAMMPS matches them, so Av/El Coulomb is exercised and agrees. Recommendation: regenerate all Ag frame sets with MOLCHARGE before END before using them for anything charged.
Av handling: **nothing needed beyond the allow-list.** No compat option was required: Av sharing sign (+), Z=3, the heavy-atom switch (the legacy `amas>31` vs Metal-Reax `mass>21` agree for Av, mass 107.9) and the single Av-Av-Av angle (theta0 = 179.99) plus its torsion record are all reproduced at <= 7e-11 relative in the terms and <= 6e-10 in forces over 150 frames (sin(theta) is small only near collinear Ag3 configurations; none of the frames shows an anomaly). The legacy sgn-table quirks (ADAPTER.md §6) did not show up because the 14 term and force comparisons are at round-off level.
Caveat: the clusters stay bound with El 0.01-0.05 A from their cation (sharing ~1): this is the strong-sharing regime, and it passes.

## A2. fix nve through LAMMPS, Ag7_1 (`tests/lammps/run_ereaxff_md.py --nonperiodic --taper 5.7 --dt 0.25 --steps 20000 --temp 300`; velocities created at 300 K with LAMMPS, seed 4928459; masses Av 107.8682, El 1.0; non-periodic box)
| run | etotal start -> end (kcal/mol) | max-min over 5 ps (41 samples, every 125 fs) | mean T |
|---|---|---|---|
| intended charges +1/-1 | -1163.0116 -> -1162.9797 | 0.064 | 205.7 K |
| q = 0 (as the ereac reference) | -290.0967 -> -290.0687 | 0.076 | 207.4 K |
The ereac reference run spans 0.14 kcal/mol over 5 ps (its velocities/seed differ and it is sampled differently), so LAMMPS conserves energy at least as well. Temperature: starts at 300 K and averages about 206 K (virial/equipartition with a 14-atom cluster started in the minimum-energy geometry: half the kinetic energy goes to potential). The ereac temperature trace was not supplied, so no point-by-point temperature comparison was made. Berendsen was not run.

## A3. fix ipi (MISC) and the i-PI client protocol
- Stock build had PLUGIN QEQ REAXFF only; I enabled `-DPKG_MISC=on` in the same build tree (incremental rebuild + install, `lmp -h` now lists MISC). Plugin rebuilt unchanged and loads; ctest 64/64 pass afterwards (64 incl. the reference-data and ffield tests).
- `project/inputs/ag_pimd/run_ipi_lammps.py`: writes the i-PI server input (template `inputs/jamil_pimd/run_ipi_jax.py`; masses and Gaussian 300 K velocities in m/s, unix socket, nve, dt 0.25 fs, optional `<classical_atoms>`), starts i-PI (fork `classical-atoms`) and P LAMMPS clients (`fix ipi <sock> 32123 unix`), and runs the same system with `fix nve` from the identical positions and velocities. LAMMPS cannot take paths with spaces, so the script copies the ffield next to the input.
- **P = 1, Ag7_1, intended charges, 10 and 100 steps:** max |x_iPI - x_LAMMPS-nve| = 5.8e-7 A and 5.9e-7 A (largest displacement 0.06 A; the residual is the 9-digit RESTART precision plus fix ipi's unit constants 0.52917721 and 3.1668152e-06/boltz). Also with `--classical-atoms`: identical (5.9e-7 A). So the unit conversion mismatch risk R9 is below 1e-6 A here.
- Gotcha: i-PI's xyz position output has only 6 significant digits (`4.00015e+01`), which looks like a 5e-5 A discrepancy at 40 A coordinates; read positions from RESTART instead (done in the script).
- i-PI conserved quantity over the 100 P=1 steps sampled every step: span 1.3e-3 eV (0.03 kcal/mol); not compared with a per-step LAMMPS trace (LAMMPS energies above were sampled every 125 fs), so I do not claim equivalence of the instantaneous fluctuation.
- Not done: P > 1 (several clients) was not run; `matching='lock'` and shm/batch transports untested; fix ipi on a periodic box not tested; no Ag/El comparison against the JAX driver.
