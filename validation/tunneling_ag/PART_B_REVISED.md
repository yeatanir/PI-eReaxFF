# Part B (revised, owner decision "Option 2") - intrinsic barrier, gap 4.0 A

Labels: **[measured]** exact rebuilt ereac single points / short i-PI runs done here; **[derived]** arithmetic on measured numbers; **[estimate]** model or guess, not measured; **[to do]**.
Data: `data/landscape_4.0.json` (all measurements below), `gap4p0_landscape.png`; scripts `inputs/tunneling_ag/{landscape_gap,reduce_landscape_gap,analyze_rings,run_ipi_barrier}.py`.

## 1. The landscape at gap 4.0 A  [measured]
System: two Ag3 triangles facing apex to apex, apex-apex distance 4.0 A (apexes at z = -2.0 and +2.0), 6 Ag+, 5 El (left unit 3, right unit 2; the transferring electron is the 3rd electron of the left unit), fixed charges, Ag fixed, other four electrons on their cations (+0.001 A). z = donor-acceptor axis, y in the Ag3 plane, x out of plane.

| quantity | value |
|---|---|
| pit minima (left / right) | z = -1.975 / +1.975 A, energies equal to 0.0000 kcal/mol (mirror symmetry) |
| **middle region vs pit** | **E(z=0) - E(pit) = +0.14 kcal/mol (0.2 kT)**: the middle is a broad minimum, practically degenerate with the pits, *not* clearly lower or higher (lowest point of the middle region is at z = 0). With the other four electrons relaxed: +0.38 kcal/mol. |
| **intrinsic barrier top** | **z = -1.700 / +1.700 A (0.275 A from each pit minimum, 0.30 A from the apex Ag), height 10.51 kcal/mol above the pit = 17.6 kT** (10.37 above the middle) |
| profile E - E(pit) (kcal/mol) at z = -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.3, -1.0, -0.5, 0 | 0.06, 3.24, 8.89, **10.51**, 9.91, 9.00, 7.37, 5.12, 1.64, 0.14 |
| shape | strongly asymmetric: half-height points at z = -1.869 and -1.018 -> only 0.17 A wide on the pit side, 0.68 A on the middle side (FWHM 0.85 A) |
| curvature at the top (quadratic fit, |u| <= 0.08 A) | k_z = **-9.44 eV/A^2**; k_x = +3.02, k_y = +2.95 eV/A^2 (a true 3-D saddle) |
| pit curvature (gap 4.0 fine scan) | k = 44.3 / 45.1 / 45.0 eV/A^2 along z / x / y |
| middle curvature | k_z = 0.56, k_x = 1.14, k_y = 1.12 eV/A^2 |

### Imaginary barrier frequency and crossover temperature  [derived from k_z = -9.44 eV/A^2]
omega_b = sqrt(|k|/m), T_c = hbar omega_b / (2 pi k_B), 1 amu = 103.6427 eV fs^2/A^2.
| m (amu) | omega_b (fs^-1) | (cm^-1) | **T_c (K)** | T_c / 300 K |
|---|---|---|---|---|
| 4.0 | 0.151 | 801 | **183** | 0.61 |
| 1.0 | 0.302 | 1602 | **367** | 1.22 |
| 0.5 | 0.427 | 2266 | **519** | 1.73 |
| 0.1 | 0.954 | 5066 | **1160** | 3.87 |
Reading: at 300 K the 4 amu electron is above its crossover (thermally activated crossing, ring is a nearly compact object at the top); 1, 0.5 and 0.1 amu are below it (deep tunneling regime, rings become unstable/stretched at the top; the instanton picture applies). **Caveat [estimate]:** the barrier is far from parabolic (0.17 A / 0.68 A asymmetric half-widths, 10.5 kcal/mol high) and the curvature was fitted in a +-0.08 A window, so T_c is a local-curvature number; an instanton calculation would give the true crossover.

### Is the 1-D axis scan with everything else fixed a fair picture?  [measured]
- **Ag atoms:** held exactly fixed in the dynamics (`<fixatoms>`): in the 25 fs NVE test their displacement from the start is 0.0 (max over beads and frames).
- **Other four electrons:** in that 25 fs P = 4 NVE run each stays within 0.07-0.16 A of its own cation (so they stay paired); relaxing them with the transferring electron held at the pit / top / middle (L-BFGS on the exact energy and gradient, 12 degrees of freedom) lowers every state by 0.87-1.17 kcal/mol with maximum motion 0.014-0.023 A. The **barrier height (top - pit) changes from 10.45 (fixed) to 10.54 kcal/mol (relaxed), middle - pit from 0.08 to 0.38**: the fixed-others profile is a fair picture to ~0.1-0.3 kcal/mol.
- **Off-axis:** (x, y) grids of +-0.6 A (step 0.15 A) at z = -2.0, -1.9, -1.8, -1.7, -1.5, -1.3, -1.0, 0: at every z the grid minimum is the on-axis point (gain 0.00 kcal/mol): the electron has no lower-energy path around the barrier within this window, and the transverse curvature at the top is positive. The grid is coarse (0.15 A vs a 0.2 A wide pit) - a finer grid near the top is a cheap follow-up [to do].
- **Not captured [important]:** if all five electrons are ring polymers of mass m, the four "other" electrons also delocalise (free-ring radius 0.64 A at 0.1 amu > the pit-to-top distance 0.275 A), so the fixed-others picture breaks down. **Default for Part B: one quantum electron** (`run_ipi_barrier.py --one-quantum-electron`: transferring electron = ring of mass m, the other four El are classical via `<classical_atoms>` at their force-field mass 1 amu, Ag fixed) [smoke-tested: starts and runs, 6 steps]. The five-ring variant is an optional extra.
- f10.5 position rounding of the ereac geo path (1e-5 A) and the 1e-3 gtol of the relaxation limit these numbers at the ~0.01 kcal/mol level.

## 2. System, masses, beads and time steps
Pit stiffness (k = 45 eV/A^2) sets P and dt if the ring visits the pit; the middle region alone would need far fewer beads (A6 table). Rule: (beta hbar omega / P)^2 / 24 < 1 %, omega dt < 0.2, with the dt rounded down. **[derived]**
| mass (amu) | label | free-ring radius (3-D rms, 300 K) | omega_pit | beta hbar omega | P (2P) | dt (fs) | T_c (K) |
|---|---|---|---|---|---|---|---|
| 4.0 | near-classical *barrier crossing* (the pit is still quantum: beta hbar omega = 8.4, zero-point energy ~ 4 kT) | 0.10 A | 0.329 | 8.4 | 18 (36) | 0.5 | 183 |
| 1.0 | force-field mass | 0.20 A | 0.659 | 16.8 | 36 (72) | 0.25 | 367 |
| 0.5 | | 0.28 A | 0.932 | 23.7 | 50 (100) | 0.20 | 519 |
| 0.1 | **zero-point delocalisation regime** (ring radius 0.64 A larger than pit width and than pit-to-top distance: label results "delocalised/zero-point", not "tunneling through the barrier") | 0.64 A | 2.084 | 53.1 | 110 (220) | 0.09 | 1160 |
Classical baseline: P = 1, 1 amu (static distributions do not depend on mass), dt 0.25 fs. [P is rounded up from 17.1/35/49/109; whether the pit numbers are really needed is decided by the P vs 2P histogram test (A6) - **[to do]**.]

## 3. The Gaussian boost on the intrinsic barrier
Placement: centred on the barrier tops, **at s = -1.70 A and s = +1.70 A from the midpoint (z = +-1.70 in the case frame; both sides)**, along the axis n = z. The current `gaussian_barrier` driver has one centre: use **two barrier components (two ffsockets)** or add a second-centre option [to do, small]. Boosts V0 = 0 (intrinsic only), +5, +10 kcal/mol, FWHM 0.30 A (sigma 0.127 A; default) - or 0.40 A.
Effect of the Gaussian on the profile (merged measured profile + Gaussian; top stays at -1.70 A in all cases) **[derived]**:
| FWHM | boost | shift of the pit energy | barrier above the shifted pit | in kT |
|---|---|---|---|---|
| 0.3 | 0 / +5 / +10 | 0 / +0.37 / +0.68 | 10.51 / 15.14 / 19.82 | 17.6 / 25.4 / 33.2 |
| 0.4 | 0 / +5 / +10 | 0 / +1.11 / +2.16 | 10.51 / 14.40 / 18.35 | 17.6 / 24.1 / 30.8 |
Note the net barrier rises by less than the boost because the Gaussian also lifts the pit (it is only 0.275 A from the top); report the net barrier, not the nominal boost. FWHM 0.3 A is the default because it disturbs the pit least; the intrinsic barrier is itself wider on the middle side, so the Gaussian mainly adds a spike at the top.

## 4. Run table
Electron: one quantum electron (default). Production by PILE-G after PILE-L equilibration (i-PI fork `classical-atoms`). 3 seeds per row. **Every PI row is run from two starts (sampling check, section 5).**
| rows | masses | boost (kcal/mol) | beads | seeds | starts | runs |
|---|---|---|---|---|---|---|
| classical | 1 (P = 1) | 0, +5, +10 | 1 | 3 | pit, middle | 18 |
| PI | 4, 1, 0.5, 0.1 | 0, +5, +10 | P and 2P (table above) | 3 | pit, middle | 4 x 3 x 2 x 3 x 2 = 144 |
| control | boost 0 is the barrier-off-Gaussian control (intrinsic landscape only); extra control: no Gaussian component present (identical to boost 0, already shown byte-identical for V0 = 0 in A4) | | | | | included |
**Staging (recommended) given the cost below:** stage 1 = boost 0, P only, 1 seed, both starts, 4 masses + classical (10 runs); stage 2 = add 2P at 1 and 0.5 amu; stage 3 = boosts +5, +10 at 1 amu and 0.5 amu; the 0.1 amu rows last.
Run length (proposal [estimate]): 2 ps PILE-L equilibration + 10 ps PILE-G production per run; extend until the block errors of the zone fractions are < 20 % of their value. Classical P = 1 crossing estimate: rate ~ (omega_pit / 2 pi) exp(-17.6) ~ 2e6 s^-1, i.e. one crossing per ~0.4 microsecond - the classical rows will show a *statistically clear zero* crossings, and two classical starts will give two non-communicating populations (this is the sampling check working, not a failure).

## 5. SAMPLING CHECK (mandatory)
Start every PI run twice: **start A** = all beads collapsed in the left pit (z = -1.975), **start B** = all beads collapsed in the middle (z = 0); other atoms identical. After the same equilibration compare zone fractions (donor side z < -1.70 | middle | far unit z > +1.70) between A and B in the production phase. Acceptance: all three fractions equal within 3 combined block errors **and** both starts show ring exchange between zones (centroid zone changes > 0 in both). If they differ: report as **leakage out of the pit / non-ergodic sampling**, state the dwell behaviour, **do not claim equilibrium**, and report each start separately (the quantity then is a *time-limited residence statistic*, not an equilibrium density). Also compare the far-unit fraction to the left-right mirror (z -> -z) expectation for start B (should be symmetric for a symmetric landscape; deviation > 3 sigma flags insufficient sampling).

## 6. Cost  (**[estimate]; per-force-call time is a GUESS - please supply the measured time per bead call for LAMMPS + Metal-Reax**)
Cost = (1000 / dt) x P x t_call / n_parallel seconds per simulated ps; n_parallel = 2 (CPU cap). Guess t_call = 10 / 20 / 30 ms. (An 11-atom system may be much faster than 10-30 ms on CPU-64 - not measured; the guess may be pessimistic.)
| mass | dt | P | wall per ps (10 / 20 / 30 ms) | 2P |
|---|---|---|---|---|
| 4 | 0.5 | 18 | 180 / 360 / 540 s | 360 / 720 / 1080 s |
| 1 | 0.25 | 36 | 720 / 1440 / 2160 s | 1440 / 2880 / 4320 s |
| 0.5 | 0.20 | 50 | 1250 / 2500 / 3750 s | 2500 / 5000 / 7500 s |
| 0.1 | 0.09 | 110 | 6111 / 12222 / 18333 s | 12222 / 24444 / 36667 s |
Full table (144 PI runs of 12 ps) at 20 ms: ~ 3000 h wall on 2 cores - **not feasible as written**; stage 1 alone (P only, 1 seed, both starts, 12 ps) is ~110 h at 20 ms and 0.1 amu is 87 % of it. Levers [to do/ask]: measured t_call, more cores, shorter 0.1 amu runs, ring-polymer contraction or a cheap reference potential with the eReaxFF correction (MTS), and checking whether P = 36/50/110 (pit rule) is really needed (the histogram convergence test may allow fewer).

## 7. What each figure shows
1. **Landscape (F1)**: on-axis profile with pits, middle, barrier tops (z = +-1.70), the +5/+10 Gaussian-boosted profiles, and the (x, y) plane at the top; table of T_c per mass. Purpose: the energy scales everything else is judged against.
2. **Bead density along the axis (F2)**, one panel per mass, classical vs PI, boosts 0/+5/+10 overlaid, zone boundaries shaded: shows density beyond the barrier top and in the middle.
3. **Zone fractions vs mass (F3)**: fraction of bead density in the middle region and on the far unit (and donor side), for each boost, with P and 2P overlaid and block errors: the monotonic mass dependence and its fall with boost (success criterion), and bead-number convergence.
4. **Ring snapshots (F4)**: a few configurations in the z-y plane (Ag atoms, other electrons, barrier tops) with beads on both sides of a top: direct evidence of straddling rings; plus the straddling-frame fraction vs mass (`analyze_rings.py --top 1.7`).
5. **Sampling check (F5)**: zone fractions for start A vs start B per mass: flags leakage / non-equilibrium.
6. **Penetration vs boost and vs mass (F6)**: the success curve (monotonic fall with boost at fixed mass; monotonic rise as mass decreases at fixed boost), with T_c marked on the mass axis (crossover from activated to deep-tunneling regime).

## 8. Software state for this protocol
- Done: barrier driver + tests (A4), ring analysis incl. zones/straddling (`--top`, 4 tests), single-quantum-electron launcher option, all landscape measurements. 
- To do: second barrier centre (or two components); production engine (LAMMPS + Metal-Reax multi-bead - workstream G); A3 conserved-quantity runs for the 4 masses; A6 P vs 2P histogram test; finer off-axis grid at the top.
- Caveat for write-up wording: electrons are distinguishable pseudo-particles of a chosen mass; the output is an equilibrium (or time-limited) density, not a current.
