# Tunneling demo - Stage 0 and Stage 1 (2026-11-11)

Labels: **[measured]** run here, number read from output; **[derived]** arithmetic on measured numbers; **[estimate]** model or guess; **[not done]**.
System: two Ag3 triangles apex to apex, gap 4.0 A, 6 Ag+ (fixed), 5 El; transferring electron = atom 11 (the left apex electron), the other four El stay on their cations. Engine: LAMMPS + Metal-Reax (`pair_style reaxff/metal ereaxff yes electron_taper 5.7`, CPU-64, `fix ipi`) through i-PI (fork `classical-atoms`), PILE-L, 300 K, tau 100 fs. All numbers are for the machine as it was loaded (another user job holding 4 cores, load average 7-9 on 18 cores); my own use was at most 4 processes.
Scripts: `inputs/tunneling_ag/{run_barrier_lammps,stage0_nve,analyze_stage1,stage1_report,schrod1d,schrod3d}.py`. Raw: `data/schrod1d.json`, `data/schrod3d_h0.04.json`, `data/schrod3d_h0.03.json`, scratch runs in `scratch/tun_stage/`.

## Headline results

1. **The cost estimate was off by a factor of about 50 to 1000.** One force call of LAMMPS + Metal-Reax for this 11-atom system takes **10.4 us** [measured]; i-PI needs **0.21 ms per bead request** [measured, P = 36-72], so i-PI (Python), not the force, is the bottleneck. The earlier guess of 10-30 ms per call was wrong.
2. **The "3 x redundant evaluation" is not in `fix ipi`.** It is in the i-PI fork (`project_classical` rewrote the internal modes of the classical atoms with zeros at every thermostat step, which invalidates the bead positions and forces a recomputation). One-line fix, trajectories byte-identical, 3 -> 1 force evaluation per bead and step with NVT; 1.9 x faster at P = 32. [measured]
3. **The physical picture is not "electron in the donor pit, tunnelling through a barrier". At 300 K the electron's equilibrium is the broad middle region, at every mass tested, classically and quantum-mechanically.** The 3-D Schroedinger solve of the exact eReaxFF landscape gives the probability of finding the electron in the pit regions (|z| > 1.7 A) as 5e-5 (4 amu), 1e-8 (1 amu), 5e-8 (0.5 amu), 4e-5 (0.1 amu), against 7e-3 classically. The pit is a high-energy state in the full 3-D problem (3-D pit ground state 5.6 / 9.7 / 12.1 / 18.7 kcal/mol above the global ground state at 4 / 1 / 0.5 / 0.1 amu). The PIMD agrees: started in the pit, the ring leaves it within about 1 ps (1 amu) or 0.5 ps (0.5 amu) and never returns. A classical (P = 1) electron stays in the pit for the whole 10 ps. **This changes what "penetration" can mean in this system; see "What this means for the demonstration" below. It needs your decision before Stage 2.**
4. The slow conserved-quantity wander of the Ag7_1 test **does return at P = 64 (same size as P = 32, 0.105 vs 0.102 kcal/mol per bead, independent of dt)**; it is not a ring-frequency x time-step effect. In the production setup (Ag fixed, 4 classical El, 1 ring) it is **absent** (slow span 3e-4 kcal/mol per bead).

---

## Stage 0a - cost of a force call, and the redundant evaluation

### Measured times (machine loaded as described above)
| what | time | note |
|---|---|---|
| LAMMPS + Metal-Reax, standalone NVE, 40 000 steps, 11 atoms, 1 core | **Pair 10.4 us/step** (0.415 s / 40 000), loop 13 us/step | 3 repeats within 2 %; the force alone |
| same through `fix ipi` with a minimal protocol server (Python, 20 000 requests, one geometry each) | 42.8 us/request; with a zero pair style 40.8 us | server-side Python and socket dominate |
| LAMMPS Pair time inside a `fix ipi` run (mini server, 20 000 steps) | 12-13 us/request, 1 LAMMPS step per request | so `fix ipi` itself does **1** evaluation per request |
| i-PI + 1 client, P = 8, 1000 steps NVT, after the fix | 0.33 ms/request | |
| i-PI + 1 client, P = 32, 1000 steps NVT, after the fix | 0.24 ms/request | before the fix: 0.46 ms |
| production runs (P = 36, 50, 64, 72; 4 jobs sharing the machine) | **0.21-0.22 ms/request** | e.g. P = 36: 7.6 ms per i-PI step |
So a bead request costs about 0.2 ms, of which about 0.012 ms is the eReaxFF force [measured]. More LAMMPS clients cannot help; the Python server is the limit. [derived]

### Where the factor 3 comes from, and the fix
A1 had seen 3 client force evaluations per i-PI step per bead. Reproduced and localised here [measured]:
| configuration (P = 4, 200 steps, dump at every client step) | client force evaluations | expected |
|---|---|---|
| NVE, classical atoms | 805 | 804 |
| NVT PILE-L, classical atoms (before) | **2405** | 804 |
| NVT PILE-L, Ag + other El **fixed** (no classical atoms) | 805 | 804 |
| NVT PILE-L, classical atoms, after the fix | **805** | 804 |
Cause: `DummyIntegrator.project_classical` (i-PI fork, branch `classical-atoms`) always executed `self.nm.qnm[1:, cd] = 0.0`; assigning to `qnm` marks the bead positions as changed, so i-PI recomputed the forces after every thermostat half-step. Fix: assign only if something is non-zero (3 added lines in `ipi/engine/motion/dynamics.py`). `fix ipi` and the plugin were not changed (the directive pointed at `forks/metal-reax-lammps`; the redundancy was not there).
Verification [measured]: the 200-step run before and after the fix, same seed: `sim.out` identical (max difference 0.0), all 4 bead trajectory files byte-identical, and every force frame of the new run matches (bitwise) a frame of the old run with the identical positions (max force difference 0.0). Timing, 1000 steps NVT: P = 8: 4.2 s -> 2.6 s; P = 32: 14.6 s -> 7.8 s. The i-PI classical-atoms test file (23 tests) passes after the fix (`ca_tests_after_fix.log`: 23 passed in 18 min).
Also confirmed [measured]: the LAMMPS + Metal-Reax profile along z (121 points, transferring electron at (0.001, 0, z)) agrees with the exact rebuilt `ereac` to 6.5e-5 kcal/mol. An isolated singular point exists: the electron exactly on a core (r = 0) gives +22.5 kcal/mol; at r = 1e-8 A the value is smooth (+0.016). It has measure zero in sampling.

### Revised cost per simulated ps (t_req = 0.215 ms measured) [derived]
Seconds of wall time per ps = (1000 / dt) x P x 0.215 ms. The i-PI process is single-threaded, so one run is about one core; the LAMMPS client is about 6 % busy.
| mass (amu) | P | dt (fs), omega_max dt = 0.63 / 0.35 / 0.31 | s per ps | 10 ps run |
|---|---|---|---|---|
| 4 | 18 | 0.44 / 0.25 / 0.22 | 8.8 / 15.5 / 17.6 | 1.5-3 min |
| 1 | 36 | 0.22 / 0.125 / 0.11 | 35 / 62 / 70 | 6-12 min |
| 0.5 | 50 | 0.16 / 0.09 / 0.08 | 67 / 119 / 134 | 11-22 min |
| 0.1 | 110 | 0.073 / 0.04 / 0.036 | 324 / 591 / 648 | 0.9-1.8 h |
(omega_max = 2 P kT / hbar = 0.0786 P fs^-1 at 300 K.) The previous table's "3000 h" becomes about 140 core hours for the whole 144-run table at 12 ps and omega_max dt = 0.35, of which 75 % is the 0.1 amu rows; without the 0.1 amu rows it is about 36 core hours. [derived, and runs of 5 ps would do, see Stage 1] P = 110 timed afterwards: 0.227 ms/request (600 steps, one other job running) [measured].

---

## Stage 0b - Schroedinger prediction

### 1-D solve on the measured axis profile (as requested)
V(z) from a 2 mA scan of the exact landscape (LAMMPS + Metal-Reax, transferring electron at (0.001, 0, z), all else fixed), finite-difference Hamiltonian, 60 lowest states, 300 K Boltzmann sum (`schrod1d.py`). Profile: pit minimum z = -1.986, top 10.57 kcal/mol at z = -1.700 (agrees with the earlier 10.51). The pit is strongly anharmonic as you said: a harmonic well of the measured stiffness (45 eV/A^2) would reach 39.2 kcal/mol at 0.275 A, against 10.57 measured. [measured/derived]
| mass (amu) | levels above ground state (kcal/mol) | 2 pit-localised states at | their weight in the middle (leakage past the top) | thermal zones L / mid / R (z < -1.7 / between / z > 1.7) | classical 1-D Boltzmann |
|---|---|---|---|---|---|
| 4 | 0, 0.55, 1.08, 1.61 | 1.95 | 0.0 | 0.022 / 0.957 / 0.022 | 0.115 / 0.771 / 0.115 |
| 1 | 0, 1.07, 2.11, 3.10 | 3.75 | 0.3 % | 0.0015 / 0.997 / 0.0015 | |
| 0.5 | 0, 1.50, 2.91, 4.24 | 4.96 | 2.5 % | 0.0002 / 0.9996 / 0.0002 | |
| 0.1 | 0, 3.12, 5.71, 7.19 | 8.1 and 26.7 (mixed with middle states) | 38-44 % | 0.0001 / 0.9999 / 0.0001 | |
The splitting of the symmetric/antisymmetric pit pair is below 1e-3 kcal/mol for 4, 1 and 0.5 amu (5e-7 at 1 amu): the two pits are decoupled in the pit-state picture; the time to go from one to the other is not a "splitting" quantity here. [measured]
**Caveat that matters:** this 1-D picture leaves out the transverse motion. It is not what the system does (next section).

### 3-D solve (needed: the transverse zero-point energy is as large as the barrier)
The electron's x and y motion is stiff in the pit (45 eV/A^2, transverse zero-point energy about 10 kcal/mol at 1 amu) and soft in the middle (1.1 eV/A^2, about 1.6 kcal/mol) and at the top (3.0). Adding these to the axis profile lowers the effective barrier from 10.5 to about 3 kcal/mol at 1 amu and puts the middle below the pit. So I solved the full 3-D problem: V(x, y, z) of the same frozen landscape from LAMMPS + Metal-Reax (cell-centred grid, x, y in [0, 1.2] A, z in [-3, 3] A, 2nd-order finite differences, four parity sectors in x and y; V is mirror-symmetric in all three coordinates) (`schrod3d.py`). Thermal density by summing all states up to 14 kT above the ground state. [measured]
| mass (amu) | states summed | ground state E0 (kcal/mol, h = 0.04 / 0.03) | 3-D pit ground state above the global ground state (h = 0.04 / 0.03) | P(\|z\| > 1.7), both pit regions (h = 0.04 / 0.03) | rms z (A) | P(\|z\| > 0.5) | P(\|z\| > 1.0) |
|---|---|---|---|---|---|---|---|
| 4 | 464 | 1.05 / not run | **5.63** kcal/mol = 9.4 kT / not run | 5.0e-5 / not run | 0.232 | 0.027 | 1.3e-4 |
| 1 | 77 | 2.114 / 2.121 | **9.74 / 9.78** = 16.3 kT | 1.0e-8 / 9.1e-9 | 0.253 | 0.045 | 2.2e-4 |
| 0.5 | 34 | 2.989 / 2.996 | **12.07 / 12.17** = 20.4 kT | 5.2e-8 / 4.4e-8 | 0.277 | 0.066 | 5.8e-4 |
| 0.1 | 7 | 6.655 / 6.663 | **18.65 / 18.89** = 31 kT | 3.8e-5 / 3.4e-5 | 0.392 | 0.19 | 1.3e-2 |
| classical 3-D Boltzmann | | | | 7.1e-3 | 0.279 | 0.032 | 7.2e-3 |
Convergence: the two grid spacings agree to 0.5 % in the energies and 20 % in the 1e-8 tail probabilities (the numbers that matter are the pit energies). Wavefunction weight in the outermost grid cells: 1e-4 (1 amu) to 1.8e-3 (0.1 amu, so the 0.1 amu rows carry a box-size error of about 1 %). 4 amu was solved at h = 0.04 only (pit rms ground width 0.05 A against h = 0.04: marginal) [read: treat 4 amu as +-10 %]. The potential is the one with the other four electrons frozen on their cations (relaxing them lowers the pit and the middle by about 1 kcal/mol each and changes their difference by 0.1-0.3 kcal/mol, from the earlier test).
**Consequences [derived]:** (i) The equilibrium electron density is in the middle region for every mass, 99.99 % at 1 amu; the donor pit carries a Boltzmann factor exp(-16) at 1 amu. (ii) The pit is a *metastable* state of the ring-polymer only if the ring is smaller than the pit: the free-ring radius at 300 K is 0.10 A at 4 amu but 0.20 A, 0.28 A, 0.64 A at 1, 0.5, 0.1 amu, against a pit-to-top distance of 0.275 A. (iii) The effect of mass on the *equilibrium* density is a width effect: rms z 0.23 -> 0.25 -> 0.28 -> 0.39 A from 4 to 0.1 amu, P(\|z\| > 1.0) from 1e-4 to 1e-2. It is monotonic, but the density is already entirely inside the barrier tops. (iv) Visible effect: yes, as a width change and as a pit-escape *time*, not as density beyond a barrier that confines it.

### Run lengths this implies [derived]
The ring leaves the pit within about 1 ps at 1 amu and 0.5 ps at 0.5 amu (Stage 1); the middle-region density converges slowly because z-fluctuations of a broad, soft basin decorrelate in ps (block errors in Stage 1). 10 ps gave +-4 % on the rms width at 1 amu; to get +-1 % one needs 100 ps (a P = 36 run: about 1.7 h). The densities beyond |z| = 0.75 A (0.4 % at 1 amu) need longer still; zone fractions beyond the barrier top are 1e-8 and cannot be measured by sampling at all.

---

## Stage 0c - does the conserved-quantity wander return at P = 64?

**What the earlier 0.10 kcal/mol was.** [read from the A3 pipeline and the A3 output] It is the span of i-PI's printed conserved quantity, which is `econs / nbeads`, i.e. **per bead** (x P for the total, about 3.3 kcal/mol at P = 32). The runs were on **Ag7_1** (7 Ag + 7 El, all El rings, **Ag classical atoms and free, not fixed**), not on the donor-acceptor system.
| case | omega_max dt (omega_max = 0.0786 P fs^-1) |
|---|---|
| P = 32, dt 0.25 fs (wandered) | 0.63 |
| P = 32, dt 0.125 fs (passed) | 0.31 |
| P = 64, dt 0.125 fs | 0.63 |
| P = 64, dt 0.0625 fs | 0.31 |
| P = 109 (planned 0.1 amu), dt 0.096 fs | 0.82 |
(Your estimates of 0.63 / 0.31 / 0.82 are confirmed.)

### Ag7_1 at P = 64 [measured] (same pipeline as A3: 6 ps PILE-L at dt 0.125, then NVE 2 ps from the same state)
| P | dt (fs) | span per bead | slow span (25 fs block means) | ripple rms (step to step) |
|---|---|---|---|---|
| 1 | 0.25 | 0.081 | 0.046 | 2.3e-3 |
| 8 | 0.125 | 0.069 | 0.029 | 1.4e-3 |
| 32 | 0.25 / 0.125 / 0.0625 | 0.297 / 0.151 / 0.122 | 0.102 / 0.102 / 0.101 | 3.5e-2 / 4.6e-3 / 9.7e-4 |
| **64** | **0.125** | **0.166** | **0.105** | 9.5e-3 |
| **64** | **0.0625** | **0.133** | **0.107** | 1.5e-3 |
(P = 1, 8 and 32 rows from `A3_results.json`.) **The wander returns at P = 64, with the same size as at P = 32 (0.105 vs 0.102), and does not change with the time step (0.105 vs 0.107).** The ripple falls with dt (6.5 x) as before. So (a) it does not scale with omega_max dt, as the time-step hypothesis would require (the P = 64, dt 0.0625 case has the P = 32, dt 0.125 value of 0.31 and wanders just as much), (b) it does not grow from P = 32 to 64 but is about 2 x its P = 1 value (0.046), (c) the per-bead span at P = 64 is 2.0 x (dt 0.125) the P = 1 span - below the factor-3 threshold. Cause still **not identified**. The Rg blocks (0.158-0.167 A, flat within 5 %) show the state was reasonably equilibrated, but the 6 ps equilibration is shorter than the 15 ps of A3.

### Production setup [measured] (Ag fixed, 4 classical El at 1 amu, 1 ring of mass 1 amu, 4 ps PILE-L, then NVE 2 ps from the same state; `stage0_nve.py`)
| P | dt (fs) | omega_max dt | span per bead | span x P | ripple rms | slow span (25 fs) | drift (kcal/mol/ps per bead) | T mean |
|---|---|---|---|---|---|---|---|---|
| 32 | 0.25 | 0.63 | 0.0399 | 1.28 | 2.2e-3 | 1.2e-3 | -1.5e-4 | 295.9 K |
| 32 | 0.125 | 0.31 | 0.0101 | 0.32 | 2.7e-4 | 2.6e-4 | +9e-6 | 295.9 K |
| 64 | 0.125 | 0.63 | 0.0097 | 0.62 | 3.1e-4 | 3.3e-4 | -1.7e-5 | 296.4 K |
| 64 | 0.0625 | 0.31 | 0.0023 | 0.15 | 3.9e-5 | 8.1e-5 | +5.7e-6 | 296.4 K |
The span falls as dt^2 (x 4.0 and x 4.2 on halving dt), and the slow part is 80-1000 times smaller than the 0.10 of Ag7_1. **No wander in the production setup at P = 32 or 64; the integration is better at the same omega_max dt than in Ag7_1 because only one electron is a ring and there are no free Ag atoms.** The P = 64 / dt 0.125 setting is therefore adequate for conserved-quantity purposes (span 0.62 kcal/mol in total, about 1 kT). P = 109 NVE was not run [not done]; P = 110 was only timed.

## Stage 0d - the production combination (Ag `fixatoms`, 4 `classical_atoms`, 1 ring)

Code check [read]: `Dynamics._setup_classical_atoms` refuses only atoms that are both fixed and classical (`"An atom cannot be both fixed and classical"`); the two sets are disjoint here (Ag fixed: ids 0-5; classical: ids 6-9; ring: id 10), so they combine. Bookkeeping [read]: for fixed atoms `pstep` applies no force (`activeatoms_mask`), so the only momentum `pconstraints` removes from them is initial or thermostat momentum, and i-PI's `fixatoms` booking (`eens += 0.5 p^2/m`) is the right one - the concern that it books differently from the classical-atom code does not apply to this combination, because the force kick is not applied to the fixed atoms in the first place. Check [measured]: the NVE runs in the table above (P = 32 and 64, Ag fixed + 4 classical + 1 ring) conserve the quantity (drift <= 1.5e-4 kcal/mol per ps per bead, i.e. <= 3e-4 kT per ps per bead), with the thermostat on during equilibration and off during the NVE segment. The booked quantity during thermostatted runs was not separately tested in a long NVT conserved-quantity run [not done].

---

## Stage 1 - runs

Setup for every row below: `run_barrier_lammps.py`, gap 4.0 A, Ag fixed, 4 classical El (1 amu) + one quantum electron of the stated mass, PILE-L (tau 100 fs, thermostat on all modes incl. centroid), 300 K, start = all beads collapsed in the left pit (z = -1.975 A), bead positions written every ~1 fs, bead density along z analysed with the first 2 ps (or 1 ps for 0.5 amu) discarded, 5-block standard errors (these underestimate the error when blocks are shorter than the correlation time; treat as indicative). P and dt: P from the pit rule of PART_B_REVISED, dt from omega_max dt = 0.35 (the NVE test above shows this is conserving).

### Results [measured]; z-density of the electron (all beads), compared with the 3-D Schroedinger solve (others frozen)
| run | length used | rms z (A) | P(\|z\|>0.5) | P(\|z\|>0.75) | P(\|z\|>1.0) | left zone (z<-1.7) |
|---|---|---|---|---|---|---|
| Schroedinger 3-D, 1 amu | | 0.2532 | 0.0448 | 0.0039 | 2.2e-4 | 4.5e-9 |
| classical P = 1, 1 amu (all five El classical), dt 0.25 | 10 ps | ring = point, stays at the pit: **1.0000 in the left zone, never beyond the top** | | | | 1.000 |
| **1 amu, P = 36, dt 0.125, pit start** | 38 ps (10 ps base + 30 ps continuation, first 2 ps dropped) | **0.2467 +- 0.0086** | 0.0437 +- 0.0066 | 0.0027 +- 0.0007 | 6.3e-5 +- 3.6e-5 | 0 |
| 1 amu, P = 72, dt 0.0625 (2P, same omega_max dt) | 12 ps | 0.2469 +- 0.0146 | 0.0459 +- 0.0119 | 0.0023 +- 0.0010 | 2.3e-6 +- 2.3e-6 | 0 |
| same P = 36 run, first 12 ps only (like for like with 2P) | 12 ps | 0.2529 +- 0.0089 | 0.0474 +- 0.0099 | 0.0026 +- 0.0012 | 6e-5 +- 6e-5 | 0 |
| 1 amu, P = 36, dt 0.0625 (dt/2) | 4 ps | 0.2246 +- 0.014 | 0.0249 +- 0.0095 | 5e-4 +- 4e-4 | 0 | 0 |
| 1 amu, P = 36, dt 0.125, **other four El frozen** (benchmark of the Schroedinger solve) | 10 ps | 0.2546 +- 0.010 | 0.0529 +- 0.0097 | 0.0022 +- 0.0006 | 8e-6 | 0 |
| Schroedinger 3-D, 0.5 amu | | 0.277 | 0.0658 | 0.0077 | 5.8e-4 | 2.2e-8 |
| **0.5 amu, P = 50, dt 0.09, pit start** | 9 ps | 0.299 +- 0.017 | 0.096 +- 0.019 | 0.0106 +- 0.0032 | 2.6e-4 +- 0.9e-4 | 0 |
Time to leave the pit [measured]: 1 amu, centroid beyond the barrier top at 0.83 ps; left-zone fraction in 0.5 ps windows 0.98, 0.64, 0, 0 ... for 10 ps (and for 38 ps). 0.5 amu: centroid beyond the top at 0.02 ps, left-zone fraction 0.04 in the first 0.5 ps, then 0 for 10 ps. (The very first beads move at once because the initial ring is collapsed and its internal momenta are thermal; the centroid time is the meaningful one.) Neither run ever returned to the pit region (zero left-zone frames after the departure; right zone also zero: nothing near the far pit at z > +1.7).

**Agreement with the prediction.** The frozen-others PIMD reproduces the Schroedinger solve (rms z 0.2546 +- 0.010 vs 0.2532, 0.1 sigma): this tests the 3-D solver, the engine and the sampling together, at one P and dt. With the other four electrons classical and moving (the production setup), 38 ps give 0.2467 +- 0.0086 (-0.8 sigma), P(\|z\|>0.5) 0.0437 vs 0.0448; the far tail P(\|z\|>1.0) is lower than predicted (6e-5 vs 2.2e-4, 4 sigma on the block errors) - in that tail the other electrons, which can move and are not frozen, differ from the frozen picture, or the 5-8 block errors are too small; not resolved. The 0.5 amu run is 1.3 sigma above the prediction in rms z and 1.6 sigma in P(\|z\|>0.5) with 9 ps and large errors.

**Convergence checks (1 amu).** P = 36 vs 72 (2P) at 12 ps each: rms z 0.253 +- 0.009 vs 0.247 +- 0.015; P(\|z\|>0.5) 0.047 +- 0.010 vs 0.046 +- 0.012: **no difference within +-6 %**. dt vs dt/2 (4 ps only): 0.2429 +- 0.005 (8 ps) vs 0.2246 +- 0.014: 1.2 sigma apart, not resolved; it would need >= 12 ps. These tests cannot show differences below about 5 %, and the 0.5 amu row has no P/2P or dt/2 check [not done].

**Cost actually spent [measured wall time on the loaded machine]:** 1 amu P = 36, 10 ps: 678 s; P = 72, 6 ps: 1470 s; 0.5 amu P = 50, 10 ps: 1226 s; 30 ps continuation of P = 36: 1900 s; frozen-others 12 ps: 730 s.

## What this means for the demonstration (needs your decision before Stage 2)

1. **In this landscape the electron is in the middle region at equilibrium for every mass (Stage 0b), and a PIMD ring started in the pit leaves it in 1 ps (1 amu) or 0.5 ps (0.5 amu).** The classical electron stays for the whole 10 ps (the classical crossing time estimate is 0.4 us). The "density beyond the barrier top" is therefore about 1 for 1 and 0.5 amu from the start - it cannot show a mass trend or a barrier-height trend, and the Gaussian boost on the tops (+5, +10) would only change how fast the ring leaves, not where it ends up (the equilibrium pit weight is 1e-8). I have **not** run 4 or 0.1 amu or any boost (Stage 2/3 not started).
2. What *is* measurable: (a) the width of the density (rms z 0.23, 0.25, 0.28, 0.39 A for 4, 1, 0.5, 0.1 amu, Schroedinger) and the tail beyond |z| = 0.75 A; (b) the pit-escape time from a pit start as a function of mass (0.02 ps at 0.5 amu, 0.83 ps at 1 amu, classical > 10 ps). (b) is a sampling-time statistic of ring-polymer dynamics with a thermostat, not a physical rate. 4 amu is the one mass where the ring (radius 0.10 A, smaller than the pit-to-top distance 0.275 A) may stay in the pit; the equilibrium pit weight there is 5e-5, but the dwell time is a Stage 2 measurement.
3. A design that would make the pit the equilibrium state (one Schroedinger-only check, no PIMD run, 1 amu, h = 0.04): add a Gaussian lift in the **middle** (z = 0, sigma 0.8 A) of height H. H = 5: pit still above the ground state by 5.2 kcal/mol; H = 10: 1.0 kcal/mol (zone fractions L / mid / R 0.055 / 0.89 / 0.057); **H = 15: pit degenerate with the ground state (0.095 kcal/mol), L / mid / R = 0.38 / 0.24 / 0.39 quantum, against 0.50 / 0.0000 / 0.50 classically.** There the density inside the barrier region (0.24 quantum vs 0 classical) is a clean penetration signal and the mass dependence should appear in it. This was not run at other masses and the height/width are untuned. The alternative of keeping the landscape as it is and reporting (a) and (b) is cheaper.

## Not done / limits
- 4 amu and 0.1 amu PIMD, boosts, middle-start sampling check (Stage 2/3, as instructed); middle start was not run at all.
- 0.5 amu has one run, 9 ps, no P/2P or dt/2 test; dt/2 at 1 amu only 4 ps.
- One seed per row. Block errors are from 5-8 blocks and are indicative only.
- Schroedinger: frozen other electrons; 4 amu only at h = 0.04 (pit width ~ h); 0.1 amu box-edge weight ~1e-3; 1-D and 3-D differ qualitatively (the 3-D one is the one to use).
- The cause of the Ag7_1 conserved-quantity wander is still unknown; it is absent in the production setup.
- P = 109 NVE conserved-quantity run not done. A long NVT conserved-quantity run of the production setup not done.

