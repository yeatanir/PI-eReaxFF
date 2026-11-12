# Higher barrier (V0 = 8 kcal/mol, sigma = 0.10 A): PIMD on the Mac Mini vs the 3-D solver

**Short verdict.** A classical electron is trapped on its start side (no crossings in 2 x 200 ps). The **0.5 amu** quantum electron crosses the barrier 20-21 times per 40 ps, passes the equilibration rule, and agrees with the solver within 5-8 % (below it). The **1 amu** electron is **not equilibrated** in 4 x 50 ps (3-6 centroid crossings per run; one run spent half its time in the left pit). Its density is reported per run and per start side, not as an equilibrium profile. Conserved-quantity noise is 5-10 x worse than for the 3 kcal/mol design, and no dt/2 or bead check was run on this design.

Labels: **[measured]** = from these runs; **[derived]** = computed from them; **[estimate]** = not measured.

## 1. What was run [measured]
Design: z-only Gaussian V0 = 8.0 kcal/mol, sigma = 0.10 A at z = 0 on the quantum electron (classical barrier B_cl = 12.19 kT, `HIGHER_BARRIER_SOLVER.md`). Ag fixed, four classical El (1 amu), one quantum El, PILE-L (tau 100 fs), 300 K, i-PI latency 1e-5 s, legacy control file (`CONTROL_USED NULL(legacy)`, logged in each `run_info.txt`), commits project 96b78d1 / i-PI 071f6b9 / metal-reax-lammps 120b609 on the Mac Mini tree `~/eb_stage3`, binaries md5-identical to the MacBook's (same-seed 1 amu check earlier: bitwise identical).

| Set | Runs | P | dt (fs) | Length | Start z (A) | Seeds | Speed (steps/s) |
|---|---|---|---|---|---|---|---|
| 1 amu | 4 | 72 | 0.125 | 50 ps | L -0.3 / R +0.3 alternating | 3001-3004 | 54-55 (83 for the last, alone) |
| 0.5 amu | 2 | 100 | 0.09 | 40 ps | L, R | 3011, 3012 | 42 |
| classical | 2 | 1 | 0.125 | 200 ps | L, R | 3021, 3022 | 1064 |

Position stride 200 (25 fs at 1 amu, 18 fs at 0.5 amu), 400 for classical (50 fs). First 2 ps discarded. **Not run:** 0.1 amu (cost), the frozen-others set (the solver freezes the other four electrons, PIMD lets them move: this stays a caveat), bead (2P) and dt/2 checks, an NVE-from-equilibrium run. Raw runs: `scratch/hb_mini/runs/` on the MacBook (copied from the Mini, checksum-verified, then removed from the Mini); the Mini is now empty of these runs.

## 2. Equilibration verdict [derived]
Rule fixed in advance (`analyze_higher_barrier.py::verdict`): equilibrated only if (a) every run has >= 10 centroid sign changes after the 2 ps skip, (b) every run first crosses within 25 % of its length, (c) left-start and right-start means of rho(0)/rho_max and of the straddle fraction agree within 20 %. Centroids are sampled every 18-25 fs, so crossings are a lower bound (fast recrossings are missed).

| Run | Start | Crossings (per ps) | First crossing (ps) | Centroid time left / right | Straddle | Pit frames (|centroid z| > 1.2 A) |
|---|---|---|---|---|---|---|
| **0.5 amu** 3011 | L | 21 (0.55) | 0.13 | 0.56 / 0.44 | 0.141 | 0.0 % |
| **0.5 amu** 3012 | R | 20 (0.53) | 0.83 | 0.54 / 0.46 | 0.139 | 0.05 % |
| **1 amu** 3001 | L | 6 (0.12) | 8.2 | 0.78 / 0.22 | 0.037 | 0.0 % |
| **1 amu** 3003 | L | 4 (0.08) | 0.13 | 0.80 / 0.20 | 0.041 | 0.0 % |
| **1 amu** 3002 | R | 3 (0.06) | 3.3 | 0.68 / 0.32 | 0.019 | **51 %** |
| **1 amu** 3004 | R | 5 (0.10) | 9.7 | 0.55 / 0.45 | 0.037 | 0.05 % |
| classical 3021 | L | 0 | never | 1.00 / 0.00 | 0 | 0.03 % |
| classical 3022 | R | 0 | never | 0.00 / 1.00 | 0 | 0 % |

- **0.5 amu: equilibrated** by the rule (21 and 20 crossings, both starts reach the other side within 1 ps, left/right 0.56/0.44 and 0.54/0.46, rho(0)/rho_max 0.126 vs 0.139, straddle 0.141 vs 0.139). Only two runs, so the run-to-run spread (n = 2) is a weak error estimate. Both runs spend slightly more time on the left (0.54-0.56); the solver is exactly symmetric, so this is either sampling or a real asymmetry from the moving classical electrons (not resolved).
- **1 amu: NOT equilibrated.** Reasons the rule reports: 3-6 crossings per run (< 10) and left-start vs right-start disagreement (rho(0)/rho_max 0.052 vs 0.039, straddle 0.039 vs 0.028). **Every 1 amu run spent more time on the left than the right (0.55-0.80), including the two right-start runs.** The 50 ps runs are too short for this barrier: the electron stays on one side for 5-20 ps at a time (Fig. 8d). I therefore report the 1 amu densities **per run and per start side** (Section 3) and **not** as an equilibrium profile.
- **Run 3002 entered the left pit.** The centroid first went beyond -1.2 A at 18.6 ps, left at 27.4 ps, re-entered at 28.4 ps and flickered in and out until 50 ps (51 % of post-2-ps frames in total, minimum z = -2.13 A). The solver puts the equilibrium weight of the pit state at ~1e-7 (the pit state is 16.4 kT above the ground state), so this is a kinetic trap of an unequilibrated run, like the 4 amu pit traps at the 3 kcal/mol design. The pit rule (fixed in advance) removes those frames from every histogram; **its entry time is 18.6 ps**. No other run has more than 0.05 % pit frames (3004 entered once at 11.6 ps for one frame, 3012 once at 33.8 ps).
- **Classical:** consistent with trapped. 0 crossings in 2 x 200 ps; the centroid came within 0.13 A (L) and 0.10 A (R) of z = 0. The classical crossing time is an **[estimate]** of >= 8 ns (`HIGHER_BARRIER_SOLVER.md`), so 400 ps total cannot prove trapping, only fail to contradict it.

## 3. Quantum vs classical and vs the solver [measured / derived]
Densities are bead-pooled z histograms (0.05 A bins, frames with |centroid z| > 1.2 A removed, renormalised over |z| < 1.5 A), compared with the 3-D solver (h = 0.03 A, frozen other electrons) for the same design (`data/higher_barrier.json`). Spread = std/sqrt(n) over independent runs (not block errors).

| Quantity | classical P=1 (2 runs) | 1 amu (4 runs) | 0.5 amu (2 runs) |
|---|---|---|---|
| rms z (A) PIMD (run-to-run) | 0.366 (0.002) | 0.457 (0.021) | 0.467 (0.0006) |
| rms z solver | 0.368 | 0.433 | 0.464 |
| P(|z|<0.1) PIMD | 1e-4 (0.0000, 0.0003 per run) | 0.0141 (0.0009) | 0.0288 (0.0001) |
| P(|z|<0.1) solver | 6.9e-5 | 0.0140 | 0.0304 |
| rho(0)/rho_max PIMD | 0.000 | 0.046 (0.008) | 0.133 (0.007) |
| rho(0)/rho_max solver | 0.0002 | 0.055 | 0.144 |

- **0.5 amu:** rms z agrees (+0.6 %). P(|z|<0.1) is 5 % below the solver and rho(0)/rho_max 8 % below. The differences are larger than the run-to-run spread, but that spread rests on two runs and on conserved-quantity noise (Section 4). I cannot say whether the deficit is time-step error, bead number, or the moving electrons the solver freezes: no dt/2 or 2P run exists for this design.
- **1 amu, per run [measured]:** rho(0)/rho_max = 0.050, 0.054, 0.057 (3001, 3003, 3004) and 0.021 (3002, which spent half its time in the pit); P(|z|<0.1) = 0.0140, 0.0150, 0.0157, 0.0116; rms z = 0.437, 0.433, 0.437 and 0.521 (3002). The three runs without pit time are within 5-10 % of the solver (0.055, 0.0140, 0.433), but this exclusion is post hoc and the crossing counts stay low, so it is not an equilibrium result. Left-start runs: rho(0)/rho_max 0.050, 0.054; right-start runs: 0.021 (pit run), 0.057.
- **Quantum vs classical [derived]:** the classical runs sample |z| < 0.1 A in only 1e-4 of frames (0 and about 1 of 3961 frames per run, resolution-limited), against 0.0141 (1 amu) and 0.0288 (0.5 amu): at least ~100x more barrier-region density for the quantum electron. The solver predicts 3.5e3 (1 amu) and 1.1e4 (0.5 amu) for the density ratio at z = 0, which the classical sampling here cannot resolve. In addition, the classical electron never crosses, while the 0.5 amu ring crosses ~0.5 times per ps and straddles the barrier 14 % of the time (1 amu: 2-4 % in the three non-pit runs).
- Straddle fraction = frames with beads on both sides of z = 0; the solver does not give this number.

## 4. Conserved quantity and speed [measured]
sim.out column 3 (per bead), after 2 ps, in kT (0.02585 eV):

| Run | std (kT) | linear drift (kT/ns) |
|---|---|---|
| 0.5 amu 3011 / 3012 | 0.26 / 0.19 | -17 / +6 |
| 1 amu 3001 / 3003 / 3002 / 3004 | 0.084 / 0.073 / 0.362 / 0.044 | +5.0 / -4.5 / +5.0 / -0.6 |
| classical 3021 / 3022 | 0.012 / 0.008 | 0.04 / 0.06 |

For comparison, the 3 kcal/mol design at 1 amu (P = 36, dt = 0.125 fs, 60 ps) had std 0.034 kT and drift 0.0015 kT/ns (`STAGE2_DESIGN.md`). At V0 = 8 kcal/mol the quantum runs are noisier, and the signed drifts (a 40-50 ps run drifts by about 0.2-0.7 kT) look like slow fluctuations rather than a steady trend. **dt = 0.125 fs (1 amu) and 0.09 fs (0.5 amu) were not re-validated for the narrower, taller Gaussian**; a dt/2 check should come before these densities are quoted to better than ~10 %. The classical runs conserve well.

Steps/s: 1 amu P = 72 about 54-55 with five jobs concurrent on the Mini (83 once alone), 0.5 amu P = 100 about 42, classical P = 1 about 1064. Wall time: 7370 s (1 amu, concurrent), 4800 s (last 1 amu run), ~10 500 s (0.5 amu), 1500 s (classical).

## 5. Disk footprint on the Mini [measured]
A finished 0.5 amu run directory is **6.2 MB** (ring.npy 2.67 MB + sim.out + logs; 7.5 MB on the MacBook with the log), a finished 1 amu directory 4.9 MB, a classical one ~12 MB. The raw bead xyz files peaked at ~150 MB per 0.5 amu run (691 B per frame per bead, 11 atoms; 112 MB at 75 % of the run) before the automatic reduction to ring.npy. For 0.1 amu (220 beads, 50 ps at dt 0.045 fs = 1.11 M steps) the raw peak would be ~845 MB per run at stride 200 or ~280 MB at stride 600 **[estimate]**, so 0.1 amu needs stride >= 600, <= 3-4 concurrent runs, or 10 ps restart segments. i-PI's trajectory output has no atom filter, so writing only the quantum electron's beads is not possible without an i-PI change.

## 6. Files
- Analysis: `inputs/tunneling_ag/analyze_higher_barrier.py` (usage in its docstring); summary `data/higher_barrier/hb_pimd.json`, centroid z(t) per run `data/higher_barrier/hb_centroids.npz`, run info `data/higher_barrier/run_info/`.
- Figure: `figures/fig08_higher_barrier.{py,png,pdf}` (a, b: density per run, solid = left start, dashed = right start, vs solver; c-e: centroid z(t)). The 1 amu panel is titled NOT equilibrated; read it as leakage from each side.
- Raw runs (not committed): `scratch/hb_mini/runs/hb_*`.
- Mac Mini regression check: `ipi_regression` hung at 0 % CPU on the Mini (reference files from cb9b4db present); left as instructed. Same-seed Mini-vs-MacBook run was done earlier (bitwise identical).
