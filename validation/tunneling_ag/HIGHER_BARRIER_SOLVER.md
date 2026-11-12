# Higher midpoint barrier: solver scan (no PIMD run)

Goal: find a z-only Gaussian barrier at z = 0 (on the quantum electron, over the frozen eReaxFF landscape) that is high enough that a **classical electron is trapped on one side** (classical barrier >= ~12 kT) while the **quantum electron still has appreciable density at z = 0 / is delocalised across both sides** at 0.1-1 amu.

Method [measured]: `inputs/tunneling_ag/schrod3d_barrier.py 0.03 14 data/higher_barrier.json 6:0.10 6:0.15 8:0.10 8:0.15 10:0.10 10:0.15 12:0.10 12:0.15 --masses 1,0.5,0.1` (3-D Schrödinger, h = 0.03 Å, 300 K, states up to 14 kT above the ground state, Ag fixed, other four El frozen; same code as the chosen-design scan). Raw results: `data/higher_barrier.json` (+ `.log`). Table made by `inputs/tunneling_ag/higher_barrier_report.py`. **4 amu was not run** (the h = 0.04 solve alone took 245 s on an idle machine; with 8 designs and a loaded machine it was not cheap).

Columns: B_cl PMF = classical barrier from the classical z-marginal incl. transverse entropy, V0/kT also given; margin = (W(0) − E0)/kT with W(0) the lowest transverse eigenvalue at z = 0 (transverse ZPE included); rho_q0/rho_c0 = thermal quantum / classical z-density at z = 0; rho_q0/rho_qmax = quantum density at z = 0 relative to its maximum (sampleability); split = lowest odd-z minus ground state; f(V>E) = thermal fraction of density in classically forbidden regions.
**Estimates (derived, not measured):** tau_xfer = πħ/split (isolated two-level transfer time of a left-localised start, ħ = 0.015179 ps·kcal/mol; ignores decoherence by the classical electrons and vibrations, so a *shortest* time); tau_cl = (1/0.42 ps) · exp(B_cl − 4.05), calibrated on the measured classical P=1 runs of the 3 kcal/mol design (0.42 centroid sign changes per ps including recrossings at B_cl = 4.05 kT), so a *lower bound* on the classical crossing time (Arrhenius scaling, equal prefactor assumed).

| V0 (kcal/mol) | sigma (A) | V0/kT | B_cl PMF (kT) | mass (amu) | E0 (kcal/mol) | margin below W(0) (kT) | rho_q0/rho_c0 | rho_q0/rho_qmax | split (kcal/mol) | split/kT | tau_xfer (ps, est) | tau_cl (ps, est) | f(V>E) thermal |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6 | 0.10 | 10.1 | 8.93 | 1 | 3.39 | 6.91 | 325 | 0.070 | 0.1209 | 0.203 | 0.394 | 313 | 0.12 |
| 6 | 0.10 | 10.1 | 8.93 | 0.5 | 4.34 | 6.41 | 794 | 0.210 | 0.3773 | 0.633 | 0.126 | 313 | 0.14 |
| 6 | 0.10 | 10.1 | 8.93 | 0.1 | 7.90 | 5.09 | 2186 | 0.840 | 2.0065 | 3.366 | 0.0238 | 313 | 0.13 |
| 6 | 0.15 | 10.1 | 8.12 | 1 | 3.85 | 6.21 | 75 | 0.033 | 0.0530 | 0.089 | 0.899 | 139 | 0.11 |
| 6 | 0.15 | 10.1 | 8.12 | 0.5 | 4.83 | 5.66 | 220 | 0.119 | 0.2214 | 0.371 | 0.215 | 139 | 0.13 |
| 6 | 0.15 | 10.1 | 8.12 | 0.1 | 8.40 | 4.31 | 898 | 0.761 | 1.6621 | 2.788 | 0.0287 | 139 | 0.13 |
| 8 | 0.10 | 13.4 | 12.19 | 1 | 3.51 | 10.03 | 3554 | 0.030 | 0.0652 | 0.109 | 0.731 | 8.14e+03 | 0.11 |
| 8 | 0.10 | 13.4 | 12.19 | 0.5 | 4.53 | 9.41 | 10632 | 0.108 | 0.2470 | 0.414 | 0.193 | 8.14e+03 | 0.12 |
| 8 | 0.10 | 13.4 | 12.19 | 0.1 | 8.21 | 7.89 | 45077 | 0.711 | 1.7264 | 2.896 | 0.0276 | 8.14e+03 | 0.15 |
| 8 | 0.15 | 13.4 | 11.31 | 1 | 4.00 | 9.29 | 564 | 0.010 | 0.0213 | 0.036 | 2.24 | 3.4e+03 | 0.11 |
| 8 | 0.15 | 13.4 | 11.31 | 0.5 | 5.05 | 8.63 | 2225 | 0.048 | 0.1185 | 0.199 | 0.402 | 3.4e+03 | 0.12 |
| 8 | 0.15 | 13.4 | 11.31 | 0.1 | 8.80 | 6.98 | 15925 | 0.571 | 1.3319 | 2.234 | 0.0358 | 3.4e+03 | 0.14 |
| 10 | 0.10 | 16.8 | 15.45 | 1 | 3.59 | 13.21 | 44696 | 0.014 | 0.0376 | 0.063 | 1.27 | 2.12e+05 | 0.10 |
| 10 | 0.10 | 16.8 | 15.45 | 0.5 | 4.66 | 12.52 | 155512 | 0.060 | 0.1675 | 0.281 | 0.285 | 2.12e+05 | 0.11 |
| 10 | 0.10 | 16.8 | 15.45 | 0.1 | 8.47 | 10.77 | 941237 | 0.584 | 1.4897 | 2.499 | 0.032 | 2.12e+05 | 0.16 |
| 10 | 0.15 | 16.8 | 14.56 | 1 | 4.11 | 12.44 | 5210 | 0.004 | 0.0095 | 0.016 | 5.05 | 8.78e+04 | 0.11 |
| 10 | 0.15 | 16.8 | 14.56 | 0.5 | 5.20 | 11.72 | 25953 | 0.022 | 0.0669 | 0.112 | 0.713 | 8.78e+04 | 0.11 |
| 10 | 0.15 | 16.8 | 14.56 | 0.1 | 9.12 | 9.78 | 288307 | 0.411 | 1.0704 | 1.795 | 0.0446 | 8.78e+04 | 0.15 |
| 12 | 0.10 | 20.1 | 18.73 | 1 | 3.65 | 16.43 | 623991 | 0.008 | 0.0228 | 0.038 | 2.09 | 5.64e+06 | 0.10 |
| 12 | 0.10 | 20.1 | 18.73 | 0.5 | 4.75 | 15.68 | 2454378 | 0.036 | 0.1172 | 0.197 | 0.407 | 5.64e+06 | 0.11 |
| 12 | 0.10 | 20.1 | 18.73 | 0.1 | 8.69 | 13.72 | 19801504 | 0.472 | 1.2902 | 2.164 | 0.037 | 5.64e+06 | 0.16 |
| 12 | 0.15 | 20.1 | 17.82 | 1 | 4.19 | 15.64 | 55741 | 0.002 | 0.0045 | 0.008 | 10.6 | 2.26e+06 | 0.11 |
| 12 | 0.15 | 20.1 | 17.82 | 0.5 | 5.30 | 14.87 | 339362 | 0.011 | 0.0396 | 0.066 | 1.2 | 2.26e+06 | 0.11 |
| 12 | 0.15 | 20.1 | 17.82 | 0.1 | 9.39 | 12.67 | 5308892 | 0.291 | 0.8650 | 1.451 | 0.0551 | 2.26e+06 | 0.15 |

For reference, the chosen 3 kcal/mol, σ = 0.10 Å design at 1 amu: B_cl = 4.05 kT, rho_q0/rho_c0 = 12, rho_q0/rho_qmax = 0.365 (`STAGE2_DESIGN.md`); classical electrons cross it thermally (~0.4 sign changes/ps).

## Reading the table
- **Classical trapping needs V0 >= 8 kcal/mol.** B_cl = 8.9 / 12.2 / 15.5 / 18.7 kT at V0 = 6 / 8 / 10 / 12 (σ = 0.10) and 8.1 / 11.3 / 14.6 / 17.8 kT (σ = 0.15). With the estimate above, tau_cl is ~0.3 ns at V0 = 6 (a 500 ps classical run would see ~1-2 crossings: not trapped), ≳ 8 ns at V0 = 8, σ = 0.10 (a 500 ps run crosses with probability ~6 %).
- **The quantum ground state is far below the effective barrier top everywhere** (margin 4-16 kT), so this is true under-barrier behaviour, not over-barrier. The splitting falls with V0 and mass: at 1 amu it is 0.12, 0.065, 0.038, 0.023 kcal/mol (V0 = 6, 8, 10, 12; σ = 0.10), i.e. 0.2-0.04 kT. The electron is **tunneling-coupled, not delocalised, at 1 amu** (split << kT); transfer time of a localised start 0.4-2 ps (est.). At 0.5 amu split = 0.38-0.12 kcal/mol (0.6-0.2 kT; tau_xfer 0.13-0.4 ps). At 0.1 amu split = 2.0-1.3 kcal/mol (3.4-2.2 kT) > kT: the state is delocalised across both sides at every V0 scanned, and the 0.1 amu rho_q0/rho_qmax stays 0.5-0.84 (σ = 0.10): the barrier is nearly invisible to the wavefunction even though it is 12-20 kT classically.
- **No value scanned localises the quantum electron with split < 1e-3 kT** (smallest split 0.0045 kcal/mol = 0.008 kT at V0 = 12, σ = 0.15, 1 amu, where tau_xfer ~ 11 ps, est.). So there is no 'no tunneling at this mass' case in this scan; but the largest V0 at 1 amu pushes rho_q0/rho_qmax to 0.002-0.008, which PIMD cannot sample.
- **Quantum / classical density at z = 0 is huge** (3.5e3 at 1 amu, V0 = 8, σ = 0.10; 1e4 at 0.5 amu; 4.5e4 at 0.1 amu) because the classical density there is exponentially small: this is the clean quantum-vs-classical contrast the 3 kcal/mol design lacks (12 at 1 amu).
- **f(V>E) is 0.10-0.16 for all entries**: the generic zero-point forbidden fraction, not a tunneling-specific number.

## Sampleability by PIMD
rho_q0/rho_qmax at V0 = 8, σ = 0.10: 0.030 (1 amu), 0.108 (0.5 amu), 0.71 (0.1 amu). The 3 kcal/mol design had 0.365 at 1 amu and was resolved to ±0.04 with ~20 ps pooled; 0.03 means about 12× fewer samples at the barrier, so 1 amu would need ~10× longer runs for the same relative error on the density at z = 0 (but centroid/ring straddling observables converge faster). V0 = 8, σ = 0.15 is worse (0.010 at 1 amu). V0 = 6, σ = 0.10 is the most sampleable classically-marginal choice (0.070, 0.21, 0.84).

## Recommendation (nothing launched)
**Run next: V0 = 8 kcal/mol, σ = 0.10 Å.** B_cl = 12.2 kT (V0/kT = 13.4), classical electron trapped (tau_cl ≳ 8 ns est.), while 0.1 amu is delocalised over both sides (split 1.7 kcal/mol = 2.9 kT, rho_q0/rho_qmax = 0.71), 0.5 amu is tunneling-coupled with a visible barrier-region density (0.108, tau_xfer ~ 0.2 ps est.) and 1 amu is the deep-tunneling case (0.03, split 0.11 kT; needs ~10× longer runs). Use P and dt per `STAGE2_DESIGN.md`, and the 0.5 amu P=100 setting. Classical P=1 comparison runs of 500 ps from both wells should now stay trapped, which is the contrast the animations of the 3 kcal/mol design lack. If sampling at 1 amu proves too costly, fall back to V0 = 6 kcal/mol, σ = 0.10 Å (B_cl = 8.9 kT, classical crossing ~0.3 ns, rho_q0/rho_qmax 0.07), which is only marginal for classical trapping.
Caveats: solver freezes the other four El (PIMD lets them move: ~9-25 % tail excess seen at 1 amu); h = 0.03 Å, emax = 14 kT (pit states lie above that and are not needed here); 0.1 amu box-size error ~1 %; tau_cl and tau_xfer are estimates.

## The system (answers to the questions)
- **Atoms:** 6 Ag (type Av) and 5 El: Ag at z = ∓2.0 Å (y = 0, the inner donor/acceptor sites) and at z = ∓4.605 Å, y = ±1.411 Å (two outer triangles' remaining atoms, four in total); the four outer Ag carry the four classical El (placed on them, +0.001 Å in x); the fifth El is the single quantum electron, starting at z = −0.3 Å (init.xyz). 11 atoms.
- **Charge, verified:** the LAMMPS `data.lmp` of the Stage 2 run (`scratch/tun_stage2/s2_m1_P36_dt125_L_s1/client0/data.lmp`) has six atoms at +1.0 (all Ag) and five at −1.0 (all El): **net charge +1**, as inferred. Charges are fixed (pair style in `ereaxff yes` mode reads them from the data file; `icharg 8`-style fixed charges). The two inner Ag are therefore bare Ag⁺ with no El of their own; the quantum electron is the only electron shared between them.
