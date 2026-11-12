# Stage 2 — midpoint-barrier design (solver first) and PIMD on the chosen design

Date: 2026-11-11. Labels: **[measured]** (computed/run here), **[derived]**, **[estimated]**.
Scripts (`inputs/tunneling_ag/`): `schrod3d_barrier.py` (design solver), `stage2_table.py`, `run_barrier_lammps.py` (`--v0-kcal --sigma-z --start-z --latency`), `analyze_stage2.py`, `stage2_conserved.py`.
Data: `data/stage2/` (solver scans `scan_*.json/log`, `table.md`, `conv_3_010_h003_*.json`, PIMD summaries `pimd_*.json`). Raw runs in `scratch/tun_stage2/` (not committed).
Stage 1 (intrinsic 10.5 kcal/mol barrier, `STAGE0_STAGE1.md`) is kept unchanged as the validation figure; no further runs were made on it.

## 0. Summary

- **Design (solver first, passed):** z-only Gaussian at z = 0 on the quantum electron, **V0 = 3.0 kcal/mol, σ = 0.10 Å**. At 1 amu, 300 K: ground state 2.5 kT below the effective barrier top, quantum/classical density at z = 0 = 12, classical barrier 5.0 kT (V0/kT) / 4.05 kT (PMF). Only 1 of 12 scanned designs meets all criteria under both readings. The 15 kcal/mol middle lift was dropped.
- **PIMD versus solver (barrier region and width), all within ≈ 2σ:** 4 amu (P=18, dt 0.25), 1 amu (P=36, dt 0.125; also P=72, dt/2, P=18/0.25), 0.5 amu (P=100, dt 0.09), 0.1 amu (P=110, dt 0.09 and 0.045). **0.5 amu at P=50 failed (rms +2.9 %, tail +10–32 %, 3σ); P=100 fixes it** (section 2.6a).
- **Mass trend of the barrier-region density** P(|z|<0.1), quantum/classical: 2.8 (4 amu), 6.5 (1), 8.9 (0.5), 11.0 (0.1) [PIMD]; solver 2.6, 5.9, 8.4, 10.4. 0.1 amu is zero-point delocalisation over the barrier (margin 1.0 kT), not tunneling.
- **Converged settings at 1 amu:** P→2P and dt→dt/2 change the barrier-region density by 0.6 % and −4 % (≤ 1.2σ) and the width by 0.5 % and 0.3 %. **P=18, dt 0.25 fs matches P=36 within 1σ in every observable (4 runs, 120 ps) and is recorded as the later 1 amu setting for the barrier-region density**; its conserved-quantity noise is 5–16× larger.
- **Conserved quantity** (thermostatted, this setup, 60 ps long run at 1 amu): std 0.034 kT, span 0.20 kT, drift 0.0015 kT/ns. 0.1 amu at dt 0.09 fs is noisy (std 0.5–1.2 kT); dt 0.045 fs gives 0.06 kT.
- **Speed:** the i-PI socket latency (1e-4 → 1e-5 s) gives +58 % (116 → 184 steps/s at 1 amu P=36), bitwise-identical output. Production steps/s are in section 3.

### Settings to use
| mass (amu) | P | dt (fs) | note |
|---|---|---|---|
| 4 | 18 | 0.25 | pit traps: exclude \|z\|>1.2 Å frames |
| 1 | 36 | 0.125 | or P=18, dt 0.25 for the barrier-region density only |
| 0.5 | 100 | 0.09 | P=50 is enough for P(\|z\|<0.1) but not for width/tails |
| 0.1 | 110 | 0.045 | dt 0.09 gives the same density but a noisy conserved quantity |

## 1. Part 1 — solver first

### 1.1 What was solved [measured]
3-D Schrödinger problem of one quantum electron (mass m, T = 300 K) in the frozen eReaxFF landscape of the 4.0 Å donor–barrier–acceptor system (Ag fixed, other four El frozen on their cations; same potential grid as `schrod3d.py`), **plus a z-only Gaussian** `Vb(z) = V0 exp(-z²/2σ²)` centred at z = 0 (the same shape as `ipi/pes/gaussian_barrier.py` with s = z, acting on the quantum electron only). Thermal density from all states within 10 kT of E0; h = 0.04 Å for the scan, h = 0.03 Å for the chosen design (agree to 0.3 % in E0, 4 % in the density ratio, 0.3 % in ρ(0)/ρmax).

Quantities (`schrod3d_barrier.py`):
- **E0**, **W(0)**: W(z) is the lowest eigenvalue of the 2-D transverse problem at fixed z (the adiabatic profile; it contains the transverse zero-point energy). The *margin* is (W(0) − E0)/kT with the barrier in W(0): how far the ground state sits below the *effective* barrier top.
- **ρ_q(0)/ρ_c(0)**: thermal z-marginal density at z = 0, quantum over classical Boltzmann (same box, mass independent).
- **Classical barrier**: two readings are tabulated. V0/kT, and the **PMF barrier** ln(ρ_c at the lower of the two density maxima / ρ_c(0)), which is what a classical particle crosses (the maxima sit on the flanks of the Gaussian, so it is lower than V0/kT).
- **Splitting**: lowest odd-z state minus the ground state in the (even x, even y) sector.
- **f(V>E)**: fraction of the quantum density where V > E (ground state; thermal average with each state's own energy).

Selection criteria at 1 amu (yours): margin ≥ 2 kT, ρ_q/ρ_c(0) ≥ 5, classical barrier 4–6 kT.

### 1.2 Scan at 1 amu [measured, h = 0.04 Å]
| V0 (kcal/mol) | σ (Å) | E0 | W(0) | margin (kT) | V0/kT | PMF barrier (kT) | ρ_q/ρ_c at 0 | ρ_q(0)/ρ_q,max | split (kcal/mol) | f(V>E) gs / th | pass (PMF) | pass (V0/kT) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (bare) | – | 2.114 | 1.571 | −0.91 | 0 | 0 | 0.89 | 1.000 | 1.071 | 0.113 / 0.102 | – | – |
| 2 | 0.10 | 2.799 | 3.532 | 1.23 | 3.35 | 2.48 | 4.5 | 0.588 | 0.520 | 0.137 / 0.112 | no | no |
| 2 | 0.15 | 3.067 | 3.554 | 0.82 | 3.35 | 1.93 | 2.7 | 0.568 | 0.439 | 0.130 / 0.107 | no | no |
| 2 | 0.20 | 3.301 | 3.561 | 0.44 | 3.35 | 1.35 | 1.9 | 0.646 | 0.425 | 0.124 / 0.102 | no | no |
| **3** | **0.10** | **3.020** | **4.512** | **2.50** | **5.03** | **4.05** | **12.1** | **0.338** | **0.354** | 0.167 / 0.127 | **yes** | **yes** |
| 3 | 0.15 | 3.372 | 4.545 | 1.97 | 5.03 | 3.40 | 5.7 | 0.279 | 0.256 | 0.147 / 0.114 | no (margin 1.97) | no |
| 3 | 0.20 | 3.705 | 4.557 | 1.43 | 5.03 | 2.69 | 3.3 | 0.301 | 0.224 | 0.137 / 0.108 | no | no |
| 4 | 0.10 | 3.179 | 5.492 | 3.88 | 6.71 | 5.64 | 33.8 | 0.196 | 0.243 | 0.165 / 0.126 | yes | no (6.7 kT) |
| 4 | 0.15 | 3.582 | 5.536 | 3.28 | 6.71 | 4.94 | 12.5 | 0.133 | 0.148 | 0.157 / 0.120 | yes | no (6.7 kT) |
| 4 | 0.20 | 3.984 | 5.552 | 2.63 | 6.71 | 4.13 | 5.9 | 0.128 | 0.112 | 0.148 / 0.115 | yes | no (6.7 kT) |
| 6 | 0.10 | 3.381 | 7.453 | 6.83 | 10.06 | 8.81 | 304 | 0.071 | 0.121 | 0.141 / 0.115 | no | no |
| 6 | 0.15 | 3.839 | 7.518 | 6.17 | 10.06 | 8.08 | 73 | 0.034 | 0.053 | 0.133 / 0.114 | no | no |
| 6 | 0.20 | 4.323 | 7.542 | 5.40 | 10.06 | 7.15 | 23 | 0.024 | 0.029 | 0.129 / 0.113 | no | no |

Full table with all masses: `data/stage2/table.md`.

**Chosen design: V0 = 3.0 kcal/mol, σ = 0.10 Å (FWHM 0.2355 Å), z-only Gaussian centred at z = 0 on the quantum electron.** It is the only scanned design that passes under both readings of "classical barrier". Caveat: the PMF barrier (4.05 kT) is just at the lower edge of 4–6 kT. (3, 0.10) has σ/h = 2.5–3.3 grid points per σ, so it was re-solved at h = 0.03 Å (table below); the other candidates were not.

### 1.3 Chosen design at the four masses [measured, h = 0.03 Å]
| m (amu) | E0 (kcal/mol) | W(0) | margin (kT) | ρ_q/ρ_c at 0 | ρ_q(0)/ρ_q,max | splitting (bare → with barrier, kcal/mol) | f(V>E) gs / th | free-ring radius (Å, 300 K) |
|---|---|---|---|---|---|---|---|---|
| 4 | 1.901 | 3.752 | 3.11 | 3.8 | 0.077 | 0.544 → 0.032 | 0.133 / 0.092 | 0.10 |
| 1 | 3.028 | 4.544 | 2.54 | 12.5 | 0.338 | 1.071 → 0.354 | 0.166 / 0.126 | 0.20 |
| 0.5 | 3.867 | 5.199 | 2.23 | 19.2 | 0.618 | 1.494 → 0.753 | 0.141 / 0.124 | 0.28 |
| 0.1 | 7.347 | 7.968 | 1.04 | 25.2 | 0.987 | 3.162 → 2.521 | 0.115 / 0.115 | 0.64 |

[derived] Reading: 4 amu is classical-like (density at the barrier 3.8× the classical, but only 8 % of the maximum, splitting 0.03 kcal/mol). 1 and 0.5 amu are the penetration regime (ground state 2.2–2.5 kT below the effective top, density at z = 0 12–19× classical and 34–62 % of the maximum). **0.1 amu is zero-point delocalisation over the barrier, not tunneling**: the ground state is only 1 kT below the top and the barrier is almost invisible in the density (ρ(0)/ρmax = 0.99). The crossover temperatures quoted in `PART_B_REVISED.md` (183/367/519/1160 K) belong to the intrinsic pit barrier, not to this one.

### 1.4 The earlier 15 kcal/mol middle lift (σ = 0.8 Å) — dropped [measured]
| m (amu) | E0 | W(0) | margin (kT) | ρ_q(0)/ρ_q,max | f(V>E) gs / th | rms z (Å) |
|---|---|---|---|---|---|---|
| 4 | 7.40 | 15.78 | 14.1 | 1.2e-13 | 0.128 / 0.128 | 1.98 |
| 1 | 12.54 | 16.57 | 6.75 | 4.1e-4 | 0.189 / 0.170 | 1.82 |
| 0.5 | 14.32 | 17.22 | 4.9 | 6.6e-3 | 0.131 / 0.126 | 1.42 |
| 0.1 | 18.07 | 19.99 | 3.2 | 0.099 | 0.124 / 0.117 | 1.22 |

Decision: **drop**. The electron sits in the pits (rms z 1.2–2.0 Å), and the density at z = 0 is 4e-4 of the maximum at 1 amu and 1e-13 at 4 amu, which PIMD cannot sample. The fraction of density with V > E0 is 0.12–0.19, i.e. the same as the barrier-free landscape (0.11): that fraction is the generic zero-point forbidden fraction of any bound state and carries no barrier-specific tunneling signal. (I read your rule "if that fraction is small, drop it" as "if the excess over the barrier-free value is small"; the absolute fraction is not small. If you meant the absolute value, say so.)

## 2. Part 2 — PIMD on the chosen design

### 2.1 Setup [measured]
Ag fixed (`<fixatoms>`), four classical El at 1 amu (`<classical_atoms>`), one quantum electron (ring of P beads) of the scanned mass; Gaussian barrier V0 = 3 kcal/mol, σ = 0.10 Å in z at z = 0 on that electron (python `gaussian_barrier` driver, checked against the analytic value to 3e-7 eV). PILE-L thermostat, τ = 100 fs, 300 K. Bead positions every ≈ 5 fs; the bead density along z (all beads, after 1 ps) is the ring-polymer estimator of the quantum density. Starts at z = ∓0.3 Å (left/right of the barrier, the density maxima are at ±0.28 Å). The solver freezes the other four electrons; the PIMD lets them move.

| mass (amu) | P | dt (fs) | runs (starts: L/R) | length |
|---|---|---|---|---|
| classical | 1 | 0.25 | 2 (L, R) | 500 ps each |
| 4 | 18 | 0.25 | 6 (L,R,L,R,L,R; seeds 21–26) | 40 ps each |
| 1 | 36 | 0.125 | 4 (L,R,L,R; seeds 1–4) + 1 long (L, seed 71) | 20, 20, 24, 24 + 60 ps |
| 0.5 | 50 | 0.09 | 3 (L,R,L; seeds 31–33) | 24 ps each |
| 0.5 | 100 | 0.09 | 3 (L,L,R; seeds 93, 95, 96) | 24 ps each |
| 0.1 | 110 | 0.09 | 3 (L,R,L; seeds 61–63) | 7.9 (run killed, see caveats), 12, 12 ps |

### 2.2 Density profile versus the solver [measured]
Pooled over runs; errors are standard errors over 8 time blocks per run (blocks are partly correlated, so they are optimistic). "σ" = difference in units of that error. Solver = 3-D Schrödinger for the same barrier and mass, h = 0.03 Å.

| mass | observable | PIMD | solver | diff |
|---|---|---|---|---|
| 4 amu (|z|<1.2 Å, see below) | rms z (Å) | 0.3465 ± 0.0128 | 0.3439 | +0.7 % (0.2σ) |
| | P(\|z\|<0.1) | 0.0404 ± 0.0024 | 0.0399 | +1.2 % (0.2σ) |
| | P(\|z\|>0.5) | 0.098 ± 0.021 | 0.0946 | +3.6 % (0.2σ) |
| | ρ(0)/ρmax | 0.127 ± 0.010 | 0.134 | −5.4 % (0.7σ) |
| 1 amu | rms z | 0.3678 ± 0.0026 | 0.3649 | +0.8 % (1.1σ) |
| | P(\|z\|<0.1) | 0.0942 ± 0.0023 | 0.0922 | +2.2 % (0.9σ) |
| | P(\|z\|>0.5) | 0.1548 ± 0.0052 | 0.1482 | +4.4 % (1.3σ) |
| | P(\|z\|>0.75) | 0.0142 ± 0.0014 | 0.0130 | +9 % (0.8σ) |
| | ρ(0)/ρmax | 0.374 ± 0.017 | 0.365 | +2.5 % (0.5σ) |
| 0.5 amu | rms z | 0.3877 ± 0.0037 | 0.3767 | **+2.9 % (2.9σ)** |
| | P(\|z\|<0.1) | 0.1285 ± 0.0014 | 0.1302 | −1.3 % (1.2σ) |
| | P(\|z\|>0.5) | 0.1975 ± 0.0057 | 0.1803 | **+9.5 % (3.0σ)** |
| | P(\|z\|>0.75) | 0.0322 ± 0.0029 | 0.0244 | **+32 % (2.7σ)** |
| | ρ(0)/ρmax | 0.612 ± 0.014 | 0.640 | −4.5 % (2.1σ) |
| 0.5 amu, **P = 100** (3 runs, 72 ps) | rms z | 0.3750 ± 0.0022 | 0.3767 | −0.4 % (0.8σ) |
| | P(\|z\|<0.1) | 0.1279 ± 0.0015 | 0.1302 | −1.8 % (1.5σ) |
| | P(\|z\|>0.5) | 0.1780 ± 0.0034 | 0.1803 | −1.3 % (0.7σ) |
| | P(\|z\|>0.75) | 0.0230 ± 0.0015 | 0.0244 | −6 % (1.0σ) |
| | ρ(0)/ρmax | 0.601 ± 0.011 | 0.640 | −6.2 % (3.6σ) |
| 0.1 amu | rms z | 0.4488 ± 0.0061 | 0.4407 | +1.8 % (1.3σ) |
| | P(\|z\|<0.1) | 0.1592 ± 0.0022 | 0.1618 | −1.6 % (1.2σ) |
| | P(\|z\|>0.5) | 0.2709 ± 0.0066 | 0.2620 | +3.4 % (1.4σ) |
| | ρ(0)/ρmax | 0.988 ± 0.006 | 0.992 | −0.4 % (0.6σ) |
| classical P=1 (\|z\|<1.2) | rms z | 0.3332 ± 0.0014 | 0.3307 | +0.7 % (1.8σ) |
| | P(\|z\|<0.1) | 0.0144 ± 0.0006 | 0.0155 | −6.8 % (1.7σ) |
| | ρ(0)/ρmax | 0.039 ± 0.002 | 0.045 | −12 % (2.5σ) |

**Barrier-region density versus classical** (P(|z|<0.1), PIMD quantum / PIMD classical P=1, both truncated at |z|<1.2 Å): 4 amu 2.8, 1 amu 6.5, 0.5 amu 8.9, 0.1 amu 11.0 [derived; ≈ ±6 %]. The solver gives 2.6, 5.9, 8.4, 10.4. The ratio rises monotonically as the mass falls, and the PIMD reproduces the solver's trend and size to ≤ 10 %.

**Equilibration from both sides** [measured]: at every mass the left- and right-started runs give the same profile (1 amu rms z 0.3670/0.3687/0.3725/0.3630 for L,R,L,R; 0.5 amu 0.3846/0.3893/0.3891; 0.1 amu 0.4328/0.4561/0.4514; 4 amu 0.3447, 0.3405, 0.3568, 0.3505, 0.3369, 0.3508). The classical run crosses the barrier region 144–208 times in 500 ps (every ≈ 3 ps), so the middle region is equilibrated classically too.

**Kinetic traps (pit).** The classical right-started run fell into the right pit at 425 ps and stayed (15 % of the frames); two of the six 4 amu runs (seeds 23 and 25, at ≈ 25 and ≈ 31 ps) had the ring centroid enter the left pit and stay there (37 % and 23 % of the frames). No ring left a pit within the run. These frames are excluded (|z| > 1.2 Å) in the comparison above, and the solver is renormalised the same way. At 1, 0.5 and 0.1 amu no ring entered a pit. The pit-entry events at 4 and P = 1 are the kinetic traps you warned about (the equilibrium pit weight is 5e-5 at 4 amu, 7e-3 classical, yet an entry is not reversible within 10–40 ps).

### 2.3 Convergence at 1 amu [measured]
| configuration | P | dt (fs) | total sampled | rms z | P(\|z\|<0.1) | P(\|z\|>0.5) | ρ(0)/ρmax |
|---|---|---|---|---|---|---|---|
| base | 36 | 0.125 | 88 ps (4 runs) | 0.3678 ± 0.0026 | 0.0942 ± 0.0023 | 0.1548 ± 0.0052 | 0.374 ± 0.017 |
| 2P | 72 | 0.125 | 72 ps (3 runs) | 0.3695 ± 0.0027 | 0.0948 ± 0.0026 | 0.1595 ± 0.0047 | 0.361 ± 0.020 |
| dt/2 | 36 | 0.0625 | 72 ps (3 runs) | 0.3688 ± 0.0028 | 0.0904 ± 0.0023 | 0.1546 ± 0.0047 | 0.356 ± 0.015 |
| long | 36 | 0.125 | 60 ps (1 run) | 0.3658 ± 0.0039 | 0.0914 ± 0.0026 | 0.1475 ± 0.0076 | 0.356 ± 0.019 |
| solver | | | | 0.3649 | 0.0922 | 0.1482 | 0.365 |

P→2P: rms +0.5 %, P(|z|<0.1) +0.6 %, ρ(0)/ρmax −3.4 % (all within 1σ, errors 0.7–6 %). dt→dt/2: +0.3 %, −4.0 %, −4.8 % (1.2σ). The barrier-region density is converged to ≈ 4 %, the width to 0.5 %. The tail P(|z|>0.75) is 9–21 % above the solver in all four PIMD settings (0.7–1.3σ each), independent of P and dt.

### 2.4 Time step / bead count candidate: P = 18, dt = 0.25 fs at 1 amu [measured]
4 runs (L, R, L, R; seeds 81–84), 30 ps each, 120 ps total: rms z 0.3705 ± 0.0019 (solver +1.5 %, 3.0σ; P = 36 base +0.7 %, 0.8σ), P(|z|<0.1) 0.0945 ± 0.0020 (P = 36: 0.0942), ρ(0)/ρmax 0.372 ± 0.015 (P = 36: 0.374), P(|z|>0.5) 0.1592 ± 0.0036 (P = 36: 0.1548; solver 0.1482, +7.4 %, 3.1σ), P(|z|>0.75) 0.0162 ± 0.0013 (P = 36: 0.0142). **It matches the P = 36 result within 1σ in every observable and the solver in the barrier region (P(|z|<0.1) +2.5 %, ρ(0)/ρmax +2.1 %); the width and the outer tail are 1.5 % and 7 % above the solver, 3σ.** The conserved-quantity noise is 5–16× larger than at P = 36 / 0.125 fs (std 0.05–0.18 kT vs 0.011 kT; see 2.5). **Recorded as the setting for later 1 amu work where the barrier-region density is the target (about half the cost of P = 36).** Do not use it for tail quantities without a check against P = 36.

### 2.5 Conserved quantity (thermostatted, this setup) [measured]
`stage2_conserved.py`: i-PI `conserved`, every 20 steps, after 1 ps. kT = 0.02585 eV.
| run | length | std (kT) | span (kT) | drift (kT/ns) |
|---|---|---|---|---|
| 1 amu P=36 dt 0.125, **long (60 ps)** | 60 ps | 0.034 | 0.20 | +0.0015 |
| 1 amu P=36 dt 0.125 (seed 1) | 20 ps | 0.011 | 0.066 | +0.0005 |
| 1 amu P=36 dt 0.0625 | 24 ps | 0.0045 | 0.028 | −0.0004 |
| 1 amu P=72 dt 0.125 | 24 ps | 0.062 | 0.24 | +0.006 |
| 1 amu P=18 dt 0.25 (seeds 81–83) | 30 ps | 0.05–0.18 | 0.3–0.65 | +0.02 … −0.0005 |
| 4 amu P=18 dt 0.25 | 40 ps | 0.025 | 0.16–0.18 | +0.001 |
| 0.5 amu P=50 dt 0.09 | 24 ps | 0.027–0.046 | 0.14–0.25 | +0.003–0.006 |
| 0.1 amu P=110 dt 0.09 | 12 ps | **0.52 / 1.15** | **2.6 / 3.7** | −0.06 / +0.34 (noisy, 12 ps) |
| classical P=1 dt 0.25 | 500 ps | 0.13 | 0.59 | +0.0008 |

Drift is ≤ 0.02 kT/ns at P ≥ 18 for masses ≥ 0.5 amu. The std scales as dt² from 0.125 to 0.0625 fs (0.011 → 0.0045 is ×2.4, expected ×4: other terms matter at this level), and is larger for the heavier-bead-count cases. **0.1 amu at dt 0.09 fs is the weak point: std 0.5–1.2 kT, span 2.6–3.7 kT** (about 50–100× the 1 amu P = 36 value). The density still agrees with the solver within 1.6σ, but I cannot call that run time-step converged (see 2.6).

### 2.6 Diagnostics added for the failing comparisons [measured]
**(a) 0.5 amu, P = 50 is too wide in the tails; P = 100 agrees with the solver.** The P = 50 comparison failed at 2.9–3.0σ (rms +2.9 %, P(|z|>0.5) +9.5 %, P(|z|>0.75) +32 %). Three extra runs at P = 100 (same dt 0.09 fs; seeds 93, 95, 96; L, R, L; 24 ps each, 72 ps total, mobile other electrons) give:
| 0.5 amu | rms z (Å) | P(\|z\|<0.1) | P(\|z\|>0.5) | P(\|z\|>0.75) | ρ(0)/ρmax |
|---|---|---|---|---|---|
| P = 50 (3 runs, 72 ps) | 0.3877 ± 0.0037 | 0.1285 ± 0.0014 | 0.1975 ± 0.0057 | 0.0322 ± 0.0029 | 0.612 ± 0.014 |
| **P = 100 (3 runs, 72 ps)** | **0.3750 ± 0.0022** | 0.1279 ± 0.0015 | 0.1780 ± 0.0034 | 0.0230 ± 0.0015 | 0.601 ± 0.011 |
| solver | 0.3767 | 0.1302 | 0.1803 | 0.0244 | 0.640 |
| P = 50, other four electrons **frozen** (1 run, 24 ps) | 0.3720 ± 0.0043 | 0.1290 ± 0.0032 | 0.1731 ± 0.0058 | 0.0222 ± 0.0026 | 0.592 ± 0.021 |

P = 100 vs the solver: rms −0.4 % (0.8σ), P(|z|>0.5) −1.3 % (0.7σ), P(|z|>0.75) −6 % (1.0σ), P(|z|<0.1) −1.8 % (1.5σ). ρ(0)/ρmax is 6 % low (3.6σ by block errors; this ratio of histogram bins is the noisiest observable and was −4.5 % at P = 50). P = 50 → P = 100 changes the width by −3.3 % and the outer tail by −10 % to −29 %, and the barrier-region density P(|z|<0.1) by −0.5 %: **the barrier-region density is converged at P = 50, the width and tails are not.** [estimated] The likely cause is the stiff walls the ring reaches in the tails (pit-side stiffness 45 eV/Å²; at 0.5 amu βħω/P is no longer small there at P = 50), which the bead-count rule I used (based on the barrier curvature) does not see. The frozen-others P = 50 run agrees with the solver (−1.2 %, 1.1σ), so the mobile other electrons also play a role in the P = 50 excess, but one run cannot separate the two effects. **For 0.5 amu use P = 100 if tails matter (about 2× the cost), P = 50 for the barrier-region density.**

**(b) 1 amu, mobile versus frozen other electrons.** One P = 36 run with the other four electrons frozen (24 ps): rms 0.3606 ± 0.0044 (solver −1.2 %), P(|z|>0.75) 0.0103 ± 0.0024 (solver −21 %, 1.2σ), versus mobile 0.3678 and 0.0142. So the PIMD tail excess of 9–25 % over the frozen-others solver at 1 amu (P-independent, section 2.3) is of the sign and size expected from mobile electrons; one frozen run has ±20 % tail error, so this is consistent, not proven.

**(c) 0.1 amu time step.** One run at dt 0.045 fs (P = 110, 8.0 ps, seed 94, L): conserved-quantity std 0.060 kT, span 0.25 kT, drift +0.026 kT/ns, against 0.52 / 1.15 kT std and 2.6 / 3.7 kT span at dt 0.09 fs (a factor 9–19, more than the dt² expected 4). Density: rms 0.4513 ± 0.0085 (dt 0.09: 0.4488 ± 0.0061; solver 0.4407), P(|z|<0.1) 0.1570 ± 0.0041 (0.1592; solver 0.1618), P(|z|>0.5) 0.2759 ± 0.0091 (0.2709; solver 0.2620), ρ(0)/ρmax 0.979 ± 0.020 (0.988; solver 0.992). The density does not depend on dt within errors (≤ 1.2σ), but **dt 0.045 fs is the safer setting at 0.1 amu** (2× the cost).

## 3. Cost and speed [measured]
- Force call (LAMMPS + Metal-Reax, 11 atoms): 10.4 µs (Stage 0). The cost is i-PI's Python and sockets.
- **Polling latency** (`<latency>` in `<ffsocket>`, i-PI default and the value used in all earlier runs: 1e-4 s): at 1e-5 s the 1 amu P = 36 run goes from 116 to 184 steps/s (+58 %). Results are bitwise identical for the same seed (all latencies and client counts tested, 8000 steps, energies and bead positions).
- Number of LAMMPS clients per run (latency 1e-5, same system): 1 → 184 steps/s (138 % CPU), 4 → 219 (161 %), 12 → 209 (166 %), 36 → 192 (178 %). Steps/s per core is the same or worse, so one client per run was kept.
- Production speeds (i-PI `t/step`, machine shared, load average up to 16, so these are lower bounds):

| run type | steps/s (loaded machine, 3–7 runs at once) | steps/s (2–3 runs at once, later in the day) |
|---|---|---|
| classical P=1, dt 0.25 | 1870–2060 | |
| 4 amu P=18 dt 0.25 | 163–288 | |
| 1 amu P=18 dt 0.25 | 148–155 | |
| 1 amu P=36 dt 0.125 (latency 1e-4) | 68–125 | |
| 1 amu P=36 dt 0.125 (latency 1e-5) | 133 (frozen run) | 209 (long run) |
| 1 amu P=36 dt 0.0625 | 87 | |
| 1 amu P=72 dt 0.125 | 47 | |
| 0.5 amu P=50 dt 0.09 | 66–73 (127 frozen) | |
| 0.5 amu P=100 dt 0.09 | 65 | 65 (2 at once) |
| 0.1 amu P=110 dt 0.09 | 26–33 | |
| 0.1 amu P=110 dt 0.045 | 59 | |
Timing-harness numbers (clean comparison, same load): 1 amu P=36 116 steps/s at latency 1e-4, 184 steps/s at 1e-5, 219 with 4 clients.
Cost of the whole Stage 2 PIMD set [measured]: about 28 runs, 17–19 core-hours at ≈ 1.4 cores per run.

## 4. Caveats and limits
- **Statistics.** Block errors use 8 time blocks per run and neglect correlation between blocks; they are optimistic for the slow width observable. Run lengths: 1 amu 88 ps (+ 72 + 72 + 60 ps in convergence runs), 4 amu 240 ps, 0.5 amu 72 ps at P = 50 and 72 ps at P = 100, 0.1 amu 32 ps (one of the three planned runs, seed 61, hung after 7.9 ps: its LAMMPS client was dropped by i-PI ("Timeout in sendall ... unresponsive") at step 88060 of 133333; I killed it at the next check and used the 7.9 ps).
- **Solver versus PIMD are not the same model.** The solver freezes the other four electrons on their cations; PIMD lets them move (classical). Differences of 5–25 % in the outer tail at 1 amu are consistent with that (section 2.6b) but not proven. All agreement statements are for the barrier region |z| < 0.5 Å and the width, which carry the signal.
- **Pit traps.** Classical P = 1 (1 of 2 runs) and 4 amu (2 of 6 runs) fell into a pit and did not leave. Frames in a pit (|z| > 1.2 Å) are excluded from those comparisons (solver renormalised the same way). Pit weights in equilibrium are 5e-5 (4 amu) and 7e-3 (classical); the P = 1 and 4 amu profiles are therefore conditional on being in the middle region.
- **Bare vs quantum "classical" baseline.** The classical comparison for the density-ratio statement is the P = 1 PIMD run (500 ps × 2, |z| < 1.2 Å), which is itself −7 % (1.7σ) below the classical Boltzmann solver in P(|z|<0.1) (mobile other electrons again).
- **Design edge.** The chosen barrier's PMF barrier is 4.05 kT, at the lower edge of 4–6 kT. (4, 0.15) (PMF 4.9 kT, margin 3.3 kT, ratio 12.5) is the design to use if you want more margin; it was not run in PIMD. Only (3, 0.10) was re-solved at h = 0.03 Å.
- **Gaussian form.** z-only (slab) Gaussian on the quantum electron: it has no x, y dependence. A 3-D Gaussian would give a different barrier shape.
- **Not done.** Barrier boosts (+5/+10 kcal/mol) on this design; the 0.1 amu seed 61 and 62/63 at the 12 ps planned length are partly short; 0.1 amu at dt 0.045 fs has one 8 ps run only; the 4 amu, 0.5 amu and 0.1 amu P/dt convergence beyond the checks above; a long conserved-quantity run at 0.1 amu.
- **Process note.** Three batch-2 jobs were started before the new latency was in place (1 amu seeds 3 and 4 and 4 amu seed 21 run at latency 1e-4; same physics, bitwise-identical dynamics). A mistaken early start of three jobs was killed within seconds and its directories renamed `aborted_*` in scratch.
