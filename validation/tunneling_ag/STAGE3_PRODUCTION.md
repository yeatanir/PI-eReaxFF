# Stage 3 — production set: Mac Mini check, feasibility estimate, status

Date: 2026-11-11. Labels: **[measured]**, **[derived]** (arithmetic from measured numbers), **[estimated]** (scaled, not run).
Scripts: `inputs/tunneling_ag/stage3_jobs.py` (job lists), `run_barrier_lammps.py` (writes `run_info.txt` per run: command, host, seed, commits, `CONTROL_USED NULL(legacy)`, design; reduces bead trajectories to `ring.npy`; adds the atom-resolved centroid-virial kinetic energy of the quantum electron as sim.out column 8), `scratch/tun_stage3/pool.py`.

## 1. Mac Mini [measured]
- Machine: arm64, 10 cores, `/opt/homebrew/bin/python3` 3.13 with numpy 2.3 and scipy (system `/usr/bin/python3` 3.9 has no numpy). No `lmp`, gfortran or cmake in PATH, but none is needed: `lmp` and the Metal-Reax plugin link only system libraries; `ereac` finds Homebrew's libgfortran.
- Copied (rsync, 380 MB) to `~/eb_stage3`: project, forks/i-pi, forks/metal-reax-lammps (tests + `build-ns/plugin`), external/install_lammps, `ereac_export_build`, `ipi_reg_diag`, `ipi_upstream`. **Same commits:** project `42cb693` (MacBook `42cb693`, clean), i-PI `071f6b9`, metal-reax-lammps `120b609`; `lmp` md5 `22ce4b09…` and plugin md5 `966f5f0f…` identical on both machines (`~/eb_stage3/COMMITS.txt`; the Mini copy has no git for the worktree fork, so the identity is by md5).
- **Cross-machine run, same seed (1 amu, P = 72, dt 0.125 fs, seed 11, barrier on, latency 1e-5): bitwise identical.** First 3 steps with force dumps (289 frames, 11 atoms): positions and forces max abs difference **0.0** (relative 0.0; first-frame |f|max = 24.9 kcal/mol/Å). 2000 steps: all 8 columns of `sim.out` (2001 rows) and the bead positions of the quantum electron (`ring.npy`, 72 × 21 × 3) identical (max difference 0.0). Required: first-step forces within 1e-10 relative → met with margin (exactly equal).
- `tools/run_all_checks.py` on the Mini: first run 9 PASS, 1 SKIP (`eparams`, needs `--mr-build`), 2 FAIL for **portability reasons, not physics**: (i) `tools/discontinuity/run_nve.py` had a hard-coded MacBook path to the LAMMPS harness (fixed in commit `42cb693`: now relative to the tree); (ii) `ipi_regression`: the Mini's pytest 9.2.0 names parametrised tests `regtest0…` instead of the path, so `-k CLASSICAL` selects nothing (fixed by installing pytest 9.0.3, the MacBook version). The re-run after both fixes: **10 of 11 checks PASS** (frames_ag, frames_graphene, fd_forces, smooth_vs_ereac, conservation_p1, classical_atoms, ipi_lammps_short, schrodinger, periodic_ag + the SKIP of `eparams`), numbers identical to the MacBook (e.g. frames_ag worst |dE| 1.6e-11 / 8.7e-9 kcal/mol, conservation ratio 3.95). **`ipi_regression` does NOT pass on the Mini: the i-PI regression runner hangs (pytest idle at 0 % CPU for > 4 min; 3 s on the MacBook, Python 3.14; the Mini has Python 3.13 and pytest 9.0.3).** This is the same runner hang seen earlier on the MacBook before the `ipi_reg_diag` fix; I did not diagnose it on the Mini. So the requirement 'all checks pass on the Mini' is **not met**; the physics checks and the bitwise cross-machine run are.
- Speed on the Mini [measured]: 1 amu P = 72: 2000 steps in 18.1 s (MacBook 20.2 s); **0.1 amu P = 220, dt 0.045 fs: 28.5 ms/step = 35.1 steps/s (MacBook 35.0)**. Same speed per run as the MacBook when the Mini is quiet.
- Mini state during the work: load 1–9 (varies; other users' jobs), **disk 3.2 GiB free and falling (4.4 → 3.2 GiB within one hour; my copy is 0.4 GiB)**. This is the real limit for the Mini.

## 2. Speeds used for the estimate [measured on the MacBook, latency 1e-5, load ≈ 3.6]
| run | steps/s | CPU (all processes of the run) |
|---|---|---|
| 1 amu P = 72, dt 0.125 | **109.5** | 133 % |
| 0.1 amu P = 220, dt 0.045 | **35.0** | 133 % |
| (earlier, other loads) P = 36: 184–209; P = 100: 65–79 | | |

## 2b. Important correction to the speeds: this MacBook has 6 performance + 12 efficiency cores [measured: `sysctl hw.perflevel0/1.logicalcpu` = 6 / 12]
The isolated timings above (109.5 and 35.0 steps/s) were taken with the machine nearly idle (on performance cores). With **5 production runs at once** (CPU usage 52 %, of which 27.6 % system time from the polling threads) the same 1 amu P = 72 run does **19.1 ms/step = 52 steps/s, i.e. 2.1× slower** (P = 144: 37.4 ms/step = 27 steps/s). Five concurrent runs deliver 262 steps/s in total, about 2.4 isolated runs. A 'cap of 12 cores' therefore corresponds to roughly 6 performance-core equivalents, and the wall times below, which use isolated speeds, are optimistic by a factor of about 2 whenever more than 2–3 runs are active. The measured concurrent speed is used in the wall estimates below as 'concurrent'. (The Mini's cores were not checked for this split.)

## 3. Feasibility estimate [derived from the above; P = 144/110/100/200 and 4 amu speeds are scaled ∝ 1/P, [estimated]]
Per mass: 8 mobile + 4 frozen runs of 50 ps, one bead check (2P; P/2 at 0.1 amu), one dt/2 run of 50 ps, one NVE run of 10 ps = 15 runs. Core-hours use 1.33 cores per run.

| mass | P, dt (fs) | steps per 50 ps | steps/s | wall per 50 ps run | run-hours for the set | core-hours |
|---|---|---|---|---|---|---|
| 4 | 36, 0.125 | 400 000 | ≈ 190 [est.] | 0.58 h | 12 × 0.58 + 1.1 (P = 72) + 1.2 (dt/2) + 0.1 = **9.4** | 12.5 |
| **1** | 72, 0.125 | 400 000 | **109.5** | 1.01 h | 12 × 1.01 + 2.0 (P = 144) + 2.0 (dt/2) + 0.2 = **16.4** | **21.8** |
| 0.5 | 100, 0.09 | 555 556 | ≈ 79 [est.] | 1.95 h | 12 × 1.95 + 4.0 (P = 200) + 3.9 (dt/2) + 0.4 = **31.7** | 42 |
| 0.1 | 220, 0.045 | 1 111 111 | **35.0** | **8.8 h** | 12 × 8.8 + 4.4 (P/2 = 110) + 17.6 (dt/2) + 1.8 = **129.4** | **172** |
| classical | 1, 0.125 | 4 000 000 (500 ps) | ≈ 1900 | 0.6 h | 6 × 0.6 = 3.5 | 3.5 |
| 1 amu, design 2 (V0 4, σ 0.15) | 72, 0.125 | | | | 16.4 | 21.8 |
Total without 0.1 amu: ≈ 85 run-hours (isolated speed) ≈ 105 core-hours. **At the measured concurrent speed (÷ 2.1) that is ≈ 180 run-hours; with 5 concurrent runs ≈ 36 h wall on the MacBook**, of which 1 amu design 1 + classical ≈ 7 h and design 2 another ≈ 7 h. With 0.1 amu: +129 run-hours isolated.

**0.1 amu at 220 beads × 50 ps × 14 runs is not practical as specified, and I would not do it as a full table.** Numbers: 8.8 h per run; 14 runs = 129 run-hours = 172 core-hours. On the Mini, at most ~5 concurrent runs (load from other users, 10 cores) → ≈ 26 h wall if the speed stays at 35 steps/s and the Mini stays free; and **the disk cannot hold it as the run script writes today**: bead xyz output for 220 beads at 5 fs stride is ≈ 1.4 GB per run (220 files × 10 000 frames × 650 B) before reduction → 5 concurrent runs need ≈ 7 GB against 3.2 GiB free. Thinning to 25 fs stride (≈ 0.3 GB per run) or chaining 10 ps segments is required before launching anything there.

**Smallest cut that keeps the convergence checks (proposal):** 0.1 amu on the Mini: 4 mobile runs (L, R, L, R) + 2 frozen at 50 ps; bead check P/2 = 110 at 50 ps (cheap, 4.4 h); **dt/2 at 25 ps only** (8.8 h); NVE 5 ps; bead-position stride 25 fs; 4 concurrent runs. = 6 × 8.8 + 4.4 + 8.8 + 0.9 = **66 run-hours (88 core-hours), ≈ 17 h wall at 4 concurrent**, about 0.7 GB of live bead output per run (2.8 GB at four concurrent: still too close to the free disk; segmenting into 10 ps pieces reduces it to ≈ 0.15 GB per run). Statistical power: 6 × 50 ps = 300 ps at 0.1 amu, against 32 ps in Stage 2. If you want the smaller disk risk, the bead check and dt/2 can run on the MacBook (cap permitting) while the Mini does the six production runs.

## 4. Status
- **Launched at 22:03 (MacBook, 5 concurrent, pool in `scratch/tun_stage3`, jobs in `jobs_m1.txt`):** 1 amu, design 1 (V0 3 kcal/mol, σ 0.10 Å), P = 72, dt 0.125 fs, 50 ps: dt/2 run, 2P (P = 144) run, 8 mobile-others runs (L/R alternating, seeds 2001–2008), 4 frozen-others runs (seeds 2011–2014), NVE 10 ps from the end of mobile run 2001; then classical P = 1 (4 mobile + 2 frozen, 500 ps); then the same 1 amu set for design 2 (V0 4, σ 0.15; seeds 3001…). Every run: PILE-L 300 K, τ 100 fs, latency 1e-5, legacy control (`pair_style reaxff/metal NULL ereaxff yes electron_taper 5.7`), Ag fixed; `run_info.txt` + RESTART kept per run. MacBook load was 3–4 before launch; my 5 runs add ≈ 6.7 cores.
- **Not started:** the other masses; the h = 0.02 Å final solver re-solve (larger box, states to 15 kT, thermal weight below W(0), kinetic energy per mass); the analysis script for the Stage 3 observables (centroid free-energy profile, radius of gyration vs centroid z, straddling fraction, centroid-virial KE vs solver, entry times, run-to-run spread errors). The pit rule (exclude frames with centroid |z| > 1.2 Å; report entry times) and the 2 ps discard are fixed for the analysis.
- **Waiting for your decision:** the 0.1 amu cut (section 3), and whether 4 amu and 0.5 amu should run at the settings you gave while 1 amu is running (they fit the MacBook cap only after 1 amu finishes: ≈ 11 more h).

## 5. Agreement check of the first finished 1 amu runs (2026-11-12) [measured]
`inputs/tunneling_ag/analyze_stage3.py` (ring.npy; first 2 ps discarded; pit rule |centroid z| > 1.2 Å; errors = spread between independent runs). 1 amu, design 1, P = 72, dt 0.125 fs, 50 ps each. Reference: Stage 2 P = 36 pooled (4 runs, 88 ps) and the h = 0.03 Å solver.

| observable | Stage 3, 3 runs (s2001–2003) | Stage 3, 6 runs (s2001–2006) | Stage 2 (P = 36) | solver |
|---|---|---|---|---|
| rms z (Å) | 0.3643 ± 0.0035 | **0.3663 ± 0.0018** | 0.3678 ± 0.0026 | 0.3649 |
| P(\|z\|<0.1) | 0.0878 ± 0.0013 | **0.0892 ± 0.0012** | 0.0942 ± 0.0023 | 0.0921 |
| P(\|z\|>0.5) | 0.1454 ± 0.0058 | 0.1497 ± 0.0033 | 0.1548 ± 0.0052 | 0.1482 |
| P(\|z\|>0.75) | 0.0122 ± 0.0025 | 0.0132 ± 0.0012 | 0.0142 ± 0.0014 | 0.0130 |
| ρ(0)/ρmax | 0.317 ± 0.016 | **0.322 ± 0.009** | 0.374 ± 0.017 | 0.365 |
| fraction of rings straddling z = 0 | 0.230 | 0.232 ± 0.002 | – | – |
| mean radius of gyration (Å) | 0.1882 | 0.1882 ± 0.0001 | – | – |

Width, outer density and tail agree with Stage 2 and the solver (rms −0.4 % / +0.4 %, P(|z|>0.5) −3 % / +1 %, tail −7 % / +2 %, all < 1σ). **The barrier-region density is lower: P(|z|<0.1) −5.3 % against Stage 2 (2.0σ) and −3.2 % against the solver; ρ(0)/ρmax −14 % (2.7σ) / −12 %.** Per-run ρ(0)/ρmax ranges 0.30–0.35 (±7 % per run), P(|z|<0.1) 0.085–0.093. The two single check runs sit closer to the references: P = 144 (2P) gives P(|z|<0.1) 0.0907, ρ(0)/ρmax 0.358; dt/2 gives 0.0953 and 0.365 (one run each, ±5 %, not a convergence test). Not yet explained and to be re-checked with all 8 + 4 runs: the Stage 2 P = 72 and P = 18 runs (3 and 4 runs) gave 0.0948 and 0.0945, so a 3–7 % difference in the barrier-region density between Stage 2 and Stage 3 at the same P and dt is either run-to-run scatter underestimated by 3–6 runs or a small systematic. No pit trap in any run; the centroid left |z| ≤ 1.2 Å briefly (a few frames) at 16 ps in s2005 and s2006 and at 47 ps in the dt/2 run. Conserved quantity (s2001–2003, 50 ps): std 0.051–0.058 kT, span 0.24–0.25 kT, drift −0.003…+0.002 kT/ns.

## 6. Queue change (2026-11-12)
Finished runs of 1 amu design 1 (other than s2001, needed for the NVE restart) were moved to `scratch/tun_stage3/done/` (outputs are already thin: bead xyz are reduced to `ring.npy`; 92 MB for the moved runs). A new gated pool (`pool2.py`, never more than 5 runs of mine at once) works through `jobs_order4.txt`: 0.5 amu (P = 100, dt 0.09 fs; design 1) first, then 4 amu (P = 36, dt 0.125), then the classical runs and 1 amu design 2. The 1 amu NVE run and the 0.5 amu dt/2 run were started by a pool that did not gate correctly for about a minute (6 of my runs at once for a short time; load 13.8), then stopped. 0.1 amu stays unstarted until the Mini disk question is settled. Polling-delay timings skipped as instructed.
