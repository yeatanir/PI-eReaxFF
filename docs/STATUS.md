# Status log

Newest first. Plan: ../../../PI_eReaxFF_BUILD_PLAN.md

## 2026-11-12: switch to Ag/El
- Owner: use the published Ag/El force field (Evangelisti et al. 2020) instead of graphene. ffield built from the SI (inputs/ag_el), cluster generator make_ag_case.py. ereac NVE at 300 K, dt 0.25: stable, Etot span 0.02-0.14 kcal/mol over 2-5 ps for Ag3..Ag7. Exact frames in reference_data/phase4_frames_ag.
- Control-key correction: ereac reads imdmet/mdtemp/iravel, not nmethod/tset; earlier graphene NVE runs started from rest (see phase0 README).
- Fork F (torsion) done: residual = truncated pi in the legacy cell construction when atoms cross the box edge (not a Metal-Reax defect). Fork E (LAMMPS adapter) done for C/El: all 51 graphene frames pass (worst energy 1.3e-5 kcal/mol); Av on allow-list; Ag/El validation requested.

## 2026-11-11 (very late): classical NVE in ereac, realistic-regime frames
- Done by me: relaxed periodic 60 C/60 El in ereac, NVE at dt 0.1/0.25 fs, 51 exact frames in reference_data/phase4_frames. Energy conserved to 0.01 kcal/mol for 200 fs, then the flat sheet buckles/collapses (T ~10 000 K): force-field/electron-model instability, not the integrator. With the ffield as written (taper 0) it collapses within 100 fs. **ffield/taper must be chosen before stability studies.** Details in validation/phase0_reference/README.md.
- Running in parallel: fork E (LAMMPS adapter, worktree forks/metal-reax-lammps), fork F (torsion baseline, worktree forks/metal-reax-torsion), fork D (i-PI loose ends).

## CPU policy (owner request, 2026-11-11 late)
Machine: 18 cores (Apple M5 Max). Measured peak load 24 while the JAX probe (all cores) and forks ran together = oversubscribed. Policy: at most 8 cores in total across everything I run (Metal-Reax fork: -j4 / ctest -j2 / OMP 2; i-PI fork: at most 3 heavy processes, OMP_NUM_THREADS=1; JAX runs must set OMP_NUM_THREADS=4 and limit XLA threads). The legacy ereac is single-threaded. High load also makes the rebuilt ereac fail to start (dyld error); run it when the machine is quiet.

## Owner directions (2026-11-11, late)
- Metal-Reax tasks 3.3/3.4: term by term with finite differences (fork C resumed).
- V2 periodicity difference (60 C/60 El: V2 near periodic ereac, 0.213 kcal/mol off; 112 C Jamil case: V2 non-periodic): **accepted and documented, to be revisited later** (see validation/phase0_reference/README.md, task 0.2b).
- Stability question: **wait for the Metal-Reax engine**; the JAX probe (P=4, NVE, dt 0.1 fs) is left to finish and is documented, no further JAX stability runs.
- Next items picked: i-PI loose ends (fork D running), LAMMPS adapter for eReaxFF, performance work (both after the engine has sharing and over/under terms).

## 2026-11-11 (workstream B: in-house F1 fix)
- `forks/inhouse-engine` (copy of the old engine; original untouched): production-path test_36 (coupled classical-atom + ring-polymer harmonic potential through the mixed driver) **fails before** the F1 fix (<X^2> off by +47/+126/+288 % at P = 4/8/16) and **passes after** (all within 0.32 %, <= 2 SE; 6 seeds × 2e6 steps). One-line fix in `half_kick_from_summary`. Commits 418e6cf, 7db9bc4, a80d240, 1e0bdb1. See `../validation/phase1_inhouse_f1/README.md`.
- Whole-system `conserved_like` implemented; in NVE it exposes the audit's **F3** defect (carbon heats; drift ~dt^1.2). A scratch experiment (electron first kick reuses start-of-step forces) removes it (500× smaller at dt 0.5 fs). Patch saved, **not committed**; decision needed whether to adopt it.

## 2026-11-11 (Workstream C, tasks 3.3 / 3.4 done)
- Electron sharing, lone pair (stlpc, C2 with legacy constants) and over/under-coordination (pxel1/pxel2 factor) implemented in the CPU-64 engine; FD tests pass (<= 2.6e-6 kcal/mol/A on gradients up to 800); ten purpose-made sharing cases and the 225-atom / 60C60El cases (non-periodic and periodic) agree with the rebuilt legacy executable: lone pair <= 1e-10, over+under 3e-9 abs, vdW/Coulomb/electron 1e-14..1e-16 relative. CTest 38/38.
- Legacy single-precision / truncated constants flagged (c1c, fothird, 0.040/0.16, pi). Open: torsion differs 3e-6 relative in small clusters (stock vs legacy, independent of electrons). See `../validation/phase3_engine/README.md`.

## 2026-11-11 (Workstream C, Metal-Reax port)
- Fork cloned to `forks/metal-reax`, branch `ereaxff`. Baseline 18/18. Task 3.1 (parser/data model) and 3.2 (electron pair term) done and committed (`5dae7b0`, `c8cb917`); 22/22 tests with the optional real-file tests.
- Electron energy of the 225-atom Jamil geometry (non-periodic, taper 10 A) = 16563.553412899309 kcal/mol vs exact 16563.553412899320 (6.6e-16 relative). Finding: legacy c1c is a single-precision literal (332.0638122558594).
- Not started: 3.3 (electron sharing, lone pair), 3.4 (over/under-coordination); ordered to-do list in `../validation/phase3_engine/README.md`.

## 2026-11-11 (night)
- i-PI branch: L-P1, L-P2, guard, multi-seed test all pass; committed (`1a47729`). ⟨X²⟩ deficit gone (+1.0 σ with 4 seeds, stiffer trap).
- Task 0.2 done: export patch + runner with retry + strict comparison; rebuilt `ereac` equals V2 to machine precision (energies 2e-13 kcal/mol, forces 8e-15 relative). Audit finding F8 (energy/force conversion mismatch) was wrong: V2's force constant is 0.043364 × 9.6485e-3.
- Phase 2: JAX driver for i-PI written and validated against exact `ereac` (force diff 6.3e-7 eV/Å); F2 probe running (P=4, NVE, dt 0.1 fs).
- Owner instruction: use the ffield that is already there for now (parameter 38 = 0 as written; JAX effectively uses a 10 Å electron taper). Authoritative ffield to be chosen later.

## 2026-11-11 (evening)
- **Task 0.3 done.** NaN cause = CRYSTX field widths (found by an independent review, verified). Rebuilt `ereac`, non-periodic with electron taper 10 Å, reproduces all 14 saved V2 energy terms (max diff 0.0035 kcal/mol). Periodic runs differ strongly. See `../validation/phase0_reference/README.md`.
- i-PI branch committed (`1b8f05e`). Added since: guard that projected modes are zero (uncommitted), L-P1 passes (byte-identical to upstream when `<classical_atoms>` is absent), negative control run (booking the kick momentum drifts 0.031 kT per 8000 steps vs 1e-6 kT), multi-seed redesign of the <X^2> test (running).

## 2026-11-11 (later)
- Phase 1.1/1.2 tests pass: see `../validation/phase1_ipi_mixed/README.md` (conservation in NVE and NVT, ring-polymer variances within 1 σ, classical-atom ⟨X²⟩ 2.2 σ low, still to be tightened).
- Phase 0.1: stock source does not build with stock `cbka.blk` (finding P0-1, see `../validation/phase0_reference/README.md`). Full-size build succeeded in scratch after a one-line declaration fix.
- Owner asked to build in scratch and copy build files to the original folder: done (`Margaret/eReaxFF_src/build_2026-11-11/`).
- Rebuilt `ereac` runs, but gives `Energy not a number` on the Jamil geometry (not the taper); open item for task 0.3.
- Nothing committed in `forks/i-pi` yet.

## 2026-11-11
- Layout created (see ../../README.md).
- Phase 1.1: `classical_atoms` patch written on `forks/i-pi` branch `classical-atoms`
  (engine/motion/dynamics.py, inputs/motion/dynamics.py, engine/properties.py), plus test PES
  `ipi/pes/coupled_harmonic.py`. Tests in progress (see validation/phase1_ipi_mixed).
- 2026-11-11 torsion baseline (fork F): cause = legacy dista2 with truncated pi at atoms crossing the cell wrap; exactcell references make all 14 reference cases pass with torsion enabled (<=2.5e-10); Ag comparison not done (no Ag ffield on disk). See validation/phase3_engine/TORSION.md
- 2026-11-11: LAMMPS adapter (E): C/El/Av validated, Ag/El 150 frames pass, fix ipi P=1 matches fix nve; see validation/phase3_engine/ADAPTER.md. Ag reference geos have END before MOLCHARGE (charges ignored).

## 2026-11-12 (workstream H: tunneling prerequisites A4/A5/A6 numbers/A7)
- A4 done: `ipi/pes/gaussian_barrier.py` (i-PI fork `f3b9e82`), unit tests pass; combination with the eReaxFF component via a 2nd ffsocket: V0 = 0 byte-identical to no barrier; barrier on the electron's start region raises the potential by the analytic 0.065 eV and the per-bead conserved span is unchanged (0.0887 vs 0.0885 kcal/mol, 100 steps).
- A5 done (exact ereac): symmetric Ag3...Ag3 pair (5 El) has degenerate sites; narrow pits (k ~ 45 eV/A^2) on the apex Ag, intrinsic barrier 9-12 kcal/mol only 0.3-0.4 A from each pit, soft middle region; cluster + Ag+ rejected. Recommended gap 4.0-4.5 A. A6 numbers: pit needs P >= 35/49/109 and dt <= 0.30/0.22/0.096 fs for 1/0.5/0.1 amu. A6 histogram convergence not done.
- A7 done (analysis script + tests). Notes: validation/tunneling_ag/A4_A7.md.
- Observation: an external SOFTEXIT 'kill signal' ended one of my i-PI runs mid-way (another workstream's cleanup is the likely cause); avoid broad `pkill -f i-pi`.

## 2026-11-12: workstream G (tunneling prerequisites A1-A3) - see validation/tunneling_ag/A1_A3.md
- A1 PASS: Ag/El P=8, 8 LAMMPS clients (each own dir); 10 bead geometries: client vs offline plugin 2e-13 rel, plugin vs exact ereac (same f10.5-rounded positions) 4.9e-11 rel. `fix ipi` triples the LAMMPS force work per request (2 of 3 steps repeat x).
- A2 PASS: El mass (i-PI 1/0.5/0.1; LAMMPS data file 1/0.5/0.1/25) leaves energy/forces bitwise identical; flags come from ffield mass only (ffield_parser.cpp:237).
- A3: printed conserved quantity is per bead (properties.py:1573). 15 ps PILE-L equilibration (R_g ~stationary: P8 last blocks 0.1505/0.1502 A, P32 0.161-0.168 A) then 2 ps NVE. Span/bead P1 0.081, P8 0.158, P32 0.297 kcal/mol at dt 0.25 (P32 = 3.66x P1: above the factor-3 flag); at dt 0.125: 0.063/0.069/0.151 (P32 2.4x), at 0.0625: 0.059/0.063/0.122 (2.1x). Fast ripple falls steeply with dt; a dt-independent slow wander (0.10 vs 0.045 kcal/mol) remains and is NOT explained. Use dt <= 0.125 fs at P >= 32.
- Gotchas: forks/metal-reax-lammps/build/ plugin is stale (use build-ns); i-PI zero-pads bead file names for P > 10; the first 3 ps P=8 equilibration was not stationary and was redone.
- Scripts: project/inputs/tunneling_ag/{run_pimd_lammps,a1_compare,a2_mass,a3_pipeline}.py. Commits: project 391f3ab (+ earlier); no source changes in forks/*.

## 2026-11-12 (coordinator request, Option 2): revised Part B written
- Gap 4.0 A landscape measured with exact ereac: pits degenerate (mirror), middle region +0.14 kcal/mol above the pits, intrinsic barrier tops at z = +-1.70 A, 10.51 kcal/mol (17.6 kT), k_top = -9.44 eV/A^2 -> T_c = 183 / 367 / 519 / 1160 K for 4 / 1 / 0.5 / 0.1 amu (local-curvature estimate). Fixed-others scan is fair (relaxing the other 4 electrons changes the barrier by +0.09 kcal/mol); no lower off-axis path within +-0.6 A. Protocol, run table, sampling check and cost (t_call is a guess: please supply the measured value) in validation/tunneling_ag/PART_B_REVISED.md. The full run table is ~3000 h wall at 20 ms/call on 2 cores: staged plan proposed.
