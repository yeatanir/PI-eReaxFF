# PI-eReaxFF

Path-integral molecular dynamics (PIMD / RPMD) for **eReaxFF**, the ReaxFF variant with explicit pseudo-electrons
(atom type `El`), built from three pieces:

1. an **eReaxFF mode of the Metal-Reax LAMMPS pair style** (`pair_style reaxff/metal ... ereaxff yes`), validated term by term against the legacy `ereac` Fortran code;
2. a **patch to i-PI** (`<classical_atoms>`) that lets classical atoms (ions, fixed-charge cores) live inside a ring-polymer run next to quantum electrons;
3. the **drivers, run scripts, solvers and validation reports** that tie them together, plus a first scientific demonstration (barrier penetration of one pseudo-electron in an Ag donor–barrier–acceptor system, compared against a 3-D Schrödinger solver).

> **Status in one paragraph.** For the **Ag/El system with fixed charges** (the published Evangelisti–Fichthorn–van Duin force field) the stack works end to end and is validated (sections 4 and 5). It is **not** a general eReaxFF engine yet: fixed charges only, `C`/`El`/`Ag(Av)` only, graphene dynamics collapse after about 110 fs (undiagnosed), the i-PI patch is not upstreamed, and one test (`lammps_charge_models`) was 63/64 at the last run. Section 6 lists everything known to be missing or wrong. This repository is **private**; the license is undecided and some content is legacy-derived (section 9).

---

## 1. Architecture

```
                        i-PI  (patched: <classical_atoms>)
        ring-polymer propagation, PILE-L / PILE-G thermostat, bead bookkeeping
                    |  unix/inet socket (one force request per bead)
        +-----------+------------------------+--------------------------+
        |                                    |                          |
   LAMMPS  `fix ipi`                  Python driver                Python driver
   + Metal-Reax plugin                ipi_ereac_driver.py          ipi_ereaxff_jax_driver.py
   pair_style reaxff/metal            (exact legacy `ereac`,       (JAX kernel, central cell,
   ... ereaxff yes                     slow, reference only)        fixed charge; probe only)
   (CPU-64 engine, 10.4 us/call)
        |
   control file: legacy.ctrl (validation)  or  smooth.ctrl (mobile atoms, Ag/El only)
   eparams file: per-element Z / sharing sign / tested flag
```

Conventions (verified in the i-PI source, see `validation/phase1_ipi_mixed/README.md`):

- i-PI uses ring frequency `omega_P = P k_B T / hbar`, thermostat temperature `P T`, and an unscaled bead potential.
- A **classical atom** is a ring polymer whose internal normal modes `k >= 1` are projected to zero (`<classical_atoms>`). The centroid then obeys `dP/dt = (1/P) sum_j f_j`: physical mass, bead-averaged force, temperature `T`.
- Energy bookkeeping: momentum removed right after a force kick is **not** booked; thermostat noise added and then removed **is** booked in the conserved quantity.

## 2. Repository layout

| Path | Contents |
|---|---|
| `pair_style/` | Metal-Reax with the eReaxFF mode, branch `ereaxff-lammps-merged` of the working fork (own upstream README: `pair_style/UPSTREAM_README.md`). C++ engine (CPU-64), LAMMPS plugin, tests, tools. |
| `ipi/patches/` | `git format-patch` series (0001-0010) for i-PI on upstream commit `93162d2`: `classical_atoms`, tests, a z-axis Gaussian barrier PES, regression layout. |
| `ipi/UPSTREAM_PR_DRAFT.md` | Draft text for an upstream pull request (**not submitted**). |
| `drivers/` | i-PI Python drivers: `ipi_ereac_driver.py` (exact `ereac`), `ipi_ereaxff_jax_driver.py` (JAX). |
| `inputs/ag_el/` | Ag/El force field (`ffield_Ag_El`, built from the SI of the 2020 paper), cluster generator `make_ag_case.py`, `make_ffield_ag_el.py`, launcher `run_ipi_ereac.py`. |
| `inputs/ag_pimd/` | `run_ipi_lammps.py`: i-PI + LAMMPS (`fix ipi`) launcher for Ag/El frames. |
| `inputs/control/` | `legacy.ctrl`, `smooth.ctrl` and their README. |
| `inputs/tunneling_ag/` | Tunneling demonstration: system builder, barrier launcher `run_barrier_lammps.py`, Schrödinger solvers (`schrod1d.py`, `schrod3d.py`, `schrod3d_barrier.py`), analysis scripts. |
| `tools/` | `run_all_checks.py` (one command for the checks), `compare_phase0*.py`, `run_ereac*.py` helpers, frame generators, `discontinuity/` (cutoff-jump hunt). |
| `validation/` | All validation reports (Markdown) and small data/figures, by phase. |
| `reference_data/` | Small reference sets: Ag/El exact frames (`phase4_frames_ag`), stored Schrödinger data, a short i-PI+LAMMPS baseline. |
| `docs/` | Build plan, audit report, status log, progress report, task prompt. |
| `legacy_derived/` | **Restricted**: small patches against the legacy `ereac` source. Read section 9 first. |

## 3. Quick start

The scripts were written for a **build-root layout** and resolve their dependencies relative to it. Create:

```
build_root/
  PI-eReaxFF/                      this repository
  forks/i-pi/                      upstream i-PI with ipi/patches applied
  forks/metal-reax-lammps/         a copy of pair_style/ (build in build-ns/)
  external/install_lammps/bin/lmp  LAMMPS built as a shared library with the `ipi` fix
```

`run_barrier_lammps.py` honours `IPI_ROOT` and `REAX_PLUGIN`; `tools/run_all_checks.py` uses `--mr-build` and `--ipi-reg-dir`. Other scripts assume the layout above. **I have not re-run the scripts from a fresh checkout of this repository; they were run in the original tree.** Expect to fix paths.

**3.1 Patch i-PI**

```
git clone https://github.com/i-pi/i-pi.git forks/i-pi
cd forks/i-pi && git checkout 93162d2 && git am ../../PI-eReaxFF/ipi/patches/*.patch
IPI_UPSTREAM_ROOT=<unpatched i-pi checkout> python3 -m pytest -q -s ipi_tests/classical_atoms/test_classical_atoms.py   # about 10 min
```

The series applies cleanly on `93162d2` (checked with `git apply`; only trailing-whitespace warnings in the stored reference files). Input: `<classical_atoms>` inside `<dynamics>` (atom indices); allowed with nve/nvt and the none/Langevin/PILE-L/PILE-G/SVR thermostats; refused with an explicit error for other ensembles, GLE thermostats, bosons, open paths, `fixatoms` overlap, and contraction.

**3.2 Build the pair style** (from `pair_style/UPSTREAM_README.md`, macOS arm64)

```
tools/fetch_lammps.sh --full /path/lammps && tools/build_lammps_reference.sh /path/lammps /path/build /path/install
cmake -S . -B build-ns -DREAXMETAL_ENABLE_METAL=OFF -DREAXMETAL_BUILD_LAMMPS_PLUGIN=ON \
      -DREAXMETAL_LAMMPS_SOURCE_DIR=/path/lammps/src -DREAXMETAL_LAMMPS_PREFIX=/path/install \
      -DREAXMETAL_FFIELD_DIR=/path/lammps/potentials
cmake --build build-ns -j && ctest --test-dir build-ns -j4 --output-on-failure
```

`-DREAXMETAL_ENABLE_METAL=OFF` selects the CPU-64 engine (the only one used in eReaxFF mode; the Metal GPU backend is not ported to eReaxFF). The working tree used `build-ns`; an older `build/` there is stale. LAMMPS usage:

```
pair_style reaxff/metal inputs/control/legacy.ctrl ereaxff yes [eparams <file>]
```

**3.3 Run an Ag/El PIMD example**

```
python3 inputs/ag_pimd/run_ipi_lammps.py --case <case dir> --lmp <lmp> --plugin <plugin .so> --out out --beads 8 --steps 400
python3 inputs/ag_el/make_ag_case.py <si_text> Ag7_1 case_dir      # builds geo/ffield/control for a neutral Ag cluster
```

**3.4 Run every check**

```
python3 tools/run_all_checks.py                  # quick set, 11 checks, PASS/FAIL/SKIP, exit 1 on any FAIL
python3 tools/run_all_checks.py --only fd_forces schrodinger
```

At the last run all 11 passed in 162 s (`validation/phase3_engine/ROBUSTNESS_2026-11-11.md`). Several checks need reference data that is **not** in this repository (see section 9); those will fail or skip here.

## 4. What is validated (labels: measured / derived)

| Check | Result | Where |
|---|---|---|
| Legacy `ereac` rebuilt vs saved V2 energy terms (non-periodic, taper 10 Å, 225 atoms) | all 14 terms equal to machine precision; largest difference 2.3e-13 kcal/mol [measured] | `validation/phase0_reference/` |
| Metal-Reax eReaxFF electron terms vs `ereac` | agree to 1e-10 … 1e-16 [measured] | `validation/phase3_engine/README.md` |
| Finite-difference tests (sharing, lone pair, over/under) | gradient error <= 2.6e-6 kcal/mol/Å [measured] | same |
| LAMMPS adapter vs `ereac` | graphene 51 frames (worst energy 1.3e-5 kcal/mol); Ag/El frames pass; 30 Ag/El frames adapter-vs-ereac max dE 6.5e-4 kcal/mol, max dF 1.1e-2, both control files (floor = 1e-5 Å geo quantisation) [measured] | `ADAPTER.md`, `SMOOTH_ADOPTION.md` |
| i-PI `classical_atoms`, P=8, 300 K | classical-atom beads identical (<1e-9); ring variance 0.08586±0.00072 vs exact 0.08591; temperature 301.2±0.7 K; NVE conserved span 7.3e-4 kT, drift 7e-8 kT; patch inactive = byte-identical to unpatched i-PI over 400 steps [measured] | `validation/phase1_ipi_mixed/README.md` |
| Restart, baoab, MTS, propagators | pass; restart agrees to 5e-9 relative (unpatched i-PI shows the same) [measured] | same |
| Cost | one force call 10.4 µs (11 atoms) [measured]; i-PI + sockets cost ≈ 0.21 ms per bead request and dominate | `STAGE0_STAGE1.md` |
| Socket polling latency | `<latency>` 1e-4 → 1e-5 s: +58 % steps/s (116 → 184 at 1 amu, P=36), output bitwise identical [measured] | `STAGE2_DESIGN.md` §3 |
| Ag7_1 energy wander | cause found: two cutoff jumps present in legacy `ereac` itself (Ag–Ag `bond_cut` 5.0 Å with BO′ ≈ 0.0021 > `bo_cut`; 3/4-body `thb_cut`); the `smooth.ctrl` setting removes it, span then falls as dt² [measured] | `DISCONTINUITY.md`, `SMOOTH_ADOPTION.md` |

### Stage 2 tunneling demonstration vs the Schrödinger solver

System: Ag donor–barrier–acceptor, gap 4.0 Å, Ag fixed, four classical El (1 amu), one quantum electron of the scanned mass, plus a z-only Gaussian barrier V0 = 3.0 kcal/mol, σ = 0.10 Å at z = 0; PILE-L, 300 K. Solver: 3-D Schrödinger on the frozen eReaxFF landscape (h = 0.03 Å). Errors are block errors over 8 blocks per run (optimistic). From `validation/tunneling_ag/STAGE2_DESIGN.md` [measured]:

| Mass (amu) | P, dt (fs) | rms z PIMD / solver (Å) | P(\|z\|<0.1) PIMD / solver | Verdict |
|---|---|---|---|---|
| 4 (\|z\|<1.2 Å) | 18, 0.25 | 0.3465 ± 0.0128 / 0.3439 | 0.0404 ± 0.0024 / 0.0399 | agree (pit-trapped frames excluded) |
| 1 | 36, 0.125 | 0.3678 ± 0.0026 / 0.3649 | 0.0942 ± 0.0023 / 0.0922 | agree within 1.3σ |
| 0.5 | 50, 0.09 | 0.3877 ± 0.0037 / 0.3767 | 0.1285 ± 0.0014 / 0.1302 | **width/tails fail at ≈3σ** |
| 0.5 | 100, 0.09 | 0.3750 ± 0.0022 / 0.3767 | 0.1279 ± 0.0015 / 0.1302 | agree (ρ(0)/ρmax 3.6σ low) |
| 0.1 | 110, 0.09 / 0.045 | 0.4488 ± 0.0061 / 0.4407 | 0.1592 ± 0.0022 / 0.1618 | agree (≤ 1.6σ); labelled zero-point delocalisation, not tunneling |

Barrier-region density, quantum / classical: PIMD 2.8, 6.5, 8.9, 11.0 at 4, 1, 0.5, 0.1 amu; solver 2.6, 5.9, 8.4, 10.4 [derived, ≈ ±6 %]. Convergence at 1 amu: P → 2P and dt → dt/2 change the barrier-region density by +0.6 % and −4 % (≤ 1.2σ) and the width by 0.5 % and 0.3 %. The conserved quantity at 1 amu, P=36, dt 0.125 fs (60 ps) has std 0.034 kT, drift 0.0015 kT/ns; at 0.1 amu dt 0.09 fs it is noisy (std 0.5–1.2 kT) and dt 0.045 fs is the safer setting. The earlier intrinsic 10.5 kcal/mol pit barrier (`STAGE0_STAGE1.md`) does **not** show barrier penetration: the 3-D solver puts the equilibrium in the broad middle region at every mass.

## 5. Control files and the species table

- **`legacy.ctrl`**: legacy thresholds (identical to `NULL`, verified). Use for **all** validation and reference comparisons against `ereac`.
- **`smooth.ctrl`**: `nbrhood_cutoff 5.2`, `thb_cutoff 1e-9`, `thb_cutoff_sq 1e-18`. Use for any run with **mobile** atoms. It changes the energy surface: on 150 Ag/El frames energy differs by up to 0.44 kcal/mol (rms 0.115), force components by up to 1.62 kcal/mol/Å (rms 0.118), mostly torsion, growing with cluster size. **Ag/El only**: on carbon it shifts the energy by −83 kcal/mol (225 atoms) and the graphene run went unstable. (The planned explicit refusal of `smooth.ctrl` for other elements is a pending follow-up.) `ereac` runs the same thresholds from its control keys `vlbora` and `cutof2`. Every run script logs the control file it used (`CONTROL_USED`); `run_barrier_lammps.py` does not yet have a `--control` option and always uses the legacy defaults.
- **`eparams <file>`** (pair-style keyword): lines `SYMBOL Z SHARE_SIGN tested|untested` override the built-in nuclear-charge and electron-sharing sign tables per element. `C`, `El`, `Av` are accepted; any other element needs a `tested` line; `untested` or absent elements are refused. Example: `pair_style/tests/data/legacy_table.eparams`.
  **Rule for the `tested` flag:** mark an element `tested` only after it passes `ereac` reference frames and finite-difference checks. (This rule is policy, not enforced by code.) The derivative-side sign tables of the legacy `poten.f` have no counterpart in the engine, so an element whose legacy energy and derivative tables disagree would not be detected.
- **Periodic cells**: in eReaxFF mode the pair style refuses cells with an edge < `2*bond_cut + 1` Å (11 Å legacy, 11.4 Å smooth) and periodic triclinic cells. Measured: the adapter agrees with `ereac` at box edges 11–22 Å (dE ≤ 9e-7 kcal/mol) and disagrees at 10 Å (4e-2) and 8 Å (5 kcal/mol).

## 6. Known issues and what is not done

- **Fixed charges only** (`icharg 8`). No charge equilibration with electrons; the charge model for other systems (e.g. Si) is unknown.
- **Elements**: `C`, `El`, `Av` only. Graphene force fields: the authoritative carbon field is unsettled.
- **Graphene dynamics collapse**: a 400 fs periodic 60 C/60 El run conserves energy for about 100 fs and then collapses (kinetic energy 170 → 3500 kcal/mol from ≈ 110 fs). Cause not diagnosed; physical vs integration failure unknown. The energy drift is attributed (unconfirmed) to a stale pair list (`nvlist=25`) plus cutoff discontinuities in `ereac`; the Ag cutoff jumps do not explain it.
- **Source 2 of the Ag7_1 wander**: the 3/4-body `thb_cut` jump is localised to torsion but the responsible bond/quadruple was not identified.
- **`lammps_charge_models`**: 63/64 in `ctest` after the merge; the last test was stopped under load and does not use eReaxFF mode. A 64/64 rerun is pending.
- **V2 periodicity mismatch**: 60 C/60 El differs from the saved V2 oracle by 0.213 kcal/mol (V2 is near periodic, the 112 C Jamil case is non-periodic). Accepted and documented, not resolved.
- **i-PI patch**: not upstreamed (PR drafted only). i-PI's own regression runner never actually hung; the reference output files were missing because `.gitignore` excludes `*.out` (fixed by force-adding them, patch 0010). NPT, GLE thermostats and bosons are refused with classical atoms by design. Quartic-coupling test shows a residual ≈ +1.2 % in bead variance (consistent with step-size error, not resolved beyond that).
- **Process-level**: scripts use paths of the original tree (section 3); `run_all_checks.py` needs reference data not stored here.
- **Tunneling demo gaps**: barrier boosts (+5/+10 kcal/mol) not run; the solver freezes the other four electrons while PIMD lets them move (suspected cause of a 9–25 % tail excess at 1 amu, unproven); 0.1 amu statistics are short (≤ 12 ps per run); pit traps at 4 amu and P=1 are excluded by a fixed `|z| > 1.2 Å` rule. A production-quality set (more beads/seeds, frozen and mobile sets) was requested and was running when this repository was assembled; it is **not** reported here.
- **Performance**: i-PI + socket overhead (≈ 0.21 ms per bead request) dominates; an in-process engine would help but is not done.

## 7. Reproducing the Stage 0–2 tunneling results

(Commands follow the scripts' real argument names. Paths assume the build-root layout; runs take hours.)

```
# 1. Solver: design scan and chosen design, all four masses
python3 inputs/tunneling_ag/schrod3d_barrier.py 0.03 10 solver_chosen.json 3:0.10 --masses 4,1,0.5,0.1

# 2. One PIMD run: 1 amu, P=36, dt 0.125 fs, 20 ps, left start, chosen barrier
python3 inputs/tunneling_ag/run_barrier_lammps.py --out run_m1 --beads 36 --dt-fs 0.125 --steps 160000 \
    --el-mass 1.0 --v0-kcal 3.0 --sigma-z 0.10 --start-z -0.3 --seed 1 --pos-stride 40 --latency 1e-5

# 3. Compare with the solver (spec = dir:P:dt_fs:pos_stride)
python3 inputs/tunneling_ag/analyze_stage2.py --runs run_m1:36:0.125:40 \
    --solver validation/tunneling_ag/data/stage2/conv_3_010_h003_m1.json --mass 1 --skip-ps 1.0
```

Options worth knowing: `--freeze-others` (the solver's model), `--all-rings`, `--no-fix-ag`, `--clients`, `--checkpoint-stride`/`--from-restart`. Stage-2 data and tables are under `validation/tunneling_ag/data/stage2/`. Settings used per mass are in `STAGE2_DESIGN.md` section 2.1 (seeds and run lengths).

## 8. Reports (read these for the evidence)

`validation/phase0_reference/README.md` (exact reference, 2×2 table, graphene findings) · `validation/phase1_ipi_mixed/README.md` (i-PI patch) · `validation/phase3_engine/{README,ADAPTER,TORSION,DISCONTINUITY,SMOOTH_ADOPTION,ROBUSTNESS_2026-11-11}.md` · `validation/tunneling_ag/{A1_A3,A4_A7,PART_B_REVISED,STAGE0_STAGE1,STAGE2_DESIGN}.md` · `docs/PI_eReaxFF_AUDIT_REPORT.md` (findings on the earlier in-house implementation) · `docs/PI_eReaxFF_BUILD_PLAN.md` · `docs/STATUS.md` (dated log).

## 9. License, legacy-derived material, what is not here

- **License: to be decided.** No license is granted by this repository. `pair_style/` carries its own upstream license files (`pair_style/LICENSE`, `LICENSES/`); i-PI patches fall under i-PI's license (GPL-3.0 upstream).
- **Legacy-derived content** (public release undecided): `legacy_derived/` holds small patches against the legacy `ereac` Fortran source (export patch, exact-cell patches, a torsion dump patch). `pair_style/tests/` contains expected values generated with the legacy code. The legacy Fortran sources themselves are **not** included.
- **Kept local, not in this repository**: legacy sources and binaries, graphene/Jamil `reference_data` (`phase0_*`, `phase4_frames`), run directories and trajectories, the in-house engine fork, paper PDFs, LAMMPS/i-PI checkouts and build outputs.
- The Ag/El `ffield` follows the published parameters (SI of the paper below).

## 10. References

- Evangelisti, Fichthorn, van Duin, *J. Chem. Phys.* **153**, 104106 (2020): eReaxFF Ag/El parameters and SI.
- i-PI: https://github.com/i-pi/i-pi · LAMMPS: https://www.lammps.org · Metal-Reax: https://github.com/yeatanir/Metal-Reax
