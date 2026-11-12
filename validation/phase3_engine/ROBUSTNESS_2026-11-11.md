# Robustness items (a)-(d), 2026-11-11

Companion to `DISCONTINUITY.md` and `SMOOTH_ADOPTION.md`. Labels: [measured] run and seen.

## (a) One command for every check: `python3 tools/run_all_checks.py` [measured]
11 checks, PASS/FAIL/SKIP per line, exit status 1 on any FAIL, every LAMMPS check logs `CONTROL_USED`. Last full run: **ALL CHECKS PASSED in 162 s** (quick set; `--full` adds the whole classical-atoms pytest file; `--only name ...`; `--mr-build <dir>` for the eparams check).
| check | what | result |
|---|---|---|
| frames_ag | Ag3/Ag7_1 charged frames, adapter vs exact ereac (legacy) | dE 1.6e-11 / 8.7e-9 kcal/mol, forces 3e-12 rel |
| frames_graphene | 8 periodic C/El frames vs stored ereac | dE 1.3e-5, torsion 6.6e-7 rel |
| fd_forces | adapter force vs central difference of its own energy (legacy and smooth) | 3.5e-5 kcal/mol/Å |
| smooth_vs_ereac | smooth vs ereac with `vlbora 5.2`, `cutof2 1e-9` | dE 5.9e-4 (floor = 1e-5 Å geo quantisation), same as legacy |
| conservation_p1 | Ag7_1 NVE mobile Ag, smooth, dt 0.25 → 0.125 | span 0.0235 → 0.0060 (ratio 3.95) |
| classical_atoms | i-PI classical-atoms quick subset (L-P1, L-P2, restart, exact ring moments, NVE conservation) | 5 passed, 90 s |
| ipi_regression | i-PI's own regression runner, CLASSICAL_ATOMS | 2 passed, 3 s |
| ipi_lammps_short | i-PI + LAMMPS, P=4 NVE 200 steps, legacy, vs stored output | stored (regression baseline only, not an independent reference) |
| schrodinger | 1-D solver re-run vs stored JSON; stored PIMD moments vs solver | diff 0; max 1.2 σ |
| eparams | Metal-Reax unit test + LAMMPS test (c) | pass |
| periodic_ag | periodic Ag7_1 straddling the boundary, 22 Å, vs ereac | dE 8.8e-7 kcal/mol |
Not covered: PIMD-versus-solver on *new* runs (only the stored 1 amu result), the full classical-atoms file (`--full`), graphene long-time conservation.

**The "hung" i-PI regression runner was not a hang.** Cause found [measured]: the reference outputs `ref_simulation.out` of the CLASSICAL_ATOMS tests were never in the branch (`.gitignore` has `*.out`), so the runner died with FileNotFoundError after the run and, when started over the whole suite, appeared to hang. With the two files added (`git add -f`, commit `cb9b4db` on fork branch `regression-runner-diag`, worktree `scratch/ipi_reg_diag`, based on `classical-atoms` 071f6b9) the runner passes both tests in 3 s through pytest (`-k CLASSICAL`) and as a script (`test_run.py -f CLASSICAL_ATOMS`). The earlier "hung for >5 min" was not reproduced; I did not find a hang.

## (b) Torsion fork folded into the main line [measured]
`metal-reax-torsion` contains no engine change: its only commit beyond `ereaxff` (08cdbf3) is a test-harness switch `REAXMETAL_LEGACY_PI` (5 lines in `tests/test_ereaxff_reference.cpp`); the torsion residual is a legacy artefact (TORSION.md). It was cherry-picked (02e20e9) onto a new branch `ereaxff-lammps-merged` (worktree `scratch/mr_merged`, from `ereaxff-lammps` 120b609). Full ctest in that tree: **64/64 pass** (before the cell check was added); 63/64 after (the one remaining test, `lammps_charge_models`, takes 2–3 min on a quiet machine and was stopped under load; it does not use eReaxFF mode). The branch used by the running jobs was not touched.

## (c) Species tables out of code [measured]
New pair_style keyword `eparams <file>` (commit e97c9a8 on `ereaxff-lammps-merged`): lines `SYMBOL Z SHARE_SIGN tested|untested` override the built-in znucel / elecpos tables per element. Accepted in eReaxFF mode: C, El, Av (validated) or any element declared `tested`; `untested` and absent elements are refused with specific messages. A new element therefore needs a parameter file, plus the validation that justifies marking it `tested`, not a code change. `tests/data/legacy_table.eparams` writes the built-in tables as a file; using it reproduces the built-in run exactly (energy and eelect identical to all printed digits). Tests: unit (`test_compat_flags`: parse, override, 6 malformed cases) and LAMMPS (`tests/lammps/run_eparams.py`: Av's Z changed 3 → 3.5 changes eelect −175.44 → −212.91; H refused / refused-as-untested / accepted-as-tested with a finite energy; malformed file rejected). Not done: the Metal GPU backend (unused in eReaxFF mode); the derivative-side sign tables of the legacy `poten.f` have no counterpart in the engine (one sign per element), so an element whose legacy energy and derivative tables disagree (the audit's R10 quirk) is not detected by the file.

## (d) Periodic Ag/El vs ereac, and the cell-size refusal [measured]
Ag7_1 (10 frames) placed straddling the periodic boundary, cubic cells, adapter vs the exact `ereac` with CRYSTX (`tools/discontinuity/periodic_ag_check.py`):
| box edge (Å) | max \|ΔE\| (kcal/mol) | terms / torsion (rel) | force (rel) |
|---|---|---|---|
| 22 | 8.8e-7 | 2.6e-8 / 1.8e-8 | 1.1e-7 |
| 16 | 6.4e-7 | 1.9e-8 / 1.3e-8 | 8.1e-8 |
| 12 | 4.7e-7 | 1.4e-8 / 9.9e-9 | 6.1e-8 |
| 11 | 4.3e-7 | 1.3e-8 / 9.1e-9 | 5.7e-8 |
| **10** | **4.0e-2** | 1.2e-8 / **1.1e-3** | **1.4e-2** |
| 9 | 2.8 | 1.1e-8 / 6.3e-2 | 0.11 |
| 8 | 5.2 | 1.0e-8 / 0.12 | 0.18 |
The adapter agrees with the legacy code down to 11 Å and disagrees (torsion/valence, i.e. the bonded terms that use the legacy single minimum image `dista2`) from 10 Å. **This is 2 × bond_cut, not 2 × the 10 Å non-bonded cutoff**: with the legacy pair list (`nvlx/ix`) the non-bonded terms handle several images correctly, as the validated 12.3 × 12.8 Å graphene cell already showed. A refusal at 2 × 10 Å = 20 Å would also have rejected that validated cell. I therefore implemented the refusal at **edge < 2·bond_cut + 1 Å** (11 Å for the legacy 5.0 Å; 11.4 Å for smooth), eReaxFF mode only, plus refusal of periodic triclinic cells (the legacy reader has no tilt): this deviates from "twice the cutoff" on evidence; change `minlen` in `init_style` if you prefer the literal rule. LAMMPS tests: edge 10 Å refused, 12 Å accepted. The 12.3 Å graphene cell and all other periodic ctest cases still pass.

## Graphene, one longer adapter run (low priority)
60 C/60 El periodic (phase4_frames), legacy thresholds, adapter, 298 K, dt 0.05 fs, 8000 steps (400 fs) [measured]. Total energy: 0–100 fs span 0.029 kcal/mol (conserved); then the system collapses: KE mean 169 → 946 (100–200 fs) → 3500 kcal/mol, total-energy span 137 kcal/mol, and from 110.7 fs on 247 steps change the total energy by more than 0.5 kcal/mol (largest −4.7 at 118.6 fs, KE changing by 1; the terms moving are Coulomb `eco`, bond `eb` and the electron term `eelect`). This reproduces the collapse at 110–120 fs already seen with dt 0.1 (ADAPTER.md §5). It is a loss of energy conservation during an electron/core collapse, not the Ag-type cutoff jumps (those are 0.017 kcal/mol); whether the collapse is physical for this force field (electrons falling onto cores) or an integration failure at 0.5 amu electrons was **not diagnosed**. Nothing more was done on graphene.
