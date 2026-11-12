# Phase 5 (Ag/El PIMD) - working notes

Engine for these tests: exact rebuilt `ereac` called through `tools/ipi_ereac_driver.py` (one single-point per bead, parallel, clean directory per task). Launcher: `inputs/ag_el/run_ipi_ereac.py`.

Bugs found and fixed on 2026-11-12 (all mine):
1. Cluster generator wrote `END` before `MOLCHARGE` -> legacy reader stopped at END -> all charges 0, Coulomb = 0. Found by the LAMMPS-adapter fork. Frames, MD runs and the first PIMD runs were invalid and were redone. Rule: `END` must be the last line of a geo file.
2. Driver reused one working directory across calls: ereac keeps state between runs and returned the same energy for different geometries (finite-difference check gave exactly 0). Now every task starts in a clean directory.
3. Driver used `idx % workers` as directory: concurrent tasks of one batch could share a directory (P = 8, 2 workers) -> 'bad ereac result'. Now one directory per task.
