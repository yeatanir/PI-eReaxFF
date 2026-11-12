# Torsion residual Metal-Reax vs legacy ereac (task 0.4) — work in progress

Worktree `forks/metal-reax-torsion`, branch `torsion-baseline` (from `ereaxff` e4b53df). Temporary debug code was saved as patches and removed from the engine.

## Cause found (carbon): the legacy code, not Metal-Reax, is inexact
Per-torsion dumps of both codes on sharing case A (42 torsions) show identical torsion sets and identical bond orders, distances and r4 (1e-14). The only difference is the valence angles `ht1`/`ht2` computed by legacy `calvalres` (poten.f:5258): they differ from the exact angle by up to 1.7e-7 rad, but only for triples containing atom 12, and are exact to 1e-16 for all other triples. Torsion = f(arg), arg built from those angles, so the torsion energy inherits it (1.7e-6 abs).
Why: legacy `itrans=1` translates atoms into the box (case A: atom 12 sits at x = 197.2064 in a 200 A box), and `dista2` (reac.f:5464-5510) takes the minimum image with `dx - anint(dx/aaxh)*aaxh`, `aaxh = sqrt(tm11²+tm21²+tm31²)` (reac.f:5453). `tm11 = axis1*sinbet*sinphi` is built from angles converted with `dgrrdn = 1/(180/pi)` and the truncated `pi = 3.14159265d0` (reac.f:45), so `aaxh` is not exactly 200 (error ~1.4e-7 A). A pair straddling the wrap therefore gets a displacement wrong by ~1.4e-7 A in `calvalres`/torsion; the bond list (Verlet) uses another path and is exact.

Evidence (sharing case A, non-periodic):
| legacy variant | torsion diff vs Metal-Reax (rel) | torsion_conj diff (rel) |
|---|---|---|
| shipped legacy | 2.9e-6 | -3.5e-8 |
| same geometry shifted +6 A in x,y (no atom crosses the wrap) | **2.5e-10** | **7.3e-11** |
| legacy rebuilt with exact pi (reac.f:45) | 2.5e-10 | 7.3e-11 |
Total gradient max diff 4.2e-5 -> 7.5e-8 kcal/mol/A for the shifted case.
So Metal-Reax's torsion formula, cut-offs and literals agree with the legacy ones; the 3e-6 residual is a legacy artefact of wrapped atoms. No Metal-Reax change is required; the torsion tolerance of the reference tests should be tightened using references that avoid the artefact.

## Other findings
- With exact pi the legacy `valence+penalty` moves by 4.6e-7 (case A) away from Metal-Reax, because Metal-Reax's `deg2rad` deliberately uses the truncated constPI (LAMMPS convention, bonded.cpp:53) while exact-pi legacy no longer does. The reference tests must therefore use references built consistently (shipped legacy: valence 7e-9; torsion artefact only for wrapped atoms).
- Reference sets generated: `reference_data/phase0_exactpi/*` (14 cases, pi exact). `phase0_exactcell` (only the cell angles use the exact conversion, pi elsewhere unchanged) was being generated when interrupted.
- Test harness: `REAXMETAL_LEGACY_PI=exact` switches the box construction in `tests/test_ereaxff_reference.cpp`.
- Patches: `reference_data/legacy_patches/{legacy_exactpi,legacy_exactcell,legacy_torsion_calvalres_dump}.patch`, `validation/phase3_engine/torsion_dump_debug.patch`.

## Final result (exactcell references)
`reference_data/phase0_exactcell/*` = legacy built with only the cell-angle conversion exact (patch `legacy_exactcell.patch`, 42 lines in reac.f; everything else, including pi at reac.f:45, unchanged). Harness run with `REAXMETAL_LEGACY_PI=exact`, **no term skipped**, tolerance 1e-8:
- 14/14 cases pass (225 and 60/60 atoms, periodic and not, sharing cases A-J). Torsion and torsion_conj <= 2.5e-10 / 1.2e-10 relative, valence <= 4.4e-10, total gradient max diff 6.6e-12 .. 2.1e-7 kcal/mol/A (largest value ~285).
- Full ctest in the worktree: 36/36 pass (-j1). (The main `ereaxff` tree reported 38; the difference is not investigated, probably the REFDIR-gated tests registered differently in this build.)
Recommendation: validate against exactcell references (or shifted geometries); keep the shipped-legacy references only as regression data with the documented 3e-6 torsion artefact for atoms crossing the cell wrap. No change to Metal-Reax's torsion code is needed.

## Not done
- Ag/El (`Av`) plain-ReaxFF comparison (coordinator request): **no Ag force field exists on this machine** (searched Margaret/*), so nothing was run. Needs the Ag ffield; then run Ag-only periodic/non-periodic with run_ereac_case.py + the exactcell binary (`scratch/phase0/ereac_exactcell_build/ereac`). The angle record with theta0 = 179.99 (sin ~ 1e-3) is the case to watch against Metal-Reax's kMinSine handling.
- The plain C/H/O molecule comparison (no electrons) was not built; the electron-free limit is covered only via pxel2 = 0 cases.
- Exact-pi references are kept in `phase0_exactpi` but are not recommended (valence moves 4.6e-7 away from Metal-Reax).
