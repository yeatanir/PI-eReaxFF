# Phase 3 — eReaxFF in Metal-Reax (CPU-64)

Fork: `../../../forks/metal-reax`, branch `ereaxff` (cloned from yeatanir/Metal-Reax @ 735ad6a). Legacy source is NOT copied into the fork;
formulas are re-implemented with comments citing `Margaret/eReaxFF_src/{reac.f,poten.f}` file:line.

## Step 1 — baseline (done)
`cmake -S . -B build -DREAXMETAL_ENABLE_METAL=OFF -DREAXMETAL_BUILD_LAMMPS_PLUGIN=OFF`, build, `ctest`: **18/18 pass** (3 s).
No LAMMPS was downloaded; tests needing `REAXMETAL_FFIELD_DIR` / the LAMMPS prefix are not registered in this configuration (ffield_tables, neighbor_fixtures, lammps_ghosts, lammps_probe, in-LAMMPS tests).

## Facts established for the port (from reading the code)
- Atom block columns (reac.f:4516-4524, format 1250 = 3x,10f9.4): line2 = alf, vop, valf, valp1, **pxel2**, chi, eta, vnphb; line3 = vnq, vlp1, **vincr**, bo131, bo132, bo133, softcut, **alphag**; line4 = vovun, vval1, **betael**, vval3, vval4, rcore2, ecore2, acore2. The Metal-Reax parser skips line2 col5, line3 col3 and col8, line4 col3.
- Bond block line 2 (reac.f:4668-4670): psp, pdp, ptp, bom, bop1, bop2, ovc, **pxel1** (8th column; the parser reads up to ovc).
- Electron pair term (poten.f:3924-4160): core-El `-c1c*betael(core)*Z(core)/r*erf(sqrt(2*alphag(core))*r)`; El-El `+c1c*betael(i1)*erf(sqrt(alphag(i1))*r)/r` (parameters of atom i1 only); c1c = 332.0638; both times the 7th-order taper built from (vpar(12), vpar(38)); pairs with r > vpar(38) or r < 1e-10 skipped; each pair once, self-images excluded. Ho (hole) branches exist; refused here.
- Z (`znucel`) table is hard-coded in reac.f:4590-4615; sharing sign table in reac.f `elecpos` (≈8815-8851 and the two later copies at ≈8900/8905-8935, whose lists differ slightly: the later copy has no Lv/Av and has `Al`, `Os` — to be audited in task 3.3).
- The shipped reference_data ffield (`np_export/ffield`) has general parameter 38 = 10 (the scratch edit used for the exact reference); `jamil_vpar38_0/ffield` has 0 as in the repo original.

## Plan decisions
- Do not add `Electron` to the 13-term `EnergyTerm` enum yet (invariant: every term owned by exactly one LAMMPS pvector slot, tested). The electron term has its own result struct; the 15th slot belongs with the LAMMPS adapter (task 3.5).
- Electron taper is an explicit option (`FfieldOptions::electron_taper`, NaN = use gp[37]); in eReaxFF mode a resolved taper <= 0 is rejected.

## Status log (2026-11-11)

Fork `forks/metal-reax`, branch `ereaxff`, commits: `5dae7b0` (3.1), `c8cb917` (3.2), `3a86567` (pure atom-energy functions), `5e01724` (engine: sharing, lone pair, over/under),
`bc633e4`, `9c9646c` + one more (exact-reference harness, legacy constants, periodic cell, sharing cases). CTest: **38/38 pass** (CPU limit respected: build -j4, ctest -j2).

| Task | State | Evidence |
|---|---|---|
| Step 1 baseline | done | 18/18 before any change |
| 3.1 parser + data model | **done** | `ereaxff_ffield` (synthetic + real Jamil-target ffields) |
| 3.2 electron pair term | **done** | closed forms 1e-12; FD 1e-7; exact reference to 6.6e-16 relative (energy) and 1e-14 (gradient, stage export) |
| 3.3a electron sharing (e2, A = avalc, B = stlpc, sign table) | **done** | `ereaxff_bonded` FD (below) and the ten sharing cases against the legacy executable |
| 3.3b modified lone pair + C2 correction | **done** | FD; legacy: lone pair <= 1e-10 relative in cases D-J (C lone-pair strength 10), 1.5e-10 on the 225-atom case |
| 3.4 modified over/under-coordination (pxel1/pxel2 factor, all channels) | **done** | FD; legacy: over+under 3e-9 absolute (2.5e-11 to 8e-11 relative) in all ten sharing cases; the pxel2 effect (A minus B = 0.0134 kcal/mol) reproduced to 1e-9 |
| Energy breakdown | partly | lone pair / over / under fill the existing `EnergyTerm` slots. The electron energy is NOT an `EnergyTerm`: adding one breaks the tested 13-term / 14-slot LAMMPS pvector contract (each term owned by exactly one slot); it needs the 15th slot of task 3.5. It is returned in `ElectronResult` and printed by the harness. |

### Finite-difference tests (`tests/test_ereaxff_bonded.cpp`, `test_ereaxff_terms.cpp`, `test_electron.cpp`)
Central differences of the energy against the analytic gradient; synthetic ffield with invented numbers; 3 random configurations each.

| Variant | cluster | periodic cell (ghost gradients folded) |
|---|---|---|
| lone pair only (A,B channels) | 8.7e-9 | 5.7e-7 |
| lone pair only, pxel2 = 0 | 2.9e-8 | 1.1e-8 |
| under only | 1.8e-7 | 3.3e-8 |
| over only (+bond) | 3.5e-7 | 2.6e-6 |
| over only, pxel2 = 0 | 9.8e-7 | 3.1e-7 |
| all terms | 3.3e-7 | 5.4e-7 |
| all terms + C2 correction | 5.4e-7 | 6.9e-7 |
(kcal/mol/A; gradients up to ~800; electrons change the bonded energy by 79.5 kcal/mol in the test system, so the test is not vacuous). Pure-function partials (lone pair, over/under, corrected atom): 3.5e-7 / 4.9e-9.

### Exact-reference results (`tests/test_ereaxff_reference.cpp`; expected terms and gradient read from the rebuilt-legacy export, kcal/mol)
Relative differences per term (largest per case); the cell for periodic cases is built exactly as the legacy code does (see constants below).

| Case | bond | over+under | lone pair | valence | torsion | torsion_conj | vdW | Coulomb | electron | total E | gradient (max abs / max ref) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 225 atoms, non-periodic | 9e-12 | 3e-10 | 1.5e-10 | 3e-11 | 7e-11 | 1e-10 | 1e-16 | 2e-14 | 7e-16 | 4e-8 abs | 2.1e-7 / 285 (electron stage 1.6e-12, nonbonded stage 1.9e-12, bonded stage 2.1e-7) |
| 60C/60El, non-periodic | 7e-12 | 2e-10 | - | 7e-12 | 3e-11 | 1e-10 | 4e-16 | 2e-14 | 1e-15 | 2e-8 abs | 1.2e-8 / 195 |
| 225 atoms, periodic | 2e-11 | 1e-9 | - | 1e-10 | **7e-8** | 1e-10 | 1e-16 | 1e-14 | 2e-15 | 2e-6 abs | 1.1e-5 / 273 |
| 60C/60El, periodic | 2e-11 | 1e-9 | - | 4e-10 | **6e-7** | 2e-9 | 6e-16 | 2e-14 | 8e-16 | 8e-6 abs | 3.5e-6 / 126 |
| sharing A-J (non-periodic) | 1e-11 | 8e-11 | <= 1e-10 | 7e-11 | **3e-6** (skipped) | 3e-8 (skipped) | 3e-15 | 1e-14 | 3e-14 | 2.6e-6 abs | 4e-5 / ~200-410 (same value in all cases that share the carbon geometry: torsion) |

Sharing cases (`project/reference_data/phase0_sharing/A..J`, generator `tools/make_sharing_case.py`): C ring + 2 carbons, 1-5 electrons 0.05-0.7 A from carbons, vincr(C) 6 or 181, pxel2(C) 0 or 0.8, C lone-pair strength 0 or 10, fixed charges. In these cases e2 reaches 0.16-0.9, so they exercise the new terms strongly (the 225-atom and 60C/60El reference geometries have electrons 0.35 A from carbons: e2 ~ 2e-10, essentially a stock-ReaxFF check).

### Legacy compatibility constants found (each flagged in the code with the legacy file:line; policy: check for REAL*4 literals first)
| Constant | Legacy | Effect if ignored | Where |
|---|---|---|---|
| `c1c = 332.0638` is REAL*4 -> 332.0638122558594 | poten.f:3290, 3940 | Coulomb and electron energies off by 3.1e-7 / 3.7e-8 relative (LAMMPS uses 332.06371) | `electron.cpp`, `NonbondedOptions::coulomb_constant` |
| `fothird = 4.0/3.0` evaluated in single precision (1.3333334) in the Coulomb derivative | poten.f:3292, 3347 | nonbonded gradient off by 6e-8 relative (the legacy force is not exactly the gradient of its energy) | `NonbondedOptions::legacy_coulomb_exponent` |
| `0.040` and `0.16` REAL*4 in the C2 correction; gate `abs(vpar(6)) > 0.001` | poten.f:2118, 2135, 2145 | lone pair off by 3e-7 relative | `ereaxff_terms::c2_correction` |
| `pi = 3.14159265d0` (truncated) | reac.f:45 | electron gradient off by 9e-10 relative; periodic energies off by ~1e-10 through the cell construction (a 90 degree cell has ~2e-8 A off-axis components) | `electron.cpp`, harness cell |
| heavy-atom threshold `amas > 31` (LAMMPS: mass > 21) | poten.f:887 | only for elements between 21 and 31 (Mg, Na, Al, Si, P); equal for C, H, O, N, El | parser, eReaxFF mode |

### Open items (not hidden)
1. **Torsion differs from the legacy value** in small clusters (3e-6 relative, 1.7e-6 absolute at 0.57 kcal/mol, with a 4e-5 kcal/mol/A gradient effect) and in periodic cells (7e-8, 6e-7). It is identical in cases with and without electron sharing, and for pxel2 = 0, so it is independent of the eReaxFF terms: a stock Metal-Reax (LAMMPS-flavour) versus legacy torsion difference. It is not caused by `thb_cut` (0.001 is right; other values are far worse) nor by the `thb_cutsq` product threshold. Unresolved; belongs to task 0.4 (plain-ReaxFF baseline). The torsion terms are explicitly skipped in the sharing tests; the 225-atom non-periodic case matches torsion to 7e-11, so it depends on geometry.
2. Bonded energies carry ~1e-11 relative (bond) to 3e-10 (over+under) residuals and the bonded gradient ~1e-9 relative against the legacy executable, also in the stock-only limit (electrons far away). Not analysed further.
3. The legacy electron-sharing distance is a minimum image over all atoms; here the sum runs over far-list pairs (r <= nonb_cut) including ghost images. Identical for vincr = 181 and 6 at the distances tested; differs for tiny vincr.
4. Three copies of `elecpos` exist in the legacy source with slightly different element lists (first copy used here: C,O,H,F,N,Os negative; Li,Lv,P,B,Av,Al positive). Only C and El were tested.
5. Not done: FP32 twin and Metal kernels for any eReaxFF term (the engine throws if built in float), per-atom energy/virial tallies in eReaxFF mode, LAMMPS adapter (task 3.5), electron term as an `EnergyTerm`.
6. Generator caveat found and fixed: the legacy ffield reader is fixed-column (`3x,10f9.4`); a one-character shift clipped `181.0000` to 81 in the legacy executable while the token parser read 181. Any tool that edits ffield lines must write the exact columns.

## Ordered to-do list for tasks 3.3 and 3.4 (what the legacy code does, where it goes in Metal-Reax)

Metal-Reax insertion points (src/cpu/bonded.cpp, 802 lines): bond list and corrected bond orders `:108-165`; atom quantities `:165-220` (`terms::atom_quantities`, include/reaxmetal/terms.hpp:137-164); Bonds `:220`; Atom_Energy lone pair `:256-263`, over/under `:290-320`; Valence_Angles `:322` (uses `aq[j].nlp`, `dDelta_lp` at `:344`); bond-derivative gather `:718`; the Metal/FP32 twins (`bonded_fp32.cpp`, `bonded_device.cpp`, shaders) are NOT to be touched until the CPU-64 version is validated.

1. **eshare pre-pass** (legacy `elecpos`, reac.f:8779-8851, later copies ~8870-8990): for every El atom e and every non-El, non-Ho atom i (ALL pairs, no cutoff except `eshare > 1e-12` in the later copies at reac.f:8872 and ~8905), `eshare = exp(-vincr(type_i) * r^2)`; `eshare2(i) += eshare`; `avalc(i) = aval(type_i) + share_sign_i * eshare` ; `stlpc(i) = stlp(type_i) + eshare` for the elements in the sign table (others unchanged but still accumulate eshare2). `r` from `dista2` = minimum image (V-10), so for periodic cells the sum is NOT over all images. Decide: reproduce minimum image (compat) or ghost-native image sum; they differ only for tiny vincr. Audit the three copies: they differ in the element lists (first copy has Lv, Av; the later copy lacks them and uses `Al`, `Os`), so identify which copy each call path uses (`calval`, `lonpar`, `ovcor`, force routines) before coding. Data: `eshare3(i,e)` and `rr3(i,e)` are stored for the derivative loops.
2. **lone pair** (poten.f:397-619, `lonpar`): uses `voptlp = 0.5*(stlpc - avalc)`, `vund = abo - stlpc`, `vlp = exp(-vpar(16)*(2+vlpex)^2) - int(vund/2)`; `dvlpdstlpc = +vpar(16)*2*(2+vlpex)*expvlp` is the new derivative channel (stock has only d/dBO). The valence-angle lone-pair adjustment keeps the UNCORRECTED `stlp` (`vlpnc`, poten.f:2351-2361, V-18). Map: `terms::atom_quantities` needs `avalc`/`stlpc` instead of `valency`/`valency_e` for Delta, Delta_e, nlp, Delta_lp (not Delta_val / Delta_boc), plus electron-position force terms from `eshare` (d eshare/dx = -2 vincr * eshare * dx).
3. **C2 correction** (JAX port: `reaxff_lone_pair_overcoordination_terms`, last block) is already in Metal-Reax for stock ReaxFF (compat.c2_correction); check that its dependence on `abo`/`avalc` is the same in the legacy source.
4. **over/under-coordination** (poten.f:862-1700, ~800 lines with derivatives): `vov1 = exp(-pxel2(i)*eshare2(i)*sumov3(i)) * (abo(i) - avalc(i)) - dfvl*(voptlp(i) - vlptemp(i))*vho` and the neighbour sum `sumov += (bopi+bopi2) * (exp(-pxel2(j)*eshare2(j)*sumov34(j)) * (abo(j) - avalc(j)) - diffvlp2)` with `sumov3(i) = sum_k pxel1(ibt)*bo/abo(i)` and `sumov34(j) = sum_m pxel1*bo/abo(j)`; heavy-atom switch is `amas > 31.0` in the legacy code vs mass > 21 in Metal-Reax (V-19; equal for C/El). Derivative channels to add: through `abo` (sumov3 and the 1/abo factor), through pxel1-weighted bond orders (sum over neighbours and neighbours of neighbours: ghost shell requirement grows), through `eshare2` (electron positions), through `avalc`. Validate each channel separately with finite differences before combining; compare against the JAX implementation as a third reference.
5. **Totals check**: with 1-4 in place, compare `Ebond ... Eelect` and the full gradient (reference_data ereac_export.out, step 0, non-periodic) in a `reaxmetal_energy_tool`-style harness; target: all 14 terms to ~1e-10 relative and gradient to ~1e-9; then the periodic case against `jamil_vpar38_10` (electron 21829.68).
6. **Adapter (task 3.5)**: keyword `ereaxff yes`, 15th pvector slot, fixed-charge mode (no polarization energy when the legacy `ech` is excluded), refuse `backend metal` and Ho. Fixed charges: legacy `icharg 8` sets charges from MOLCHARGE and the energy total adds no charge energy (`ncha2*ech`, reac.f:5383).
