# Phase 0 — ground truth for eReaxFF

Source: `/Users/yeatanir/PhD/Literature/eReaxFF/Duin_group/Margaret/eReaxFF_src` (read-only; never copied into `forks/`).
Build area: `../../../scratch/phase0/ereac_build` (throwaway).

## Task 0.1 — rebuild `ereac` (partly done)

Finding **P0-1**: the checked-in `reac.f` / `poten.f` do **not** compile against the checked-in full-size `cbka.blk`.
The `electron` setup in `reac.f` (≈ line 4585, "Populate znucel(nsort)") uses `znucel`, which is declared only in
`cbka_pimd.blk`. With the stock `makefile`, gfortran fails with "Syntax error in IF-clause" because `znucel(i1)=...`
parses as a statement function. So the shipped full-size binaries (`ereac*`, built 2026-04-10) were built from an
earlier version of the source, and the only configuration of the current source that builds as shipped is the
small PIMD variant (`makefile_pimd`: nat = 1024, nbomax = 2000, nvamax = 5000, ntomax = 5000).

Workaround used in scratch only: add `$       znucel(nsort),` after the line `$       vkap(nvatym),` of a copy of `cbka.blk`.
Result: stock `makefile` then builds with 0 errors (gfortran 15.2, `-O3 -march=native`, arm64).
Binary hash (first 16 hex digits of SHA-256): `c9963253bb299d74` (the shipped `ereac` is `32a4bf18b4478de1`, so they differ, as expected for a different source state).

Source hashes (first 16 hex digits): `reac.f 17ab34a5637a8a9f`, `poten.f 1f0f9274c5c1e39b`, `cbka.blk 533eeff2d9760d6a`.

## Task 0.3 — the 2×2 table (done, 2026-11-11)

Root cause of the earlier `Energy not a number`: the standalone code reads the cell line with a fixed format
`(8x,6f11.5)` (`reac.f:11481`), but the Jamil `geo` writes `CRYSTX` with 13-character fields, so every field after the
first is misread (angles come out as 90, 0, 90 and the distance routine divides 0/0). The V2 oracle had its own geometry
reader, which is why it worked. Fix (scratch inputs only): write the line as `CRYSTX` + 2 spaces + six `f11.5`.
The generator `projects/jamil_replication/generate_jamil_graphene.py` of the old tree probably needs the same fix (not checked).
Credit: found by an independent review of the progress report; I verified the format against the source and the result below.

Inputs: V2 `control`/`ffield`/`geo` of the 225-atom Jamil geometry (`icharg 8`), force field supplied as `fort.4`, geometry as `fort.3`.
"Non-periodic" = `CRYSTX` line removed and `axis1..3` set to 80 Å. `vpar38` = general parameter 38 (line 40 of the ffield),
which the base code uses as the electron taper radius.

Reproduce: `python3 tools/compare_phase0.py` (data in `reference_data/phase0_jamil_225atom/`, kcal/mol, fort.73 prints 2 decimals):

| term | V2 saved | periodic, vpar38=0 | periodic, vpar38=10 | non-periodic, vpar38=0 | non-periodic, vpar38=10 |
|---|---:|---:|---:|---:|---:|
| bond | -20248.81 | -22065.17 | -22065.17 | -20248.81 | -20248.81 |
| overcoordination | -693.20 | -341.27 | -341.27 | -693.20 | -693.20 |
| lone_pair | 0.73 | 0.00 | 0.00 | 0.73 | 0.73 |
| valence_penalty | 646.02 | 172.46 | 172.46 | 646.02 | 646.02 |
| torsion | 20.79 | 27.20 | 27.20 | 20.79 | 20.79 |
| torsion_conjugation | -1105.34 | -1419.27 | -1419.27 | -1105.34 | -1105.34 |
| vdw | 3150.05 | 3501.92 | 3501.92 | 3150.05 | 3150.05 |
| coulomb | -10823.98 | -11926.59 | -11926.59 | -10823.98 | -10823.98 |
| electron | 16563.55 | 0.00 | 21829.68 | 0.00 | **16563.55** |

(molecule, angle_conjugation, hydrogen_bond, charge, external_field are 0.00 in all columns.)

**Result.** The rebuilt `ereac`, run **non-periodic with an electron taper radius of 10 Å**, reproduces all 14 saved V2 energy terms
to the printed precision (largest difference 0.0035 kcal/mol, i.e. rounding of a 2-decimal print). This settles two plan findings by reproduction:
- **V-12 confirmed:** the saved V2 data are non-periodic (the `CRYSTX` cell is ignored). Physical periodic energies differ by about 9 % in bond energy and by a factor 3.7 in the valence term.
- **V-13 resolved:** V2 effectively used a 10 Å electron taper; the Jamil-target ffield's `0.0` in parameter 38 gives zero electron energy in the base code. The ffield that V2 actually used had 10 Å in effect (via its wrapper's `tap7th` call). Which ffield is authoritative is still an owner question.

Caveats: single geometry; energies only (forces not compared; `fort.73` prints 2 decimals). Task 0.2 (full-precision export of energies and forces) is the next step to make this a strict test.
Run-to-run note: the rebuilt binary sometimes fails to start with a `dyld ... syscall to map cache into shared region failed` error in this environment; rerunning succeeds. The cause is not understood.

## Task 0.2 — full-precision export (done, 2026-11-11)

Observation-only patch `reference_data/ereac_export.patch` (13 lines, applied to a copy of `reac.f` in `scratch/phase0/ereac_export_build`):
after `estrc` is computed it appends, per MD step, the 14 terms, `estrc` and the gradient `d(3,na)` to `ereac_export.out` in `es24.15`
(Fortran unit 197). It changes no computation. Patched `reac.f` hash `2e3c647c3f46727f`, rebuilt binary `d9b08ccb0d032a60`.
Runner with exit-status checks, freshness checks and retry of start-up failures: `tools/run_ereac.py`.
Strict comparison: `python3 tools/compare_phase0_full.py reference_data/phase0_jamil_225atom/np_export/ereac_export.out`.

Result (225 atoms, non-periodic, electron taper 10 Å, step 0):
- All 14 energy terms and the total equal the saved V2 values to machine precision: **largest difference 2.3e-13 kcal/mol** (total: identical to all printed digits, -12490.195321194020).
- Forces: saved V2 `force_work` = -c × gradient with a single constant c = 4.183975540e-4, **maximum residual 9.4e-16 (7.9e-15 relative)**. Largest gradient 284.6 kcal/mol/Å.
- So the rebuilt `ereac` is a **bit-level reference for V2 on this geometry** (fixed charges, non-periodic).

**Correction to the audit report (F8):** c = 0.043364 × 9.6485e-3 = 4.18397e-4, so the V2 energy and force conversions are consistent with each other
(both use the truncated kcal→eV factor 0.043364). The audit's claim of a ~5e-6 mismatch between energy and force conversion was wrong.
What remains true: 0.043364 differs from the exact 0.0433641153 by 2.7e-6 relative (absolute energy scale only).

Limits: one geometry, non-periodic, fixed charges. Still wanted (per review): the periodic case and a second geometry (60 C / 60 El).

## Still to do

- 0.2b export for the periodic case and the 60 C / 60 El geometry; then force comparison there (no V2 reference exists for periodic).
- 0.1 sanity check on `Margaret/eReaxFF_Jan3_2024/` (a complete input set with outputs, per the review; not tried).
- 0.4 plain-ReaxFF comparison of Metal-Reax CPU-64 with `ereac`.
- 0.5 owner decisions (§8 of the plan): authoritative Jamil ffield and its electron taper; periodic vs non-periodic target.

## Task 0.2b — periodic Jamil case and the 60 C / 60 El geometry (2026-11-11, last step of the session)

Data in `reference_data/phase0_jamil_225atom/periodic_export/` and `reference_data/phase0_p60c60el/{np,periodic}_export/` (exact `ereac` export, electron taper 10 Å, CRYSTX line in the 11-character format).
- Jamil 225 atoms, periodic: total -10221.046 kcal/mol (non-periodic: -12490.195).
- 60 C / 60 El (Track B2 paired graphene, 12.306 x 12.789 x 20 Å, first call of the old P=16 benchmark): periodic -5501.178, non-periodic -6924.582 kcal/mol.
- Saved V2 reference for this geometry (energy -5500.965262, forces) is **close to the periodic** `ereac` result but not identical: energy difference -0.213 kcal/mol; forces fit gradient x c with c = 4.18409e-4 and max residual 4.3e-5 (largest force 0.053 in V2 units). The non-periodic run is far off (-1424 kcal/mol).
- So for this geometry V2 was periodic, whereas for the 112 C Jamil geometry the saved V2 data were non-periodic. The 0.2 kcal/mol residual and the different periodicity between the two V2 runs are **unexplained** (candidates: V2's own geometry/cell reader, image handling of the electron terms, differences between V2 and the base source). To be investigated before the periodic Metal-Reax port is judged against V2.

## Per-stage gradient export (added after the independent review)

`reference_data/ereac_export.patch` now also writes, at step 0, the **cumulative** gradient after each routine of `encalc`
(`STAGE <name>` then `SGRAD` lines): electron, lonpar, covbon, ovcor, valang, torang, hbond, nonbon. Differences between consecutive stages are per-routine gradients; the last stage equals the total step-0 gradient exactly (max diff 0.0).
Data: `reference_data/phase0_jamil_225atom/np_stage/ereac_export.out` (non-periodic, taper 10 A). Binary hash `068efe660b8a19a4`, patched `reac.f` `227981a82c144595`.
Caveat: bond-order derivatives are applied inside `covbon` and `ovcor`, so the `lonpar` stage shows zero gradient and its contribution lands in the `covbon` stage; bond, lone-pair and over/under channels therefore share the covbon+ovcor stages and need parameter knockouts (e.g. pxel2 = 0, vover/vovun = 0) to be separated. Electron, valence-angle, torsion and nonbonded gradients are cleanly isolated.
Review hypothesis to check when the V2 periodicity mismatch is revisited: the 60 C / 60 El adapter control seen in old test directories had `icharg 7` (ACKS2), so the 0.213 kcal/mol difference may come from the charge solve or the 2026 compact-ACKS2 edit rather than periodicity. The control file of the 60 C / 60 El saved run used for these cases has `icharg 8` (fixed charges), so for these inputs the ACKS2 hypothesis is not supported; the V2 run itself may have used another control, which is not recoverable.

## Classical NVE check and realistic-regime frames (done, 2026-11-11)

Method: rebuilt ereac (export build with exact positions, hash `987afd2e2ab97bc1`), 60 C / 60 El periodic, fixed charges. Relaxation by four rounds of ereac minimization (`imetho 1`, 2500 iterations each; energy -5501 -> -6803 kcal/mol, still creeping down slowly, i.e. not a true minimum). Then NVE at dt 0.1 fs (10 000 steps) and 0.25 fs (4000 steps). **Correction (found later the same day):** the control keys `nmethod`/`tset` inherited from the Jamil control file are NOT read by ereac (it reads `imdmet`, `mdtemp`, `iravel`); the runs therefore used the defaults, i.e. NVE (nmethod default 3) **started from rest** (nrand = 0, no random velocities), not at 298 K. The 'T = 298' printed at step 0 is only the label `tsetor`. So the graphene sheet collapsed spontaneously from zero velocity, at the same time scale (~200 fs), which makes the result stronger, not weaker: the flat sheet is not a stable structure in this force field/electron setup. Tools: `tools/ereac_retry.py` (retries only on the dyld start-up error; minimization runs end with exit 2 / `Read past ENDFILE`, which is normal), `tools/geo_from_export.py`, `tools/make_frames.py`.

Findings:
1. **Relaxed electrons sit 0.11-0.12 A from the carbons** (start: 0.35 A) and the sharing sum exp(-181 r^2) per electron is 0.09-0.11 (start: 2e-10), so production-regime geometries exercise the sharing/over-under terms strongly; only the purpose-made cases A..J covered this regime before. Electron relaxation releases about 0.9 eV per electron (-1302 kcal/mol for 60 electrons).
2. **Classical energy conservation of the force field is excellent in the normal regime:** over the first 200 fs (T < 800 K, Epot changing by ~280 kcal/mol) the total energy stays within 0.01 kcal/mol at dt 0.1 fs and 0.06 kcal/mol at dt 0.25 fs (fort.71 prints 2 decimals).
3. **Collapse at ~200-300 fs:** the flat sheet is unstable (carbon z-range 0 -> 2.5 Å at 200 fs -> 11 Å at 300 fs -> 26 Å), Epot falls by a further ~3300 kcal/mol (143 eV), T reaches ~10 000 K, electrons detach (mean El-C distance 0.12 -> 2.5 Å). Energy conservation then degrades (offsets up to -115 kcal/mol at dt 0.1, -141 at dt 0.25; roughly first-order in dt), but this is during an explosive event with unphysical temperatures.
4. **With the ffield exactly as written (parameter 38 = 0, electron term off) the system collapses within 100 fs** (T 298 -> 8600 K, Epot -20023 start) from the original start geometry.
5. Conclusion: neither variant of the ffield in the repo gives a stable graphene sheet at 300 K in the legacy code itself (taper 10: delayed buckling; taper 0: immediate). The earlier heating seen in the in-house engine and in the JAX probe (carbon T 1156-2139 K, 1500 K) is therefore **not an integrator defect**; it is this force-field/electron model (together with an incompletely relaxed start). **The ffield / electron-taper choice must be settled before any stability or physics study.** The Jamil paper ran stable MD, so the ffield in the repo is probably not the one used there.
6. Frames for engine validation: `reference_data/phase4_frames/` (51 exact frames from the taper-10 NVE run, 0-300 fs; energies and total gradient per frame, README inside). Not covered: two electrons near one carbon, an electron midway between carbons.

## Graphene NVE rerun with the correct control keys (2026-11-12)

Same relaxed 60 C / 60 El periodic start, now with `imdmet 3`, `mdtemp 300`, `iravel 1` (random 300 K velocities), dt 0.1 fs, 1 ps, rebuilt ereac (export build). T at step 0 is 300 K as intended.
- Taper 10 A: stable only for ~100 fs; T 300 -> 1270 K in 100 fs, then the sheet collapses (T 7000-11000 K, Epot -7152 -> -10700, ~140 eV released). Total-energy offsets reach -6.6 kcal/mol in the first 100 fs and -153 kcal/mol during the collapse.
- Taper 0 (ffield as written, electron term off): collapses within 100 fs (T 300 -> 8400 K).
- Conclusion unchanged and stronger: with a real 300 K start the flat graphene/El system collapses faster (~100 fs vs ~200 fs from rest). It is not a usable stability test case with either ffield variant; the Ag/El system (published ffield) is.

## Other graphene/El force fields on disk (searched 2026-11-12)

Found with `find` over `eReaxFF PIMD` and `Literature/eReaxFF` (all `ffield*`/`fort.4` files with an eReaxFF header, de-duplicated by content). Families with C + El:
| Family (header) | electron taper (par 38) | C: vincr / pxel2 / alphag / betael | Where |
|---|---|---|---|
| `May6_2016_training` (JPCC 2016 training) | **0.0** | 181 / 0.8 / 22.0 / 0.0027 | `Jamil/eReaxFF/eReaxFF/trainset_JPCC_2016/fort.4`; this is the one our Jamil-geometry work used |
| `Feb 20 2019 with Li/C` (Jamil's Li/B/graphene fits) | **5.7** | 181 / 0.8 / **12.76** / **0.0002** | `Jamil/eReaxFF/B_graphene/B_graphene/Graphene_trainset/**/fort.4` (many copies), `Margaret/eReaxFF_src/test_co2/fort.4` |
| `eReaxFF HAN July 14 2022` (C/H/N/O-type, Yun/Matin) | see files | no C | `Yun/eReaxFF_CO2/...`, `Matin/...` |
| `eReaxFF C/H/Cu Jan 3 2024` | see file | (Al = Cu, B = C, Li = H) | `Margaret/eReaxFF_Jan3_2024/` |
| Ag/El (Evangelisti 2020, from the SI) | 5.7 | Av: see `inputs/ag_el/ffield_Ag_El` | built by us |

Graphene test with the **Feb 2019** force field (60 C / 60 El periodic, rebuilt ereac, dt 0.1 fs, relaxed start (five minimisation rounds, -14537 -> -14861 kcal/mol; electrons 0.155-0.341 A from carbons), random 300 K velocities, NVE, 1 ps):
- No collapse (the May-2016 variant collapsed within 100-200 fs). Temperature stays bounded, 640-1300 K (the relaxed start still releases energy), Epot stays within -15176..-15019.
- **Energy is not conserved**: Etot rises by +106 kcal/mol over 1 ps at dt 0.1 (+66 at dt 0.25), starting after ~100 fs. Not yet diagnosed (candidates: threshold discontinuities in the force field, the truncated-pi cell effect for atoms crossing the box edge, bond-order cutoff events). A dt scan and a no-wrap test are the next steps.
- So the Feb-2019 force field is a better graphene candidate than the May-2016 one, but conservation must be understood before using it for dynamics. Inputs: `scratch/phase0/p60_f2019_*` (not committed).

Graphene (Feb 2019 ffield) energy drift vs timestep (same relaxed start, same random seed; Etot change at t = 400 fs): dt 0.25 fs +17.9, dt 0.1 +28.0, dt 0.05 +25.0, dt 0.025 +12.5 kcal/mol (at 200 fs: +3.5, +2.4, +1.2, +1.4). The drift does **not** vanish as dt -> 0, so it is not a discretisation error of velocity Verlet: it is a force-field discontinuity (energy jumps at bond-order/valence thresholds, e.g. int(vund/2) in the lone pair, bond-order cut-offs) acting while the unrelaxed sheet equilibrates. Locating the discontinuity is still to do (candidate method: log the per-term energy jumps between consecutive steps).

### Graphene energy drift: cause (2026-11-12, Feb 2019 ffield, 60 C / 60 El periodic, dt 0.1 fs, 400 fs, per-step log)
Two separate mechanisms; the rest of the trajectory conserves energy to 0.06 kcal/mol per 400 fs.
1. **Stale pair list.** ereac rebuilds its Verlet pair list every `nvlist` steps (hard-coded default 25, `reac.f:8412`; pairs flagged for bond orders within `vlbora` = 5 A; list radius swb + `vrange` = 2.5 A skin). With the default, all 21 energy jumps fall on multiples of 25 and carry 27.9 of the 27.95 kcal/mol drift. `nvlist 1` (control key `nvlist`) cuts the drift to +9.75; `nvlist 5` puts every jump on multiples of 5. (`irecon` is a different thing - control-file re-read - and has no effect.) **Use `nvlist 1` for ereac reference MD.**
2. **Real discontinuities.** With `nvlist 1`, 15 jumps of 0.5-2.9 kcal/mol remain at random steps (9.6 kcal/mol in total, 0.14 of the rest is smooth). At the largest ones the potential changes are in **Ecoul and Eelect with opposite sign** (e.g. step 1103: Ecoul -4.12, Eelect +3.99; step 3918: Ecoul +3.21): an El-core pair is crossing a cutoff/image boundary where the electron channel (taper 5.7 A, `vpar(38)`) and the Coulomb channel (taper 10 A) handle the pair differently. Pinning the exact pair needs a fine path scan with per-term energies (Metal-Reax harness can do it); not done yet.
Implication for PIMD: ring-polymer beads spread electrons over >0.1 A, so such cutoff discontinuities will matter there; the Ag/El clusters (non-periodic, all pairs within the list) are not affected (energy spans 0.02-0.19 kcal/mol).
