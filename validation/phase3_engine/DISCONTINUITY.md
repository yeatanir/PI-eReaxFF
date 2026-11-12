# Energy discontinuities in eReaxFF dynamics (item 1 of the robustness work) — 2026-11-11

Question: are the dt-independent conserved-quantity wander of Ag7_1 (~0.10 kcal/mol) and the graphene NVE drift one bug, an energy jump at a cutoff?
Labels: **[measured]** run and seen; **[derived]** computed by hand from measured inputs; **[not done]**.
Nothing was changed in any engine, branch, input or running job. All runs used the shared `lmp` and plugin read-only (`build-ns`), 1 core, in `scratch/disc_hunt/`.
Scripts: `project/tools/discontinuity/` (`run_nve.py` per-step NVE + 15 pair terms + x,F dump; `analyze.py` work-residual; `pairs2.py` crossing finder; `scan_pair.py`, `scan_path.py` static scans).

## Answer
**The Ag7_1 wander is a sum of energy jumps, not integration error, and there are two sources, both present in the legacy `ereac` itself (the adapter reproduces them, it does not create them).**
**The graphene drift is NOT explained by them** (see the graphene section).

| # | Source | Term(s) | Jump | Location |
|---|---|---|---|---|
| 1 | Ag–Ag pair leaves the bond list at `bond_cut` = 5.0 Å with BO′ still above `bo_cut` | bond +0.042, over/under −0.032, torsion +0.008 (net +0.0171) | **+0.0173 kcal/mol** per outward crossing (−0.0173 inward); force on the moving atom jumps ~0.3 kcal/mol/Å | d(Ag–Ag) = 5.000 Å |
| 2 | 3-/4-body threshold `thb_cut` (BO > 0.001) | torsion (+0.0145 in the located case); sizes 0.002–0.0145 in the run | up to **0.0145 kcal/mol** | not located to a specific bond (see below) |

## Evidence
**a. Per-step NVE** [measured]. Ag7_1_md frame 0 (charges ±1, taper 5.7), Ag mobile, 300 K, seed 4928459, dt 0.25 fs, 4000 steps (1 ps). Residual r_t = ΔPE + ½(F_t+F_{t+1})·Δx; median |r| 9e-5, but 29 steps with |r| ≥ 0.010 (17 positive, 12 negative).
- 24 of them have |r| = 0.0170–0.0181 and each coincides with one Ag–Ag pair crossing exactly 5.000 Å (positive when the distance increases, negative when it decreases). Same event magnitude, quantised.
- 5 (steps 3355–3413) have |r| = 0.010–0.0147 with no pair at 5.000 Å.
- Sum of r over the events +0.064; the net drift of the total energy over the run is +0.038 (span 0.082).

**b. dt independence** [measured]. Total-energy span for dt = 0.25 / 0.125 / 0.0625 fs: 0.082 / 0.072 / 0.070 kcal/mol (bond_cut 5.0); the cumulative residual at the end is +0.029 at all three. The wander does not shrink with dt, as observed in the P=32 runs.

**c. Static scan, source 1** [measured]. Frame at step 98, Ag0 moved along Ag0–Ag1 through 5.000 Å in 5e-5 Å steps (41 single points, LAMMPS): a single jump at d = 4.99995 Å, PE +0.0171; eb +0.0420, ea −0.0323, et +0.0083, ew −0.0004, eelect +0.0004; force on the moving atom jumps fx −0.09, fy +0.29 kcal/mol/Å.
**Exact reference** [measured]: the rebuilt `ereac` (export build) at d = 4.9998 and 5.0002 gives ΔE = +0.0145 against ≈ −0.003 expected for a smooth curve (slopes −8.0 and −7.1 kcal/mol/Å either side): jump ≈ +0.0174, terms bond +0.0410, over/under −0.0318, torsion +0.0082. Same values as LAMMPS. Also all 61 frames of a 4.0–5.2 Å ladder (step 0.02) and 22 frames of 4.980–5.001 Å (step 0.001) agree between `ereac` and the adapter to 1e-8 kcal/mol, so the function is the same on both sides of the jump.
**Cause** [derived from the ffield]: Av–Av has pbo1 = −0.0139, pbo2 = 8.586, r_σ = 2.4574 Å, so BO′(5.0 Å) = exp(−0.0139·(5.0/2.4574)^8.586) = 0.0021, above `bo_cut` = 0.001 (general parameter 30 = 0.1 /100). BO′ falls to 0.001 only at ≈ 5.06 Å. A bond with BO′ ≈ 0.002 is therefore dropped when the pair passes 5.0 Å. Metal-Reax defaults `bond_cut` to 5.0 (`nbrhood_cutoff`); where the legacy `ereac` sets its 5.0 was not checked, but its energies match the adapter on both sides.

**d. Static scan, source 2** [measured]. dt = 0.0625 run, step 13534 → 13535, 400 single points along the straight line: one jump of +0.01444 at s ∈ [0.3175, 0.3200], entirely in torsion (et +0.01447). In the exact `ereac`, eight points along the same line give torsion 39.2833, 39.2979, 39.2834, 39.2980, 39.2980, 39.2836, 39.2981, 39.2985 (a flip of 0.0145 that toggles with the 1e-5 Å quantisation of the geo file near the threshold); the adapter agrees with `ereac` on all eight to 1e-10. So source 2 is also in the legacy code.
Which bond crosses the threshold was **not located**: the pair-distance list in the jump bin is not discriminating (many pairs at 2.7–5.0 Å).

**e. Removing the thresholds** [measured; diagnostic only, not a fix]. Same trajectory start, via a LAMMPS control file:

| `nbrhood_cutoff` | `thb_cutoff` | dt (fs) | total-energy span (kcal/mol) | max \|r\| | end drift |
|---|---|---|---|---|---|
| 5.0 (default) | 0.001 | 0.25 | 0.082 | 0.018 | +0.038 |
| 5.2 (5.5 identical) | 0.001 | 0.25 | 0.092 | 0.010 | −0.029 |
| 5.2 | 1e-9 | 0.25 | **0.0235** | 8e-4 | +0.010 |
| 5.2 | 1e-9 | 0.125 | **0.0060** | 2e-4 | +0.0006 |
| 5.2 | 0.001, `thb_cutoff_sq` 1e-18 only | 0.25 | 0.092 | 0.010 | −0.029 |

Raising `nbrhood_cutoff` alone removes the 0.0173 events but leaves the wander (its end value changes sign). Lowering `thb_cutoff` as well removes every large event; the remaining span falls by 3.9× when dt is halved, which is the expected dt² integration error. (`thb_cutoff_sq` has no effect here.) The `thb_cutoff` change alters the energy surface (E(0) shifts by +0.28 kcal/mol), so this is a diagnostic, not a recommended setting.

**Scale check** [derived]: ~29 jumps of 0.017 per ps give a random-walk scale √29 × 0.017 ≈ 0.09 kcal/mol, the size of the observed wander (0.10 kcal/mol per bead).

## Graphene [measured, short]
60 C/60 El periodic (phase4_frames, taper 10, El 0.5 amu), 298 K, dt 0.05 fs, 2000 steps (100 fs), adapter: total-energy span 0.029 kcal/mol; no residual above 7e-4 (excluding the box-wrap steps, which my residual does not handle). **No 0.017-sized events.** This is expected: C–C BO′ at 5.0 Å is ≈ exp(−200) [derived], so source 1 cannot occur. Source 2 may exist for carbon (second-neighbour BO near 0.001) but is not visible in 100 fs. Setting `thb_cutoff` → 1e-9 made the graphene run unstable (span 38 kcal/mol) so it is not usable there. The earlier graphene drift (ereac, dt-independent) is therefore **not** explained by these two mechanisms; the previous attribution (stale `nvlist=25` list in `ereac`, which the adapter does not have, plus Ecoul/Eelect cutoff discontinuities) stands, unconfirmed. A 400 fs+ graphene run with the adapter has not been done.

## Impact
- Production tunneling setup (Ag fixed, 4 classical El, one ring): Ag–Ag distances are constant, so source 1 cannot occur; consistent with the absent wander (3e-4 kcal/mol per bead). Source 2 involves El positions and could appear; no event was seen in those runs' conserved quantity.
- Any run with mobile Ag (Ag7_1, silicon, graphene) will have jumps of ~0.01–0.02 kcal/mol per threshold crossing and force jumps ~0.3 kcal/mol/Å; the rate scales with temperature.
- Caveat: one early 4-frame `ereac` set had two stale reference frames (an artefact of my first call into the work directory); it was repeated in fresh directories and only the repeated results are quoted.

## Not done / how I would approach items 2–5
- Source 2: which bond and quadruple (needs BO per bond from the adapter; ask for a debug pvector or run the BO routine in Python from the ffield).
- Fix options for the jumps (decision needed): a) raise `bond_cut` to ≥ 5.1 Å for Ag (removes source 1 and breaks no exactness against `ereac` only if `ereac`'s cutoff is changed too — keep a legacy flag); b) taper/shift BO near `bo_cut`/`thb_cut` (changes the force field; the published Ag/El field was parametrised with the jumps).
- Item 2 (check runner): wrap the existing scripts (frames, finite differences, classical-atoms tests, one short i-PI+LAMMPS run against stored output, Schrödinger comparison) and add these per-step residual checks (events > 0.005 kcal/mol fail).
- Item 3 (torsion fork merge): independent of this; a compat flag.
- Item 4 (element tables): independent; note `kEreaxffElements` allow-list and `znucel`/sharing signs.
- Item 5 (periodic Ag/El): needs the `ereac` periodic path; ghost shell ≥ 2·bond_cut + sharing range (bond_cut 5.2 needs 10.4 + sharing range < 12 comm cutoff).
