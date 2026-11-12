# Prompt — proof-of-concept: barrier penetration of path-integral electrons (Ag / eReaxFF)

## Goal

Demonstrate, as fast as is honestly possible (target: 1–2 days), that the eReaxFF path-integral machinery in `ereaxff_build/` captures **tunneling of the pseudo-electron through a barrier, and its dependence on electron mass**.

This is a stepping stone. The eventual target is electron leakage through a barrier in a small silicon system; there is no silicon eReaxFF yet, so silver stands in. The geometry must therefore be **donor – barrier – acceptor**, not a closed well around one cluster (with nothing outside the wall the electron either stays bound or pays its binding energy to leave, and neither is tunneling).

Everything not needed for this demonstration is deferred: graphene and its energy discontinuity, the ffield/taper choice for the Jamil benchmark, Berendsen, periodic references, the torsion baseline, the i-PI upstream pull request, ACKS2, GPU/batching work.

## What you may assume (already done in this project; re-check only if something disagrees)

- i-PI fork, branch `classical-atoms`: `<classical_atoms>` works and is tested.
- Metal-Reax fork, branch `ereaxff`: all electron terms implemented on CPU-64; LAMMPS adapter passes the Ag/El reference frames at 1e-10 relative; `fix ipi` at P = 1 matches `fix nve`.
- Rebuilt `ereac` with export patch is the exact reference; corrected Ag/El frames are in `project/reference_data/phase4_frames_ag/`.
- Fixed charges are used (no charge solve).

## Part A — prerequisites, in priority order. Do these first; do not start Part B until A1–A5 pass.

**A1. Several beads through LAMMPS + Metal-Reax.** This is the one untested link and everything depends on it.
- Run Ag/El at P = 8 with one LAMMPS client per bead (unix sockets).
- Pass: for 10 bead geometries taken from that run, forces equal the `ereac` reference to the tolerance already used for the Ag frames. Each client uses its own working directory.

**A2. Mass handling.** The demonstration scans the electron mass (1, 0.5, 0.1 amu).
- Set the El mass in i-PI. Confirm in source or by test that the force and energy from the engine do not change when only the dynamical mass changes (the ffield mass column must stay untouched, because the engine derives flags from it).
- Pass: identical forces for the same geometry at the three masses.

**A3. Conserved quantity from an equilibrated start.**
- First establish whether i-PI's printed conserved quantity is divided by the number of beads, and report spans in the same units for every P.
- Protocol: PILE-L equilibration until the bead spread is stationary, then NVE from that state.
- Pass: at P = 8 and P = 32, the span per bead over 2 ps is comparable to the P = 1 span (state the numbers; a factor above about 3 needs an explanation before continuing). Repeat once at half the timestep; the span must fall.
- If any run still uses the `ereac` driver, positions must go in at full precision, not through the `f10.5` file.

**A4. External barrier as a second force component.**
- A small i-PI Python potential (same pattern as `ipi/pes/coupled_harmonic.py`) acting only on El atoms, per bead, unscaled:
  `V(r) = V0 · exp(−s² / 2σ²)`, `s = (r − r_mid) · n̂`
  with `r_mid` the midpoint between donor and acceptor and `n̂` the unit vector along the donor–acceptor axis.
- Combine it with the eReaxFF force component at equal bead count.
- Pass: forces equal finite differences to 1e-8 relative; with `V0 = 0` the trajectory is identical to a run without the component; the NVE check of A3 still passes with the barrier on.

**A5. Map the landscape without the barrier.**
- Static scan with the engine: move one electron along the donor–acceptor axis, all else fixed, and record the energy.
- Purpose: know the intrinsic wells and any intrinsic barrier, so the artificial barrier is sized against real numbers and the two sites are close in energy.
- Deliverable: the profile, with well positions, depths and curvature at each minimum. From the curvature, compute the electron frequency ω at each mass.

**A6. Beads and timestep per mass.**
- From ω in A5: choose P so that (βħω/P)²/24 is below about 1 %, and choose the timestep so that ω·dt is below about 0.2. Expect roughly 3× more beads and 0.3× the timestep at 0.1 amu compared with 1 amu.
- Pass: at each mass, the bead-position histogram along the axis is unchanged within error between P and 2P.

**A7. Analysis script.** Histogram of bead positions along the axis; histogram of centroid positions; fraction of bead density beyond the barrier midpoint; count of centroid crossings; radius of gyration per electron. One figure per mass, classical and path-integral overlaid.

## Part B — the demonstration

**System.** Donor and acceptor made from silver, 5–6 Ag atoms in total. Recommended: two identical small clusters (for example two Ag₃ units from the validated set) facing each other, which gives a symmetric double well for the transferring electron. A cluster plus a single Ag⁺ is acceptable if A5 shows the two sites are within a few kT of each other. Hold the Ag atoms fixed (`fixatoms`) for this first demonstration so the two units cannot merge; electrons are ring polymers. Non-periodic.

**Barrier.** Start from height 3 kcal/mol (about 5 kT at 300 K) and full width at half maximum 0.3–0.5 Å (σ = FWHM/2.355), then adjust using A5. Reason for the width: the free ring radius at 300 K is about 0.20 Å at 1 amu, 0.28 Å at 0.5 amu and 0.64 Å at 0.1 amu, so this range goes from mostly blocking to easily penetrated.

**Runs.** 300 K, PILE-L equilibration then PILE-G production.

| Run | Electron mass | Beads |
|---|---|---|
| classical | 1 amu | 1 |
| PIMD | 1 amu | converged P from A6, and 2P |
| PIMD | 0.5 amu | converged P, and 2P |
| PIMD | 0.1 amu | converged P, and 2P |
| control | each mass | converged P, barrier off |

At least three seeds per row; lengths chosen so that the classical run shows either a handful of crossings or a statistically clear zero.

**Evidence to produce.**
1. Bead density inside and beyond the barrier, classical against path-integral, at each mass.
2. Snapshots of a ring polymer with beads on both sides of the barrier.
3. Crossing counts and the fraction of density on the acceptor side, as a function of mass.
4. The barrier-off controls, to separate the effect of the artificial barrier from the intrinsic landscape.

**Success criterion.** Path-integral density in the classically forbidden region exceeds the classical value by more than the statistical error, and increases monotonically as the mass decreases, with bead count converged at each mass.

## Rules

- Report honestly: mark each result as run, read only, or not done. If a prerequisite fails, stop and report rather than working around it.
- Do not modify `ereaxff_pimd_modular` or the legacy source. New work goes in `ereaxff_build/` (inputs under `project/inputs/tunneling_ag/`, results under `project/validation/tunneling_ag/`).
- Keep machine load within the cap already set.
- Refuse, with an explicit error, any element or option that has not been tested in eReaxFF mode.

## Wording for any write-up

State what this is: a demonstration that eReaxFF path-integral MD captures barrier penetration and its mass dependence for the eReaxFF pseudo-electron. It is not a prediction for real electrons (the mass is a chosen parameter, electrons are treated as distinguishable, and the result is an equilibrium distribution, not a current). A leakage rate for the silicon work will need a centroid free-energy barrier or an instanton calculation on top of this.

## Deliverable

A short report in `project/docs/` with: pass/fail and numbers for A1–A7, the A5 landscape, the run table with what was actually run, the figures, and anything that did not behave as expected.
