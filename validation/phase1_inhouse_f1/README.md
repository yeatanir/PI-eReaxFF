# Phase 1 / task 1.4 — F1 fix in the in-house engine copy

Repo: `forks/inhouse-engine` (copy of `ereaxff_pimd_modular` 7ff9d90; the original is untouched). Status: **F1 fix done and verified**; F3 experiment documented (not committed); i-PI cross-check see below.

## Plan
1. Add a test hook so the *production* driver (`ereaxff_mixed_baoab_main`: bead_force_manager + mixed_baoab_quantum_step + classical step + thermostats)
   can run a coupled classical/ring-polymer harmonic potential V = 1/2 K|X|^2 + 1/2 k|x-X|^2 (config `force_provider.type: coupled_harmonic`,
   parameters from env `EREAXFF_HARM_K`, `EREAXFF_HARM_k` in eV/A^2). No oracle needed.
2. Test harness (`tests/test_36_production_path_coupled_harmonic/`) compares ring R_g^2, centroid-to-classical distance and classical <X^2> with the exact
   finite-P values at P = 4, 8, 16. Commit it **before** the fix: it must fail.
3. Fix: remove the 1/P on quantum bead kicks (`half_kick_from_summary`), keep mean force for classical atoms. Test must pass.

## Results (updated as work proceeds)

Repo `forks/inhouse-engine` commits: `418e6cf` import (unchanged) · `7db9bc4` test hook + test_36 (expected to fail) · `a80d240` F1 fix · `1e0bdb1` whole-system conserved_like.
Test: `make test_36_production_path_coupled_harmonic` (needs only `build/ereaxff_mixed_baoab_main`, python3+numpy; no oracle).
Setup: classical C (12.011 amu, Langevin, gamma 0.05/fs) bound by K = 1.0 eV/A^2 to the origin; ring-polymer "El" (1.0 amu, PILE-L, centroid gamma 0.05/fs) bound to it by k = 0.26 eV/A^2; T = 300 K; dt = 0.5 fs;
6 seeds × 2e6 steps per P, first 10 % discarded, standard error over seeds. Pass criterion: within 3 SE **and** within 3 %.
Exact values: <X^2> = 3kT/K, <|x_c-X|^2> = 3kT/k, <R_g^2> = 3kT sum_{l>=1} 1/(m w_l^2 + k), w_l = 2 (P kT/hbar) sin(pi l/P) (standard_ipi).

### Before the fix (commit `7db9bc4`; full log: `before_fix_full.log`) — FAILS

| P | <X^2> deviation | <\|x_c-X\|^2> deviation | <R_g^2> deviation |
|---|---:|---:|---:|
| 4 | +46.8 % (+64 SE) | +253 % (+209 SE) | +3.3 % (+53 SE) |
| 8 | +125.8 % (+194 SE) | +577 % (+204 SE) | +2.7 % (+66 SE) |
| 16 | +287.5 % (+180 SE) | +1205 % (+205 SE) | +2.6 % (+73 SE) |

The error grows with P, as expected when the ring's bead potential is weighted P times too weakly: the electron wanders far from the classical atom.

### After the fix (commit `a80d240`; full log: `after_fix_full.log`) — PASSES

| P | <X^2> deviation | <\|x_c-X\|^2> deviation | <R_g^2> deviation |
|---|---:|---:|---:|
| 4 | +0.06 % (+0.2 SE) | +0.05 % (+0.3 SE) | -0.05 % (-0.8 SE) |
| 8 | -0.32 % (-1.0 SE) | -0.25 % (-2.0 SE) | -0.01 % (-0.2 SE) |
| 16 | +0.04 % (+0.2 SE) | -0.15 % (-0.4 SE) | +0.01 % (+0.3 SE) |

All within 3 SE and 3 % (largest: -0.32 %). The centroid distance at P = 8 sits at -2.0 SE; no other value exceeds 1 SE. A small negative bias of a few tenths of a percent is expected from the O(dt^2) integration error at dt = 0.5 fs (not separately measured).

The change itself: one line in `half_kick_from_summary` (`force_rescale = P_max` instead of `P_max / P_active`), i.e. quantum beads get the unscaled F_j, classical atoms keep the mean force (the driver sums the manager's F_j/P_max over beads). Thermostat temperature P·T, omega_P = P kT/hbar and the velocity initialisation were already in the i-PI convention. With mixed bead counts (P_active < P_max for some particles) the classical mean force is still taken over P_max replicas with centroid-filled slots; that case is not covered by this test.

Other existing tests rebuilt and run on the fixed code (no oracle needed): test_00, 01, 02, 03, 04, 05, 06, 07, 09, 11 all PASS.

### Whole-system conserved quantity and a second defect (F3)
NVE (no thermostats), P = 8, 5 ps, new conserved_like (ring KE + spring + P·KE_classical + ΣE_j):

| dt (fs) | change of conserved_like, F1 fix only | with F3 experiment (electron first kick reuses start-of-step forces) |
|---|---:|---:|
| 0.5 | +0.184 eV (carbon heats to ~10^4 K over 50 ps) | +3.6e-4 eV |
| 0.25 | +0.073 eV | — |
| 0.1 | +0.025 eV | +1.4e-5 eV |

So after the F1 fix the NVE drift that remains is the audit's **F3** (electron pre-kick uses forces at a mixed configuration: new carbons, old electron). Reusing the cached start-of-step forces removes it (500× smaller at dt 0.5, ~dt^2 scaling).
That change is **not committed** (outside this task): `f3_reuse_forces_experiment.patch` (19 lines; it also saves one force batch per step). It also confirms that the new conserved_like formula is the right invariant.
Sampling tests pass anyway with thermostats on at dt = 0.5 fs, because the thermostats remove the heating; the NVT bias from F3 is small but not zero.

## Cross-check against i-PI `classical_atoms` (same physical potential, P = 8)
Parameters converted to atomic units (K = 1.0, k = 0.26 eV/A^2; dt 0.5 fs; centroid tau 20 fs; PILE-L for both atoms; quantum atom = i-PI "H", 1.00794 amu), 4 seeds × 60 000 steps
(`ipi_crosscheck/`): <X^2> +2.7 % (+1.0 SE), <|x_c-X|^2> +1.7 % (+1.1 SE), <R_g^2> +0.12 % (+0.2 SE) against the exact finite-P values.
The i-PI errors are larger (SE 2.5 % / 1.6 % / 0.5 %, fewer samples), so this is a consistency check, not a precision comparison. Both codes agree with the analytic values; the fixed in-house path and the i-PI patch therefore describe the same ensemble.
(Reproduce: `python3 run_seed.py <seed> 60000` ×4, then `python3 analyze.py`.)

## Limits and not done
- Single potential (two harmonic couplings), one electron, uniform P, dt = 0.5 fs, one temperature. Not a test of anharmonic forces or of the oracle path.
- The production driver paths that need the legacy oracle (test_08, 12-35 except those using `type: dummy`) were not run (the oracle is gone). The fix touches only the shared integrator, which those tests also use; their expected values were produced with the old scaling and would need regeneration.
- `conserved_like` assumes uniform P and is only verified in NVE (see F3 table), not with thermostats.
- F3 is not fixed in the repo. It is the main remaining integrator defect (NVE heating, first order in dt).
- The old `harmonic_pimd_main` path was already correct (scale 1.0) and is untouched.
