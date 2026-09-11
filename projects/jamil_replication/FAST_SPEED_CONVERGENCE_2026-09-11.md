# Fast PIMD eReaxFF Speed And Bead-Convergence Check

Date: 2026-09-11

System: Jamil 112C/113El graphene case, fast JAX backend, corrected mixed MD integrator, `dt_fs = 0.05`, 10 MD steps, 41 force stages per case.

Run roots:

- `build/jamil_replication/fast_vs_oracle_timing_20260911_134123`
- `build/jamil_replication/fast_vs_oracle_timing_20260911_134733`
- `build/jamil_replication/fast_stability_20260911_122719`
- `build/jamil_replication/fast_stability_20260911_124426`

## Speed

| P | wall_s | s/step | s/force stage | bead-force cases | s/bead-force case | status |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 4 | 75.872 | 7.587 | 1.851 | 164 | 0.4626 | PASS |
| 8 | 80.233 | 8.023 | 1.957 | 328 | 0.2446 | PASS |
| 16 | 85.283 | 8.528 | 2.080 | 656 | 0.1300 | PASS |
| 32 | 97.420 | 9.742 | 2.376 | 1312 | 0.0743 | PASS |
| 64 | 116.989 | 11.699 | 2.853 | 2624 | 0.0446 | PASS |

Interpretation:

- Current P=16 speed is about `8.5 s/step` for the 10-step diagnostic and about `7.8-9.1 s/step` for 100-step diagnostics, depending on warmup/cache and thermostat variant.
- Scaling is sublinear in bead count because all beads are batched through the JAX worker. P=16 is only `1.12x` slower than P=4 and `1.37x` faster than P=64.
- The bottleneck remains one JAX force-stage call, about `2.1 s` at P=16.

## Short-Run P Convergence

| P | physical electron T_K | carbon T_K | physical energy eV | bead RMS spread A | conserved span eV |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 356.889 | 302.696 | -510.639526 | 0.025466 | 19.898 |
| 8 | 314.604 | 303.005 | -513.631205 | 0.035666 | 13.946 |
| 16 | 283.339 | 302.349 | -507.673698 | 0.050581 | 88.635 |
| 32 | 258.665 | 302.745 | -507.950111 | 0.058196 | 1202.758 |
| 64 | 279.750 | 303.060 | -510.640537 | 0.067932 | 10157.514 |

Interpretation:

- Carbon temperature is stable and P-independent in this short diagnostic.
- P=16 electron physical temperature is close to P=64 in this specific 10-step seeded run: `+3.59 K`.
- Bead RMS spread is not fully converged at P=16. P=16 is `0.01735 A` below P=64 and `0.00762 A` below P=32.
- The physical energy and conserved-like span from these short stochastic trajectories are not clean bead-convergence observables. They are useful regression diagnostics, but not enough to certify equilibrium convergence.

## Current Assessment

P=16 is a reasonable near-term production/debug setting because it is much faster than P=32/P64 and gives stable carbon/electron temperatures in the short fast diagnostic. It is not yet proven converged for structural/path observables because bead spread continues to increase from P=16 to P=32 and P=64.

Next convergence test should be a fixed-cost production diagnostic:

- run P=16 and P=32 with the same `dt_fs`, thermostat settings, and wall-clock budget;
- discard an initial equilibration segment;
- compare block averages of bead spread, electron centroid distribution, electron-C radial/localization metrics, energy terms, and carbon temperature;
- accept P=16 only if P16/P32 block-average differences are below predefined tolerances.

## 100-Step P16/P32 Follow-Up

Matched fast-backend runs with `dt_fs = 0.05`, 100 steps, corrected mixed integrator, and carbon Berendsen `tau_fs = 10.0`.

Run roots:

- P16: `build/jamil_replication/fast_stability_20260911_124426/fast_p16_steps100_dt0p05`
- P32: `build/jamil_replication/fast_stability_20260911_153202/fast_p32_steps100_dt0p05`

Speed:

| P | wall_s | s/step | s/force stage | final electron T_K | final carbon T_K | final bead spread A | status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 16 | 778.890 | 7.789 | 1.942 | 251.527 | 1016.449 | 0.148863 | PASS |
| 32 | 981.945 | 9.819 | 2.449 | 283.283 | 1034.006 | 0.181736 | PASS |

Last-half block averages, steps 50-100:

| P | electron T_K mean | electron T_K sd | carbon T_K mean | carbon T_K sd | bead spread A mean | bead spread A sd | physical energy eV mean | max force eV/A mean |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 273.598 | 47.180 | 781.029 | 142.692 | 0.159614 | 0.005843 | -526.683 | 0.304647 |
| 32 | 296.437 | 41.178 | 787.576 | 146.745 | 0.160096 | 0.009295 | -513.398 | 0.155825 |

Interpretation:

- P32 is only `1.26x` slower than P16 for the 100-step run because bead batching remains effective.
- Bead spread is effectively converged between P16 and P32 on the last-half block average: difference is `0.00048 A`, much smaller than the block fluctuations.
- Carbon temperature behavior is also effectively identical between P16 and P32 for this diagnostic.
- Electron physical temperature differs by about `23 K` in the last-half block average, which is within one stochastic block standard deviation for these short trajectories.
- Force/energy diagnostics are not directly converged between P16 and P32; the P32 max force is lower and the physical energy offset differs. Treat these as trajectory-dependent short-run diagnostics, not equilibrium free-energy convergence.

Current conclusion: for this 100-step diagnostic, P16 is good enough for fast production/debug MD against P32 for temperature and bead-spread observables. A publishable convergence claim still needs longer block averages of electron localization/radial observables and energy terms.
