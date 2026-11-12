# phase4_frames_ag — exact ereac reference frames, Ag/El (published force field), WITH charges

Force field: Evangelisti, Fichthorn, van Duin, J. Chem. Phys. 153, 104106 (2020), SI (`ffield_Ag_El`, general parameter 38 = 5.7 A electron taper; Av = Ag+, El).
System: neutral clusters = Ag+ (type Av, charge +1) + El (-1) per atom, non-periodic, fixed charges (`icharg 8`); El at the cation +0.001 A. Geometries from the SI (Ag3, Ag4, Ag5_2, Ag6_2, Ag7_2, Ag7_1).
Source runs: rebuilt ereac (export build), NVE, 300 K random velocities (`imdmet 3, iravel 1, mdtemp 300`), dt 0.25 fs. Etot span: Ag3 0.02, Ag4 0.03, Ag5_2 0.08, Ag6_2 0.16, Ag7_2 0.19 kcal/mol over 2 ps; Ag7_1 0.19 over 5 ps. Coulomb at step 0: Ag3 -365.58, Ag7_1 -872.91 kcal/mol.
Per case: `frames.json` (20 frames; Ag7_1: 50; positions, 14 terms, ESTRC, total gradient; kcal/mol, kcal/mol/A), `control`, `ffield`, `geo`.
History: the first version of this set (kept as `phase4_frames_ag_INVALID_zero_charges`) had `END` before the `MOLCHARGE` lines, so the legacy reader assigned zero charges and the Coulomb term was 0 (found by the LAMMPS-adapter fork). The generator now writes `END` last.
