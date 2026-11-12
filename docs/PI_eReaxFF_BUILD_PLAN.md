# PI-eReaxFF — Build Plan

Date: 2026-11-11. Brief: `FABLE_PLANNING_PROMPT.md`. Primary input: `PI_eReaxFF_AUDIT_REPORT.md`.

Sources read for this plan (nothing was executed from them, no simulation was run):

| Source | Revision |
|---|---|
| In-house engine `ereaxff_pimd_modular/` | git `7ff9d90` |
| Legacy eReaxFF `…/Margaret/eReaxFF_src/` (`reac.f`, `poten.f`) | files dated 2026-04-20/23 (see V-14) |
| Metal-Reax | `735ad6a` (2026-10-09) |
| i-PI | `93162d2` (2026-10-24) |
| LAMMPS `stable` branch | `098ffde` (2026-10-29), `LAMMPS_VERSION "30 Sep 2026"`. Metal-Reax pins `8de817dd`; I did not diff the two. |

Evidence labels: **[src]** read in source at the cited line; **[data]** read in a saved output file; **[derived]** my own derivation or estimate; **[unverified]** not checked, with a note on how to check.

---

## 1. Decision summary

**Recommended architecture.** i-PI drives the path integral. Forces come from eReaxFF implemented in the Metal-Reax CPU-64 engine, running as a LAMMPS pair style, with one LAMMPS client per bead connected through `fix ipi`. Classical carbons with ring-polymer electrons are obtained with a small per-atom patch to i-PI. Until the Metal-Reax port is validated, the existing JAX kernel is wrapped as an i-PI driver so that path-integral validation does not wait for the engine.

**Key decisions**

1. **i-PI is the integrator; the in-house integrator is retired from production.** i-PI already has the consistent P·T convention, a one-force-evaluation BAOAB/OBABO step, PILE-L, PILE-G, SVR, a Cayley propagator, a whole-system conserved quantity and centroid-virial estimators (V-1 to V-4). These are exactly findings F1 to F5 of the audit. Runner-up: fix the in-house engine. I would switch only if the i-PI patch (decision 2) is rejected in review *and* cannot be carried as a fork.
2. **Mixed classical/quantum by projecting out internal modes in i-PI.** A classical atom is a ring polymer whose internal normal modes are held at zero. In i-PI's convention this gives physical mass, mean bead force and temperature T for that atom exactly (Appendix C). i-PI calls `pconstraints()` after every momentum update and thermostat step, so the patch is one block there plus input plumbing and one property fix (V-5, V-6; Appendix B). Runner-up: all atoms as ring polymers, which needs no patch and costs nothing extra in force calls (§3.4). Run it anyway as a comparison.
3. **Force engine: extend Metal-Reax CPU-64, not the legacy Fortran.** The legacy V2 oracle cannot be rebuilt, and Metal-Reax is a tested engine under the owner's control. The eReaxFF additions are larger than the audit describes: the over-coordination term carries an extra electron factor with two more parameters (V-11). Runner-up: make the JAX kernel the production engine. Switch if the Metal-Reax port exceeds roughly twice the upper effort estimate, or if Phase 0 shows that LAMMPS-flavour ReaxFF cannot be brought to parity with the standalone code.
4. **Parallelism by processes first; no batching API yet.** One LAMMPS client per bead over Unix sockets uses only existing, validated code and keeps per-bead solver history in each client for free. The multi-geometry batch API is specified (§3.3) but built only when GPU work starts.
5. **Reference data must be re-established before any port is validated.** Three facts found in this review undermine the current references (V-12, V-13, V-14): the saved V2 outputs are central-cell (no periodic images), the Jamil-target ffield has a zero electron taper radius in the slot the base code reads, and the Margaret source carries local 2026 edits made after the shipped binaries were built. Phase 0 settles these.

**Rejected**

| Option | Why |
|---|---|
| LAMMPS `fix pimd/langevin` for the mixed system | Only partly group-aware (V-7). Needs a new cross-partition force-averaging fix plus synchronised thermostat noise, and still lacks PILE-G. More new code than the i-PI patch, with fewer validated parts. Kept as a cross-check for the all-ring-polymer case. |
| Physical electron mass | Needs P of order 10³ and double-counts quantum effects already folded into a classically fitted potential (§3.5). |
| GPU or batching first | At 225 atoms the README of Metal-Reax itself reports no GPU benefit below about 1000 atoms. |
| Wrapping the V2 binaries as a long-term engine | Cannot be rebuilt, file-based IPC, seconds per call. |

---

## 2. Verified facts

| # | Claim | How verified | Confidence |
|---|---|---|---|
| V-1 | i-PI ring frequency ω_P = P·k_B·T/ħ | `ipi/engine/normalmodes.py:455-462` [src] | High |
| V-2 | i-PI thermostat temperature is P·T | `ipi/engine/motion/dynamics.py:235-238`, piped to the thermostat at `:169` [src] | High |
| V-3 | i-PI kicks beads with the unscaled bead force | `dynamics.py:427-443` (`beads.p += forces.f * pdt`) [src] | High |
| V-4 | i-PI conserved quantity = spring + kinetic + Σ bead potentials + thermostat/constraint terms | `ipi/engine/ensembles.py:254-266` [src] | High |
| V-5 | i-PI has one global `nbeads`. Per-atom controls that exist: `fixatoms_dof` (all modes of an atom frozen), `open_paths` and `bosons` (per-atom lists inside the normal-mode object). No per-atom "classical" flag. | `ipi/inputs/beads.py:27-46`; `dynamics.py:116-130,303-306`; `normalmodes.py:87-112,766,793` [src]. Absence is from grep, so medium. | Medium-high |
| V-6 | `pconstraints()` runs after every `pstep` and after every thermostat step in NVE/NVT, and already zeroes normal-mode momenta of fixed atoms with the removed kinetic energy booked into `eens` | `dynamics.py:368-401` (constraint), `:453-513` (MTS loops), `:542-562` (NVT step) [src] | High |
| V-7 | LAMMPS `fix pimd/nvt` never uses the group mask. `fix pimd/langevin` masks kicks, thermostat and estimators but gathers and normal-mode-transforms all atoms. | `grep -c groupbit`: 0 in `fix_pimd_nvt.cpp`, 18 in `fix_pimd_langevin.cpp`; masked `b_step` `:772-791`; unmasked `nmpimd_transform` `:1197-1226` and `inter_replica_comm` `:1361-1390`; `ntotal = atom->natoms` `:344` [src]. Not run. | High (source), behaviour untested |
| V-8 | F1 is real: the in-house engine uses ω_P = P·kT/ħ and thermostat temperature P·T, but divides bead forces by P | `src/pimd/normal_modes.f90:75`; `src/drivers/ereaxff_mixed_baoab_main.f90:953` and `:634`; `src/force/bead_force_manager.f90:241`; `src/integrator/mixed_baoab_integrator.f90:163-164` (the rescale restores F/P_active, not F) [src] | High |
| V-9 | The harmonic test bypasses the production force path | `physical_force_scale=1.0` only in `src/drivers/harmonic_pimd_main.f90:113` and `src/integrator/baoab_integrator.f90` [src] | High |
| V-10 | `dista2` **does** apply minimum image: `anint` for orthogonal cells, nearest of 27 images otherwise | `reac.f:5464-5495` [src] | High |
| V-11 | eReaxFF couples electrons into bonded terms through three channels: `avalc`, `stlpc`, and a factor `exp(-pxel2·eshare2·Σ pxel1·BO/abo)` in over-coordination. `pxel2` is atom line 2 column 5, `pxel1` bond line 2 column 8, `vincr` atom line 3 column 3, `alphag` line 3 column 8, `betael` line 4 column 3. `znucel` is a hard-coded element table, not an ffield entry. | `poten.f:423-432` (lone pair), `:887-949` (over-coordination), `reac.f:4514-4524,4670-4671` (parser), `reac.f:4590-4615` (`znucel`), `reac.f:8779-8851` (`elecpos`) [src] | High |
| V-12 | The saved V2 reference data behave as central-cell (no periodic images) for nonbonded, bond-order, electron-pair and electron-sharing terms | Comments recording what the JAX port had to do to match the oracle: `fast_backend/jax_force_kernel.py:464-468, 525-527, 1249-1252, 1516-1520` [src]. Indirect: I did not recompute the oracle energies. | Medium |
| V-13 | In the Jamil-target ffield, general parameter 38 is `0.0000 !Molecular energy (not used)`. The base code sets the electron taper radius `swb = vpar(38)` and skips every pair with `r > swb`. V2 evidently used 10 Å instead. Genuine eReaxFF ffields on disk have 3.05 and 5.70 labelled "Electron taper". | `build/jamil_replication/finite_difference_oracle_force/central/ffield` line 40; `poten.f:3949,3972`; `…/central/tap7th_calls.out` (enter with 0.0, 10.0); `…/central/single_point_energy_terms.out` (electron = +16563.55 kcal/mol); `Margaret/eReaxFF_Jan3_2024/ffield` and `eReaxFF_src/test_co2/ffield` line 40 [src, data] | High for the facts; the cause in V2 is unknown |
| V-14 | The Margaret source is not pristine 2016 code. `chargesa` contains a "compact C-only ACKS2" block whose comment cites the 60El+60C system; `reac.f`/`poten.f` are dated 2026-04-20/23, the `ereac*` binaries 2026-04-10. | `reac.f:4041-4130`; `ls -la`; `compile.log` [src] | High |
| V-15 | The base `electron` routine sums over periodic images through the Verlet list (at least ±1 image per direction), and excludes self-image pairs | `reac.f:5584-5650` (`kxt = max(kx,1)`), `poten.f:3954-3966` [src] | High |
| V-16 | `icharg 8` sets charges from `MOLCHARGE` and returns without a solve. The Jamil control uses `icharg 8`, `tstep 0.25`, `tdamp1 100`. | `reac.f:2971 ff.`; `…/central/control` [src, data] | High |
| V-17 | Electron-pair formulas: core–El `-332.0638·betael·znucel·erf(√(2·alphag)·r)/r`, El–El `+332.0638·betael_El·erf(√alphag_El·r)/r`, both times the 7th-order taper. The El–El term uses the parameters of atom `i1` only, not a geometric mean. | `poten.f:3987-3992`, `:4119-4127` [src] | High |
| V-18 | Valence-angle lone-pair adjustment uses the **uncorrected** `stlp` (`vlpnc`), and bond-order correction does not use `avalc` | `poten.f:2351-2361`; `avalc` occurs only in `lonpar`/`ovcor` (grep) [src] | High |
| V-19 | Legacy heavy-atom switch is `amas > 31.0`; LAMMPS and Metal-Reax use `mass > 21` | `poten.f:887`; Metal-Reax `compat_flags.hpp` (`heavy_atom_terms`) [src] | High |
| V-20 | Metal-Reax keys a triple-bond compatibility flag on `gp[37] == 2`, which is the same slot as the eReaxFF electron taper | Metal-Reax `docs/FEATURE_MATRIX.md` §3 (`compat.triple_bond_stabilisation`) [src] | High |
| V-21 | i-PI offers batched requests (`batch_size`), a shared-memory transport (`mode='shm'`), an MPI interface and client matching modes `auto`/`any`/`lock`. Batched mode drops per-replica client matching. | `ipi/inputs/forcefields.py:245-274`; `ipi/interfaces/sockets.py:761-905, 914-1000, 1483-1487, 1753-1765`; `ipi/interfaces/mpi.py:1-30`; `drivers/f90/driver.f90:45-60, 596-640` [src] | High |
| V-22 | LAMMPS `fix ipi` warns for any group other than `all`, forces a neighbour check every step and wraps atoms into the box | `MISC/fix_ipi.cpp:185, 275-276, 383-388` [src] | High |
| V-23 | i-PI is offered under GPL or MIT at the user's choice | `licenses/LICENSE.md` [src]. GPL version not checked. | High |
| V-24 | Metal-Reax implements all ReaxFF terms, fixed charges (`checkqeq no`), QEq and ACKS2 on CPU-64 and validates against stock LAMMPS, not against the standalone Fortran code. Batching is "Planned (M8)". | `docs/FEATURE_MATRIX.md` §1, §4, §5, §9; `README.md` [src] | High |
| V-25 | Local machine: Apple M5 Max, 18 cores, 64 GB; `gfortran` 15.2 and `mpirun` present; no `lmp` on PATH; `ereac` is arm64 | `sysctl`, `which`, `file` | High |

Audit claims I did **not** re-verify: the detailed step ordering behind F3, and F4, F5, F7, F8, F9. They become moot for production under decision 1, so I did not spend time on them.

---

## 3. Target architecture

### 3.1 Components

```
                    +-----------------------------------------------+
                    |  i-PI  (fork branch: classical-atoms)         |
                    |  beads, normal modes, BAOAB, PILE-L/G, SVR,   |
                    |  estimators, restart, RNG                     |
                    +---------------------+-------------------------+
                                          | unix socket (or shm), one request per bead
              +---------------------------+---------------------------+
              |                           |                           |
      +-------v-------+           +-------v-------+           +-------v-------+
      | LAMMPS client |   ...     | LAMMPS client |           | interim:      |
      | fix ipi       |           | fix ipi       |           | JAX driver    |
      | pair reaxff/  |           | pair reaxff/  |           | (Python,      |
      |  metal        |           |  metal        |           |  batched)     |
      |  ereaxff yes  |           |  ereaxff yes  |           +---------------+
      +-------+-------+           +---------------+
              |
      Metal-Reax core (CPU-64): ffield parser, bond orders, ReaxFF terms,
      + electron pair term, electron sharing, modified lone pair / over-under
```

Classical eReaxFF MD (the Jamil benchmark) is the same LAMMPS input without `fix ipi`, using `fix nve` and `fix temp/berendsen`.

### 3.2 Data flow per MD step (BAOAB, one force evaluation)

1. i-PI: half kick with stored forces, then project internal modes of classical atoms (B).
2. i-PI: exact free-ring propagation for dt/2 (A).
3. i-PI: thermostat for dt, then project (O).
4. i-PI: free-ring propagation for dt/2 (A).
5. i-PI sends P geometries. Each client wraps positions, updates ghosts, evaluates energy, forces and virial, and replies.
6. i-PI: half kick, then project (B). Properties and checkpoints are written.

**Where state lives.** All dynamical state (bead positions and momenta, thermostat state, RNG, conserved-quantity ledger) is in i-PI and its checkpoint. A client holds only neighbour lists and, when a charge solver is active, its own solver history. With fixed charges a client is stateless, so a crashed client can be restarted without affecting the trajectory.

### 3.3 Interfaces

**i-PI to LAMMPS (now).** Standard i-PI socket protocol, `<ffsocket mode='unix' matching='auto'>`. With `matching='lock'` a bead always goes to the same client, which is what a per-bead charge-solver warm start needs.

**Batch API (specified now, built in Phase 6).** A standalone C++ driver in Metal-Reax, because a LAMMPS instance holds one system and cannot evaluate a batch.

```cpp
struct BatchInput {
  int nbatch, natoms;
  const double* x;        // [nbatch][natoms][3], Angstrom, unwrapped allowed
  const double* cell;     // [nbatch][9]
  const int*    type;     // [natoms], shared by all structures
  const double* q_fixed;  // [natoms] or nullptr when a solver is used
};
struct BatchOutput {
  double* energy;         // [nbatch]
  double* force;          // [nbatch][natoms][3]
  double* virial;         // [nbatch][9]
  int*    status;         // [nbatch]: OK, NOT_CONVERGED, NONFINITE, LIST_OVERFLOW
};
class BatchEngine {       // owns per-structure neighbour lists and solver history, indexed by batch slot
 public:
  void evaluate(const BatchInput&, BatchOutput&);
};
```

- CPU: a parallel loop over structures calling the existing single-geometry kernels with the standalone image expander. Requires that `compute_bonded_core` and `compute_nonbonded_core` hold no shared mutable state [unverified: check before relying on it].
- GPU: concatenate the P structures as disjoint subsystems of one P·N-atom problem, with the replica index added to the cell key so that no cross-replica neighbours form. At P = 32 to 64 this is 7 200 to 14 400 atoms, the size range where Metal-Reax reports a GPU benefit.
- Transport: `<ffsocket mode='shm' batch_size='P'>` (V-21). Batched mode drops client matching, so slot-to-bead identity must be shown stable by a test before any per-slot solver history is trusted [unverified].
- Failure handling: any non-OK status fails the whole batch. The driver exits non-zero and i-PI stops at its last checkpoint. No force is ever returned from a non-converged or non-finite evaluation.

### 3.4 Mixed system versus all-ring-polymer (questions 1 to 3)

**i-PI.** No existing mechanism does this (V-5). The smallest correct change is the projection in Appendix B: about 40 lines of logic in two files, plus input plumbing and tests. Thermostats, the normal-mode propagator, forces, ring-polymer contraction and MTS need no change, because a mode that starts at zero position and zero momentum stays there under the free-ring propagator, and contraction of identical beads returns identical beads. Not supported in the first version, and to be refused with an explicit error: NPT/NST, Suzuki-Chin, constrained-centroid, GLE thermostats, the `bab` propagator, and classical atoms that are also listed in `bosons` or `open_paths`. Effort 4 to 8 person-days, low risk. Maintenance: carry as a fork branch pinned to an i-PI commit, and open an upstream pull request once Phase 1 tests pass. A subclass or plugin is not simpler, because i-PI's input layer maps mode names to classes and would need the same plumbing.

**LAMMPS.** V-7 confirms the audit. A correct mixed scheme would need (a) a new fix that replaces the force on each classical atom by the average over partitions, (b) identical integration and identical thermostat noise for those atoms in every partition, and (c) periodic re-synchronisation of their positions against round-off divergence. I estimate 6 to 10 days with medium risk, and the result still has only the thermostats `fix pimd/langevin` offers. The i-PI route is smaller and better tested.

**All atoms as ring polymers.** Force cost is unchanged, since every bead geometry is a full-system evaluation in either scheme. The integrator cost of 112 extra ring polymers is negligible. Physics: carbon nuclear quantum effects are real at 300 K (C–C stretch near 1600 cm⁻¹ gives βħω ≈ 7.7), and the bead count is set by the electrons anyway (§3.5). The drawbacks are scientific rather than computational: the comparison with the classical Jamil benchmark then changes two things at once, and ReaxFF carbon parameters are themselves fitted for classical nuclei. Use it as a no-patch smoke test in Phase 1 and as a reported comparison, not as the primary result.

### 3.5 Electron model (question 8)

Numbers below are estimates from the electron–core term alone with the carbon parameters of the Jamil-target ffield (`betael` 0.0027, `alphag` 22.0, `znucel` 6), harmonic at the well bottom. They assume that term is active, which depends on the taper question in V-13. All are [derived].

| Quantity | 1 amu | Physical electron mass |
|---|---|---|
| Well depth of the electron–core term | ≈ 40 kcal/mol | same |
| Curvature | ≈ 51 eV/Å² | same |
| ω | ≈ 0.70 fs⁻¹ (≈ 3700 cm⁻¹) | ≈ 30 fs⁻¹ |
| βħω at 300 K | ≈ 18 | ≈ 770 |
| Zero-point energy 3ħω/2 | ≈ 0.69 eV (≈ 16 kcal/mol) | exceeds the well depth |
| Free-ring radius of gyration | 0.20 Å | 8.6 Å |
| P for about 1 % error in a harmonic ⟨x²⟩, from (βħω/P)²/24 | 32 to 64 | above 1000 |

Consequences:

- **Classical equilibrium statistics do not depend on mass, path-integral statistics do.** The 1 amu choice in the Jamil paper is harmless for classical static properties. In a path integral the mass is a physical parameter of the model: it sets the ring size. A PIMD result with 1 amu is a property of "the eReaxFF pseudo-particle quantised with mass 1 amu", not of an electron.
- **The effect is not small.** A zero-point energy near 16 kcal/mol against a 40 kcal/mol well, and a quantum width (≈ 0.07 Å per direction) three times the classical thermal width (≈ 0.02 Å) and comparable to the 0.2 to 0.4 Å range of the sharing function exp(−181 r²), mean path-integral electrons will bind and share differently from classical ones. This is the scientific signal, and also a warning that eReaxFF was fitted with classical electrons, so quantising them double-counts whatever quantum character the fit absorbed.
- **Defensible:** static observables (bead spread, radius of gyration, El–C distributions, host occupancy, centroid-virial kinetic energy) reported as a sensitivity study in P and in mass. **Not defensible:** electron dynamics or hopping rates from RPMD/TRPMD of a fictitious-mass particle, and any claim about physical electronic delocalisation. Exchange is ignored; that is acceptable only while electrons stay localised on separate hosts.

Experiments to pin down P and mass (Phase 5): (E1) classical velocity spectrum of El to measure the real ω_max; (E2) P = 8, 16, 32, 64, 128 on the 60C/60El system; (E3) mass scan 0.25, 1, 4, 16 amu at converged P, checking the approach to the classical distribution at large mass; (E4) compare the measured ring size with the free-particle value 0.20 Å·√(1 − 1/P²) to quantify confinement.

---

## 4. Work breakdown

Effort is in person-days for one developer who knows the codebases. Thresholds marked † are proposals of this plan, not literature standards.

### Phase 0 — Establish ground truth (4 to 7 d). Blocks Phases 3 and 4.

| Task | Purpose and deliverable | Files | Effort | Risk | Acceptance test |
|---|---|---|---|---|---|
| 0.1 Rebuild `ereac` from the Margaret source | Know what the source computes today. Deliverable: build log, binary hash. | `eReaxFF_src/makefile` (copy; do not edit in place) | 0.5 | Low | Reproduces `test_co2/fort.73` final line (Eelect −25.68, Ebond −702.63) to printed precision |
| 0.2 Observation-only export patch | Full-precision energies per term and forces from the rebuilt code. Deliverable: a patch file with hash, in the style of Metal-Reax ADR-019. | new patch against `reac.f` (`encalc` exit) | 1 | Low | Patched and unpatched builds give identical `fort.73` |
| 0.3 Resolve the taper and periodicity questions | Explain V-12 and V-13. Run the 225-atom Jamil geometry through rebuilt `ereac`, periodic and with a large non-periodic cell, with `vpar(38)` 0 and 10. Deliverable: a 2×2 table of 14 energy terms against `single_point_energy_terms.out`. | none | 1 to 2 | Medium: no combination may match | One combination reproduces the saved V2 terms to 1e-4 kcal/mol†. If none does, V2 saved outputs are demoted to "regression data for the JAX port only". |
| 0.4 Plain ReaxFF baseline: standalone versus LAMMPS flavour | Measure how far Metal-Reax CPU-64 is from `ereac` with no electrons present (112 C graphene, a C/H molecule, periodic and not). Deliverable: per-term difference table and a list of causes. Known candidates: V-19, V-20, bond-order and angle cut-offs from `control`, the polarisation energy that LAMMPS always adds. | none | 2 to 3 | Medium | Every per-term difference above 1e-5 kcal/mol/atom† is traced to a named cause |
| 0.5 Owner decisions | See §8 items 1 to 3. | — | — | — | Written answers |

### Phase 1 — Path-integral core in i-PI (5 to 10 d). Independent of Phase 0.

| Task | Purpose and deliverable | Files | Effort | Risk | Acceptance test |
|---|---|---|---|---|---|
| 1.1 `classical_atoms` patch | Appendix B | `ipi/engine/motion/{motion,dynamics}.py`, `ipi/inputs/motion/{motion,dynamics}.py`, `ipi/engine/properties.py` | 2 to 4 | Low | L-P1 to L-P5 in §5 |
| 1.2 Coupled harmonic test potential | Analytic target for a classical particle bound to a ring polymer. Deliverable: a small PES for the i-PI Python driver. | `drivers/py/pes/` (new file) | 0.5 | Low | Forces match finite differences to 1e-8 relative |
| 1.3 Regression tests in i-PI's layout | Deliverable: test directories with reference outputs | `ipi_tests/regression_tests/tests/NVT/…` | 1 to 2 | Low | i-PI test suite passes with the new cases |
| 1.4 F1 fix and production-path harmonic test in the in-house engine | Independent second implementation of the mixed scheme. Remove the 1/P on electron kicks; run the same coupled-harmonic case through `mixed_baoab_quantum_step` and `bead_force_manager`. | `src/force/bead_force_manager.f90`, `src/integrator/mixed_baoab_integrator.f90`, one new test | 1 to 2 | Low | Test fails before the fix and passes after, at P = 4, 8, 16 within 2 %† of the analytic finite-P variance; agrees with i-PI within combined statistical error |
| 1.5 All-ring-polymer smoke run | Proves the i-PI plus driver chain with no patch | inputs only | 0.5 | Low | Runs 1000 steps; conserved quantity criterion L-P4 |

### Phase 2 — Interim engine: JAX kernel as an i-PI driver (2 to 4 d). Needs Phase 1.

| Task | Purpose and deliverable | Files | Effort | Risk | Acceptance test |
|---|---|---|---|---|---|
| 2.1 Python driver wrapping `jax_force_kernel` | Correct eReaxFF PIMD at low speed while Phase 3 runs. Batched over beads. | new `fast_backend/ipi_driver.py`; reuse `jax_live_worker.py` logic | 2 to 3 | Low. Known limits: central-cell only, fixed charges only, about 2 s per batch. | Forces equal the stored oracle forces to 1e-6 eV/Å on the existing parity frames; P = 1 NVE drift criterion L-E1 |
| 2.2 First real diagnosis of F2 | With a correct integrator, is eReaxFF PIMD stable at dt = 0.1 to 0.25 fs? | inputs only | 0.5 to 1 | This is the main early risk probe | 2000 steps at P = 16, dt = 0.1 fs, NVE after equilibration: conserved quantity within L-E2 |

### Phase 3 — eReaxFF in Metal-Reax CPU-64 (17 to 30 d). Needs Phase 0.

Change list in Appendix A. Each task ends with fixtures and a `FEATURE_MATRIX` row, as the repository's `docs_sync` test requires.

| Task | Files | Effort | Risk | Acceptance test |
|---|---|---|---|---|
| 3.1 Parser and data model: `pxel1`, `pxel2`, `vincr`, `alphag`, `betael`, electron taper, element table for `znucel` and sharing sign, `ereaxff` option, exemption of El pairs from the strict missing-pair rule, disabling the `gp[37]==2` predicate in eReaxFF mode | `include/reaxmetal/{forcefield,compat_flags}.hpp`, `src/io/ffield_parser.cpp`, `src/core/compat_flags.cpp` | 2 to 3 | Low | Table dump equals values read by `reac.f` for three ffields (Jamil target, `test_co2`, Jan 2024) |
| 3.2 Electron pair term (core–El, El–El) with its own taper | `include/reaxmetal/terms.hpp`, `src/cpu/nonbonded.cpp`, `include/reaxmetal/energy_terms.hpp` | 2 to 3 | Low | L-F1, L-F2 for the `Electron` term |
| 3.3 Electron sharing and modified lone pair | `terms.hpp` (`atom_quantities`, `lone_pair`), `src/cpu/bonded.cpp` (new pre-pass before the atom-quantities block, lone-pair block) | 3 to 5 | Medium | L-F1, L-F2 for `LonePair` |
| 3.4 Modified over/under-coordination with all derivative channels | `terms.hpp` (`over_under`), `src/cpu/bonded.cpp` (over/under block and the bond-derivative gather) | 5 to 10 | **High**: about 800 lines of legacy derivative code (`poten.f:862-1700`), neighbour-of-neighbour sums, ghost-shell requirements | L-F1, L-F2 for `Over`/`Under`; ghost-shell check fails loudly when the shell is too thin |
| 3.5 LAMMPS adapter: keyword `ereaxff yes`, 15th energy slot, fixed-charge polarisation switch, refusal of `backend metal` and of hole particles | `plugin/adapter/pair_reaxff_metal.{h,cpp}`, `src/io/pair_settings.cpp`, `src/core/capabilities.cpp` | 2 to 3 | Low | In-LAMMPS energy and forces equal the standalone harness to 1e-12 relative |
| 3.6 Reference fixture set (§5.1) and frozen tolerances | `tests/fixtures/`, `tolerances/` | 3 to 6 | Medium: depends on Phase 0 outcome | L-F3 |

### Phase 4 — Classical eReaxFF validation (5 to 10 d). Needs Phase 3.

| Task | Effort | Acceptance test |
|---|---|---|
| 4.1 NVE conservation, 112C/113El and 60C/60El, dt = 0.25 and 0.1 fs | 1 to 2 | L-E1 |
| 4.2 Trajectory-level comparison with `ereac` (same start, NVE, 200 steps) | 1 to 2 | Positions agree to 1e-4 Å† at step 200, or divergence is traced to a documented discontinuity |
| 4.3 Jamil benchmark: Berendsen 300 K, τ = 100 fs, dt = 0.25 fs, 1 amu electrons | 2 to 4 | Mean temperature within 2 K†; electron-host statistics consistent with the paper's figure within block error. The quantitative target needs the owner's reading of the paper (§8 item 4). |
| 4.4 P = 1 through i-PI equals LAMMPS classical | 1 to 2 | L-P0 |

### Phase 5 — PIMD validation and production protocol (8 to 15 d). Needs Phases 1 and 4 (or 1 and 2 for an early pass).

| Task | Effort | Acceptance test |
|---|---|---|
| 5.1 Ring equilibration protocol: start from collapsed rings, PILE-L with short centroid τ, monitor spread until stationary | 1 | Bead spread block means stationary within error over the last half |
| 5.2 Thermostat set: PILE-L for equilibration, PILE-G (SVR on centroids) for production; both must agree on static observables | 1 to 2 | L-S2 |
| 5.3 Timestep study at P = 32: dt = 0.05, 0.1, 0.25 fs; exact and Cayley propagators | 1 to 2 | L-E2, and static observables independent of dt within error |
| 5.4 Experiments E1 to E4 of §3.5 | 3 to 6 | L-S3, L-S4 |
| 5.5 Mixed versus all-ring-polymer versus classical comparison | 1 to 2 | Reported, no pass/fail |
| 5.6 Restart exactness and reproducibility | 1 | L-R1 |

### Phase 6 — Performance (optional; 10 to 25 d). Needs Phase 4.

Batch driver on CPU, then the disjoint-replica GPU path (§3.3). Start only if the Phase 4 timings and the benchmarks listed in §9.1 show that the process-per-bead layout is the bottleneck for the target system size.

### Phase 7 — Charges that respond to electrons (optional; 8 to 15 d).

Compact carbon-only ACKS2 with electron charges frozen and their field folded into the right-hand side. Whether LAMMPS `fix acks2/reaxff` on a carbon subgroup already includes the field of out-of-group charges is [unverified]; read `fix_acks2_reaxff.cpp` and test against rebuilt `ereac` with `icharg 7` before designing anything. Per-bead history comes from `matching='lock'`. Add a residual-versus-force-error test before trusting any tolerance.

---

## 5. Validation ladder

### 5.1 Reference data

Three independent implementations, used as a triangle:

- **A. Rebuilt `ereac` with the export patch** (Phase 0). Primary reference for periodic physics.
- **B. JAX port.** Independent implementation; matches V2 to 3.7×10⁻⁷ eV/Å on sampled frames according to the project's own logs. Central-cell only.
- **C. Saved V2 outputs** (`adapter_call_*`, `finite_difference_oracle_force/`). Regression data whose physical setup is settled by task 0.3.

Larger off-equilibrium set (target 200 frames, both periodic and non-periodic, 60C/60El and 112C/113El): frames every 50 fs from classical MD at 300 K and 1000 K; bead configurations from Phase 2 PIMD runs, since these are the geometries production will actually visit; random electron displacements of 0.05 to 0.4 Å from the host nucleus; electrons placed midway between two carbons; compressed and stretched cells. Each frame is evaluated by A, and by B where applicable.

### 5.2 Ladder

| Rung | Test | Reference | Pass criterion |
|---|---|---|---|
| L-F1 | Per-term finite differences (central, h = 1e-4 and 1e-5 Å) on 20 frames, every new term in isolation and in total | The engine's own energy | \|F_FD − F\| ≤ 1e-5 kcal/mol/Å† away from documented discontinuities; any larger value is explained by a named threshold |
| L-F2 | Per-term parity, eReaxFF-specific terms | A, and B in non-periodic mode | Energy ≤ 1e-6 kcal/mol/atom†, force max ≤ 1e-5 kcal/mol/Å† |
| L-F3 | Total energy and force parity on the 200-frame set | A | Threshold frozen after task 0.4; proposal: equal to the plain-ReaxFF baseline difference measured there, with no additional error from electron terms |
| L-E1 | Classical NVE, 20 ps, dt = 0.25 fs | — | Drift ≤ 1e-4 kcal/mol/atom/ps†, no secular heating of either species |
| L-P0 | i-PI P = 1 NVE versus LAMMPS `fix nve`, same start | LAMMPS | Positions agree to 1e-8 Å† after 100 steps (checks unit conversion in `fix ipi`) |
| L-P1 | i-PI patch at P = 1 | Unpatched i-PI | Bitwise identical output |
| L-P2 | All atoms classical at P = 8, NVE | P = 1 run | Positions agree to 1e-10 relative† over 1000 steps |
| L-P3 | Coupled harmonic pair: classical X (spring K to origin), quantum x (spring k to X), P = 4, 8, 16 | Analytic: ⟨X²⟩ = k_BT/K exactly; ⟨(x_j − X)²⟩ = finite-P harmonic variance (the formula already in `harmonic_quantum_variance.f90`) | Both within 3 standard errors and within 2 %† |
| L-P4 | Mixed-system conserved quantity with thermostat on, harmonic and LJ drivers, 10⁵ steps | — | \|Δecons\| ≤ 1e-5 of the kinetic energy scale per 10⁴ steps†, no drift trend |
| L-P5 | Reported temperature of each species in the mixed system | Target T | Within 3 standard errors |
| L-E2 | eReaxFF PIMD, NVE segment after equilibration, P = 16 and 32, dt = 0.1 fs | — | Conserved quantity drift ≤ 1e-3 kcal/mol/atom/ps†; carbon temperature stays within 5 %† over 5 ps |
| L-S1 | Free-particle limit: forces off, P = 16 | R_g² = (ħ²/4mk_BT)(1 − 1/P²) | Within 2 %† |
| L-S2 | PILE-L versus PILE-G | Each other | Static observables agree within combined error |
| L-S3 | P convergence of bead spread, R_g, K_cv, El–C g(r) at P = 8 to 128 | Largest P | Change from P to 2P below 2 %† for the reported P |
| L-S4 | Centroid-virial versus primitive kinetic estimator (quantum atoms only) | Each other | Agree within combined error |
| L-R1 | Split run equals continuous run; two runs with the same seed | Each other | Bitwise with fixed charges; otherwise 1e-10 relative† |

None of the † thresholds is taken from a published standard. They are sized from the precision each stage can actually reach and should be frozen in a tolerance file after first measurement, as Metal-Reax already does.

---

## 6. Critical path and milestones

```
Phase 0 (truth) ──────────────► Phase 3 (engine) ──► Phase 4 (classical) ──┐
                                                                            ├─► Phase 5 (PIMD production)
Phase 1 (i-PI core) ──► Phase 2 (JAX interim, F2 probe) ───────────────────┘
```

- **Start in parallel on day 1:** Phase 0 and Phase 1. They share nothing.
- **Longest chain:** 0 → 3 → 4 → 5, about 34 to 62 days. Task 3.4 is the single largest item.
- **Earliest scientific signal:** task 2.2, after roughly 8 to 14 days. It answers whether the instabilities in the audit (F2) were integrator defects or belong to the force field.

| Milestone | Content | Gate |
|---|---|---|
| M1 | References settled | Task 0.3 table and owner answers to §8 items 1 to 3 |
| M2 | Mixed path integral correct | L-P1 to L-P5 pass in i-PI; in-house engine agrees (task 1.4) |
| M3 | First trustworthy eReaxFF PIMD numbers (slow engine) | L-E2 with the JAX driver |
| M4 | eReaxFF engine validated | L-F1 to L-F3, L-E1 |
| M5 | Classical benchmark reproduced | Phase 4 complete |
| M6 | Production protocol frozen | L-S1 to L-S4, L-R1; chosen P and thermostat documented |

---

## 7. Risk register

| # | Risk | Early warning | Mitigation | Fallback |
|---|---|---|---|---|
| R1 | No configuration of rebuilt `ereac` reproduces the saved V2 data (V-12, V-13) | Task 0.3 table has no matching column | Treat A as truth for physics and C as regression data only | Validate the new engine against A and B; state in any paper that V2 results are not reproducible |
| R2 | The Jamil-target ffield is not the published one | Owner cannot source the ffield; `vpar(38)` = 0 remains unexplained | Obtain the ffield from the paper's SI or the authors | Run with a documented electron taper and report it as an assumption |
| R3 | LAMMPS-flavour ReaxFF differs from standalone beyond tolerance for carbon | Task 0.4 shows differences without a traceable cause | Add a `standalone` compatibility flag set (mass threshold 31, cut-offs, polarisation) | Decision 3 runner-up: JAX kernel as production engine |
| R4 | Task 3.4 forces fail finite differences | L-F1 errors that do not shrink with h | Build term by term with FD at each step; cross-check channel by channel against B | Ship 3.1 to 3.3 first and compare partial sums against A with `pxel2` set to zero |
| R5 | Discontinuities in the force field break conservation under PIMD: `int(vund/2)` in the lone pair (`poten.f:426`), the 1e-12 sharing cut (`reac.f:8872`), bond-order cut-offs, the hard skip at `r > swb` | L-E1 passes classically but L-E2 fails, with jumps coinciding with threshold crossings | Census of threshold crossings as Metal-Reax does; shorten dt | Report NVT results only, with the drift quantified |
| R6 | eReaxFF PIMD unstable even with a correct integrator | Task 2.2 shows heating at dt = 0.1 fs | Cayley propagator, smaller dt, stronger internal-mode friction | Restrict to smaller P or to the compact system while investigating |
| R7 | i-PI overhead becomes comparable to the force time once the engine is fast | Profile shows more than 30 % of wall time in i-PI | shm transport; fewer property outputs | Phase 6 batch driver |
| R8 | i-PI patch not accepted upstream | Review feedback | Keep a pinned fork; the patch is small enough to rebase | All-ring-polymer runs need no patch |
| R9 | Unit conversion mismatch between `fix ipi` and i-PI constants | L-P0 fails at the 1e-6 level | Measure once, document the factor | Rescale in a wrapper only if the mismatch exceeds the L-E tolerances |
| R10 | Stale-variable and list quirks in legacy force code (the real variable `sgn` is not reset per atom, and the element lists in `lonpar` and `elecpos` differ, e.g. `Ag`) | Parity fails only for elements outside C/H/O | Reproduce by default and flag, following Metal-Reax ADR-021 | Restrict supported elements to those tested and refuse others |
| R11 | Charge-solver noise once ACKS2 is enabled | Conserved-quantity drift scales with solver tolerance | Per-bead history via `matching='lock'`; residual-versus-force-error test | Stay with fixed charges |
| R12 | Electrons delocalise over several hosts, making the neglect of exchange invalid | Ring spread approaches the C–C distance | Report and stop interpreting those runs | Out of scope for this project |
| R13 | Legacy-derived code cannot be redistributed | Before any public release | Ask the van Duin group (§8 item 7) | Keep the eReaxFF branch private; publish results and inputs only |

**Audit findings reviewed.** F1 confirmed (V-8, V-9). Corrections and additions to the audit:

1. The addendum infers that `elecpos` uses central-cell distances because it calls `dista2`. `dista2` applies minimum image (V-10). The central-cell behaviour is a property of the V2 reference data, not of the base routine.
2. The audit's description of electron coupling into bonded terms omits the `pxel1`/`pxel2` exponential factor in over-coordination (V-11). The JAX port has it.
3. The El–El formula is stated with geometric means; the base code uses one atom's parameters (V-17). Identical for a single electron type.
4. i-PI has closer per-atom precedents than `fixatoms`: `open_paths` and `bosons` (V-5).
5. "Socket overhead does not matter" holds only for the current slow engine (R7).
6. Not in the audit: V-12, V-13, V-14, V-19, V-20.

---

## 8. Open questions for the owner

| # | Question | Recommended default |
|---|---|---|
| 1 | Which ffield is authoritative for the Jamil benchmark, and what is its electron taper radius? The file in the repo has 0.0 in that slot (V-13). | Obtain the published ffield. Until then use the repo ffield with an explicit 10 Å electron taper, matching what V2 did. |
| 2 | Is the target physics a periodic graphene sheet (image sums) or the central-cell behaviour of the V2 data (V-12)? | Periodic. Reproduce the central-cell numbers once, as a regression check only. |
| 3 | Electron mass in the path integral. | 1 amu, static observables only, with the mass scan of §3.5 reported. |
| 4 | Which quantitative result of the Jamil paper defines "reproduced"? | Host-occupancy statistics and temperature trace at 300 K; owner to name the figure. |
| 5 | Must carbons be classical in the primary result? | Yes (mixed patch), with all-ring-polymer reported alongside. |
| 6 | Charges for production. | Fixed (`icharg 8` equivalent). ACKS2 deferred to Phase 7. |
| 7 | May code transcribed from the legacy Fortran (JAX port, the eReaxFF branch of Metal-Reax) be made public? | Keep private until the van Duin group agrees. |
| 8 | Hardware beyond the M5 Max. | Plan for the M5 Max and a Linux CPU cluster. No GPU cluster work until Phase 6 is justified by measurement. |
| 9 | Are hole particles (`Ho`) needed? | No; refuse them with an explicit error. |

---

## 9. Remaining questions from the brief

### 9.1 Performance and layout (question 9)

Measured today [data]: about 2.1 s per P = 16 batch with JAX and four batches per step (8.5 s/step); about 20 s per batch with the four-worker Fortran pool.

Estimate for the target stack [derived, to be benchmarked first]: a compiled ReaxFF-class evaluation of 225 atoms with fixed charges should cost 3 to 15 ms on one core. With one batch per step and 16 clients for P = 32, that is two rounds, so 10 to 30 ms per step plus a few ms of i-PI overhead. One nanosecond at dt = 0.25 fs (4×10⁶ steps) then takes roughly one to two days, against about a year at the current rate. The gain factors are 4 from one batch per step and of order 100 from the engine.

| Platform | Layout |
|---|---|
| M5 Max (18 cores) | i-PI on one core; up to 16 single-threaded LAMMPS clients over Unix sockets. P = 32 or 64 runs in 2 or 4 rounds. No GPU at this system size. |
| Linux CPU cluster | Same, one node: clients = min(P, cores − 1). Several independent trajectories (different seeds, masses or P) are the efficient use of more nodes; a single trajectory does not scale past P clients. `inet` sockets or i-PI's MPI interface only if one trajectory must span nodes. |
| Future GPU cluster | Requires a non-Metal backend (Metal is Apple-only). The batch API of §3.3 is backend-neutral, but a CUDA or Kokkos kernel is a separate project. |

Benchmark first, in this order: single-geometry CPU-64 time for 225 atoms after Phase 3; i-PI overhead per step with a dummy driver at P = 32; scaling with the number of clients on the M5 Max, since performance and efficiency cores differ.

**Must have for a correct result:** Phases 0, 1, 3 (or 2 as engine), 4, 5. **Speed only:** Phase 6, shm transport, GPU.

### 9.2 Licensing and repositories (question 10)

- Metal-Reax and LAMMPS are GPL-2.0; the plugin is a derivative of LAMMPS and stays GPL-2.0-only. The eReaxFF additions live on a branch of Metal-Reax under the same licence.
- i-PI and its clients are separate processes that exchange data over a socket, so no licence obligation crosses that boundary. The i-PI patch is offered under i-PI's own terms (GPL or MIT, V-23) and is upstreamable.
- The legacy Fortran source carries no licence. Do not redistribute it, the export patch applied to it, or the binaries. Reference outputs generated from it are data and can be kept in a private repository; publishing them is covered by §8 item 7.
- Layout: three places, no monorepo. (1) Metal-Reax, branch `ereaxff`. (2) i-PI fork, branch `classical-atoms`. (3) A small private project repository for inputs, validation scripts, reference data and this plan. The components have different upstreams and licences, and a monorepo would only add vendoring.

### 9.3 Existing Fortran engine and JAX backend (question 11)

- **Fortran engine:** apply the F1 fix and the production-path harmonic test (task 1.4), use it once as the independent check of the mixed scheme, then archive it at a tag. Salvage: `harmonic_quantum_variance.f90` (analytic targets), the normal-mode and free-ring tests as specifications for i-PI cross-checks, the geometry generators and host-analysis scripts in `projects/jamil_replication/`, and the lessons file. Do not carry over the oracle adapter, worker pool or file protocol. Remove `.venv-fast` and build outputs from the 6.3 GB tree before archiving.
- **JAX backend:** keep as the interim engine (Phase 2) and as reference implementation B. It stays central-cell and fixed-charge; do not extend it unless the decision-3 runner-up is triggered.

---

## Appendix A — Change list for Metal-Reax

Design rule kept from the repository: element knowledge is data (ADR-003), reference quirks are reproduced and flagged (ADR-021), unsupported inputs are refused with an error.

**Parser and tables**

| Item | Change |
|---|---|
| `FfieldOptions` | add `bool ereaxff` |
| `SingleBody` | add `pxel2` (atom line 2, col 5), `vincr` (line 3, col 3), `alphag` (line 3, col 8), `betael` (line 4, col 3); these columns are currently skipped as unused |
| `TwoBody` | add `pxel1` (bond line 2, col 8) |
| `GlobalParams` | expose `gp[37]` as `electron_taper` in eReaxFF mode; in that mode the `gp[37]==2` triple-bond predicate must not fire (V-20) |
| `SpeciesFlags` | add `is_electron` (symbol `El`), `is_hole` (`Ho`, refused), `core_charge` (table from `reac.f:4590-4615`), `share_sign` (−1 for C, O, H, F, N, Os; +1 for Li, Lv, P, B, Av, Al; from `reac.f:8815-8851`) |
| Strict missing-pair rule (Q-12) | exempt pairs involving `El`; their bond order is zero because the El radii are negative, as in both codes |
| Heavy-atom threshold | compatibility option for 31 instead of 21 (V-19), only if task 0.4 shows it matters for a supported element |

**Terms (`include/reaxmetal/terms.hpp`)**

- `electron_core(betael, Z, alphag, r, Tap, dTap)` and `electron_electron(betael_El, alphag_El, r, Tap, dTap)` per V-17, with a second taper instance built from (`gp[11]`, `electron_taper`).
- `eshare(vincr, r²) = exp(−vincr·r²)`, summed per core atom into `eshare2`, with the legacy inclusion cut `eshare > 1e-12`.
- `atom_quantities`: take per-atom `valency_c = valency + share_sign·eshare2` and `valency_e_c = valency_e + eshare2` for the lone-pair and over/under quantities only. Bond-order correction and the valence-angle lone-pair adjustment keep the uncorrected values (V-18).
- `over_under`: replace `Δ_i` by `exp(−pxel2_i·eshare2_i·S_i)·(abo_i − valency_c,i)` with `S_i = Σ_bonds pxel1·BO/abo_i`, and the same for each neighbour inside the π-bond sum (`poten.f:920-949`). Derivatives with respect to bond orders, to `abo`, and to the electron positions through all three channels.

**Drivers**

- `src/cpu/nonbonded.cpp`: electron pair loop over the existing far list, pairs with exactly one or two `El` members, cutoff `electron_taper`. Shielded Coulomb and van der Waals for El pairs are unchanged.
- `src/cpu/bonded.cpp`: a pre-pass before the atom-quantities block builds `eshare2` and a short core–El pair list for owned and ghost cores; lone-pair and over/under blocks use the corrected quantities; force contributions on electrons are accumulated like any other atom.
- Ghost shell: over-coordination of an owned atom needs `eshare2` and the full bond list of its bonded neighbours. Required shell = existing `2·bond_cut` requirement plus the sharing range √(ln 10¹²/vincr). Checked in `init_style`, error if violated.
- Periodicity: sharing sums over all images in the ghost set. The legacy code takes the nearest image only (V-10). With `vincr` = 181 Å⁻² the two differ by less than 10⁻³⁰⁰ for any cell larger than 1 Å, so this is not a compatibility issue for the ffields on disk; a small-`vincr` ffield would need a note.

**Adapter and outputs**

- `pair_style reaxff/metal NULL ereaxff yes`. The style name is unchanged so that charge fixes still match `^reaxff`.
- `EnergyBreakdown` gains `Electron`; `nextra` becomes 15 in eReaxFF mode.
- Fixed charges: `checkqeq no` with charges set in the data file is the `MOLCHARGE`/`icharg 8` equivalent. The polarisation energy that LAMMPS always adds must be switched off in this mode, because the legacy total includes the charge energy only when `ncha2 = 1` and the saved V2 terms show `charge = 0` (V-13 data file).
- `backend metal` with `ereaxff yes`, and any `Ho` atom, raise `UnsupportedFeatureError`.
- New rows in `capabilities.cpp` and `FEATURE_MATRIX.md`: `term.electron`, `term.electron_sharing`, `ffield.ereaxff`, `compat.ereaxff_share_sign`, `compat.core_charge_table`.

---

## Appendix B — Change list for i-PI

| File | Change |
|---|---|
| `ipi/engine/motion/motion.py`, `dynamics.py` (`Dynamics.__init__`) | new optional argument `classical_atoms` (atom indices), stored as `classical_dof` (3 entries per atom) |
| `ipi/inputs/motion/motion.py`, `ipi/inputs/motion/dynamics.py` | input field `<classical_atoms>`, store and fetch so that it survives restart |
| `dynamics.py`, `DummyIntegrator.bind` (≈ `:276-306`) | keep `classical_dof`; set `nm.qnm[1:, classical_dof] = 0` and `nm.pnm[1:, classical_dof] = 0`; raise for NPT/NST/SC/constrained-centroid ensembles, GLE thermostats, the `bab` propagator, and overlap with `bosons`, `open_paths` or `fixatoms_dof` |
| `dynamics.py`, `DummyIntegrator.pconstraints` (`:368-401`) | after the existing blocks: add `0.5·Σ pnm[1:, cd]² / dynm3[1:, cd]` to `ensemble.eens`, then set `nm.pnm[1:, cd] = 0` |
| `ipi/engine/properties.py`, `get_temp` (`:1116-1175`) | count (P − 1) frozen modes per classical degree of freedom as fixed, following the `fixatoms_dof` branch |
| `ipi/engine/properties.py`, `get_kintd` | exclude classical atoms from the primitive estimator's constant term, or refuse the property when classical atoms are selected; the centroid-virial estimator is already correct, since `q_j − q̄ = 0` gives 3k_BT/2 per classical atom |
| `ipi_tests/regression_tests/tests/` | new cases for L-P1 to L-P5 |
| docs | one section stating the scheme, its equivalence proof (Appendix C) and its restrictions |

No change is needed in `normalmodes.py`, `thermostats.py`, `forces.py` or `ensembles.py`. The centroid degrees of freedom of classical atoms remain thermostatted, so thermostat degree-of-freedom counts are unchanged; noise applied to their internal modes is removed by the projection and booked in `eens`, exactly as i-PI already does for fixed atoms.

---

## Appendix C — Equations and algorithm

**Convention (i-PI).** N atoms, P beads, physical temperature T, T_P = P·T, ω_P = P·k_B·T/ħ. Let 𝒬 be the quantum atoms and 𝒞 the classical ones.

```
H_P = Σ_i Σ_j  p_ij² / (2 m_i)
    + Σ_{i∈𝒬} Σ_j  ½ m_i ω_P² |q_ij − q_i,j+1|²
    + Σ_j  V(q_1j, …, q_Nj)

constraint for i ∈ 𝒞:   q_ij = Q_i  and  p_ij = P_i   for all j
sampled density:         exp( −H_P / (k_B T_P) )
```

**Equivalence for classical atoms.** For i ∈ 𝒞 the kinetic term is P·P_i²/(2m_i), so the density in (Q_i, P_i) is exp(−β[P_i²/(2m_i) + (1/P)Σ_j V_j]) with β = 1/(k_BT): physical mass, temperature T, and the bead-averaged potential. In normal modes, q̃_0 = √P·Q, p̃_0 = √P·P, f̃_0 = Σ_j f_j/√P, and the equation dp̃_0/dt = f̃_0 reduces to dP/dt = (1/P)Σ_j f_j. Quantum beads feel the unscaled force f_ij. This is the fix for F1: no 1/P on quantum beads, mean force on classical atoms.

**Normal modes.** Orthonormal closed-ring transform C. ω_k = 2ω_P·sin(πk/P), k = 0 … P−1. Projection Π: set p̃_k = 0 (and q̃_k = 0) for k ≥ 1 and i ∈ 𝒞.

**One step of length dt (BAOAB; forces f from the previous step are reused, so there is one force batch per step).**

```
B:  p ← p + (dt/2) f                         ; Π
A:  for each k: exact rotation of (p̃_k, q̃_k) by ω_k dt/2
        k = 0:  q̃_0 ← q̃_0 + (dt/2) p̃_0 / m
        k ≥ 1:  p̃_k ←  cos(ω_k dt/2) p̃_k − m ω_k sin(ω_k dt/2) q̃_k
                q̃_k ←  sin(ω_k dt/2) p̃_k /(m ω_k) + cos(ω_k dt/2) q̃_k     (old values on the right)
O:  p̃_k ← c1_k p̃_k + sqrt(m k_B T_P (1 − c1_k²)) ξ ,   c1_k = exp(−γ_k dt)   ; Π
        γ_k = 2 λ ω_k for k ≥ 1 (λ = 1 is critical damping and the i-PI default `pile_lambda`; configurable)
        k = 0: PILE-L  γ_0 = 1/τ_0 ;  PILE-G  stochastic velocity rescaling of the total centroid kinetic energy
A:  as above
F:  f_j = −∇V(q_·j) for all j in one batch
B:  p ← p + (dt/2) f                         ; Π
```

The Cayley variant replaces the rotation in A by its Cayley transform and removes the internal-mode resonance at ω_k·dt = π; it is available in i-PI (`normalmodes.py:540-547`).

**Conserved quantity.** H_P + (heat removed by thermostats) + (kinetic energy removed by Π and by centre-of-mass removal). In i-PI this is `econs` with the projection booked in `eens`.

**Estimators.**

```
⟨V⟩      = (1/P) Σ_j V_j
K_cv     = (3N/2) k_B T − (1/2P) Σ_j Σ_i (q_ij − q̄_i) · f_ij          (classical atoms contribute 3 k_B T / 2)
R_g,i²   = (1/P) Σ_j |q_ij − q̄_i|²          bead spread = sqrt(R_g²)
free ring: R_g² = (ħ² / 4 m k_B T)(1 − 1/P²)     (0.0404 Å² at 300 K, 1 amu, large P)
```

**Timestep guidance [derived].** With exact free-ring propagation the limit comes from the physical forces: ω ≈ 0.7 fs⁻¹ for a 1 amu electron in its well gives ω·dt = 0.18 at 0.25 fs. Internal-mode resonances need ω_max·dt = 2ω_P·dt < π; at 300 K and dt = 0.25 fs this is 0.63 for P = 32 and 1.26 for P = 64. The 0.01 fs stability window in the audit is therefore not expected from the physics and is consistent with the integrator defects F1 and F3.
