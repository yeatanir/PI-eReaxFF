# PI-eReaxFF — Implementation Audit Report

Audited tree: `/Users/yeatanir/PhD/eReaxFF PIMD/ereaxff_pimd_modular` (git HEAD `7ff9d90`)
Audit date: 2026-11-10
Audience: an external expert reviewer of the path-integral (PIMD/RPMD) engine and its eReaxFF coupling.

---

## 0. Scope, method, and limits of this audit

**What was read:** every Fortran module under `src/` (core, pimd, thermostat, integrator, force, output, drivers), the two Fortran tools that bind to the legacy code (`tools/persistent_oracle_worker.f90`, `tools/ereaxff_energy_terms_export.f`, `tools/oracle_single_point_force_export.f90`), the Makefile, the YAML configs, the project docs (`README`, `DESIGN`, `PROJECT_STATUS`, `VALIDATION`, `LESSONS_LEARNT`, `docs/*`), the Jamil-replication run notes and summaries, the JAX fast backend (`projects/jamil_replication/fast_backend/*.py`), and the earlier i-PI comparison report.

**What was NOT possible, and why it matters:**

1. **The legacy eReaxFF V1/V2 Fortran source is not in this directory and not on this machine.** `Makefile` and every config point to `/Users/yeatanir/PhD/PIMD/ereaxff_pimd/{V1,V2}`; that path does not exist here. The i-PI checkout the earlier comparison cites (`/Users/yeatanir/PhD/PIMD/i-pi`) is also absent. Consequently:
   - Statements about the **legacy side** (ACKS2 solver, how `charge_in`/`solver_guess` are consumed, exact energy-term formulas) come from (a) the thin Fortran wrappers that call it, (b) the `ffield`/`control` files in `build/`, and (c) the JAX re-implementation, whose docstrings say its formulas are transcribed from the legacy `poten.f` / `reac.f`. They are flagged **[indirect]** below.
   - The code cannot be rebuilt on this machine as it stands (`make persistent_oracle_worker` needs `$(V2_ROOT)/build/*.o`). Pre-built binaries exist in `build/`.
2. **I did not run the test suite or any simulation.** All "PASS" statements are quoted from the project's own docs, logs, and summaries.
3. **I did not compare against i-PI source.** Where I compare with i-PI it is from the project's own earlier report plus my derivation.

Confidence labels used: **[code]** = read directly in source; **[doc]** = claimed by project documentation/logs; **[indirect]** = inferred via wrapper/JAX port; **[derived]** = my own derivation, should be checked.

---

## 1. Executive summary

**What exists.** A modular Fortran 2018 engine ("PI-eReaxFF") that propagates selected explicit eReaxFF pseudo-electrons (`El`) as closed ring polymers (P beads, standard-i-PI ω_P convention, exact free-ring-polymer propagation, BAOAB splitting, PILE-L Langevin thermostat in normal-mode space), while the remaining atoms (graphene carbons) move classically with Berendsen or Langevin thermostats. Forces/energies/charges come from the unmodified legacy eReaxFF V2 force oracle, called once per bead per force stage through a file-based protocol to persistent worker processes. An experimental JAX re-implementation of the force field ("fast backend") exists for speed.

**What is solid.** Normal-mode transform, ω_k, free-ring propagation, OU/PILE-L thermostat statistics, heat ledger, harmonic finite-P variance benchmark (≈2% agreement at P=16), group selectors, restart round-trip, and the plumbing to the legacy oracle (P=1 reproduction, repeatability, NaN rejection, parallel dispatch) are implemented and tested.

**What is not established, in order of importance for a reviewer:**

| # | Severity | Finding |
|---|---|---|
| F1 | **High** | **Force scaling, temperature, and spring conventions appear mutually inconsistent in the production path** (§6.4). Springs use ω_P = P·kT/ħ and the thermostat runs at P·T (the "P·T convention", which requires *unscaled* bead potential V), but `bead_force_manager` scales bead forces by 1/P. The harmonic validation bypasses the production path (`physical_force_scale=1.0`). The combination actually used for eReaxFF was never validated. [code + derived] |
| F2 | **High** | **No demonstration of energy conservation or correct sampling with eReaxFF forces on.** Logged runs show carbon heating to 1156–2139 K at dt=0.05 fs (stable only at 0.01 fs), "conserved-like" spans of 300 → 10 000 eV, and one run labelled PASS with a 303 eV span (§12). [doc] |
| F3 | **High** | **Integrator is not a clean BAOAB for the coupled system**: 4 force batches per step (≈3× redundant), and electrons' first half-kick uses forces at a mixed configuration (new carbons, old electron). First-order, non-symplectic (§3.3). [code] |
| F4 | **Medium** | PILE-G (global SVR centroid thermostat) is **not implemented** (hard `error stop`); no SVR; carbons use Berendsen (non-canonical). The two baths never exchange anything directly (§7). [code] |
| F5 | **Medium** | The "conserved-like" diagnostic uses only the *first* electron's KE/spring, excludes all classical KE and all other electrons, and adds extended-temperature KE to V/P (§9.1). [code] |
| F6 | **Medium** | Charge equilibration is a black box in this repo: Jamil production uses `icharg=8` with fixed `MOLCHARGE`-style charges; bead-local vs centroid ACKS2 policies are plumbing that is physically close to trivial for that target (§8). [indirect] |
| F7 | **Medium** | Reproducibility: legacy source and i-PI absent; absolute hard-coded paths; 6.3 GB tree with 36 591 files including `.venv-fast` and build outputs (§11). |
| F8 | **Low–Med** | Unit constants duplicated at different precision (kcal→eV 0.043364 vs 0.0433641153; force factor 9.6485e-3) → ~3–5×10⁻⁶ relative energy/force inconsistency (§4.3). |
| F9 | **Low–Med** | Fictitious electron mass (1.0 amu in Jamil configs; physical mass in others), initial beads collapsed to a point, single Park–Miller LCG RNG, no quantum estimators beyond RMS bead spread, no exchange (§5, §10). |

---

## 2. Repository map (what lives where)

```
ereaxff_pimd_modular/
  src/core/        kinds, constants, units, errors, rng, system_state, selectors,
                   particle_groups, legacy_geo_loader, input_config (YAML subset), restart
  src/pimd/        normal_modes, ring_polymer, spring_energy, free_ring_propagator,
                   harmonic_quantum_variance
  src/thermostat/  thermostat_base, _none, _langevin, _pile_l, _pile_g (stub), _registry
  src/integrator/  baoab_integrator (harmonic test only), mixed_baoab_integrator (production quantum step)
  src/force/       force_provider_base, dummy/harmonic providers, bead_force_manager,
                   charge_coupling_policy, ereaxff_oracle_adapter (legacy bridge), fast_jax_adapter
  src/output/      pimd_output, mixed_restart, mixed_baoab_restart, reaxff_compat_output (fort.7/8, molfra, xmolout)
  src/drivers/     ereaxff_mixed_baoab_main (production driver, 1239 lines), harmonic_pimd_main, ereaxff_pimd_main
  src/ipi_socket/  8-line stub (not implemented)
  tools/           persistent_oracle_worker.f90, oracle_single_point_force_export.f90, ereaxff_energy_terms_export.f
  projects/jamil_replication/   Jamil et al. graphene benchmark configs, geometries, run notes, JAX fast backend
  tests/           test_00 … test_35
```

Production entry point: `build/ereaxff_mixed_baoab_main <config.yaml> [out_dir]`.

---

## 3. The PIMD loop

### 3.1 Initialisation (driver, `ereaxff_mixed_baoab_main.f90`)

1. Parse YAML subset → `input_config` (`core/input_config.f90`).
2. Load geometry (`user_geo_structure`: legacy `geo` HETATM + `CRYSTX`, `legacy_geo_loader.f90`) or built-in 60-atom phase-4 graphene fixture; apply cell from `CRYSTX`.
3. Assign particles to groups with selectors (`element == C`, `element == El`, `id == N`, `id in [...]`, one `and`); each group has `dynamics ∈ {classical, pimd, frozen}`, per-group thermostat, temperature, beads, mass. A particle matching zero or multiple groups is a hard error.
4. One `ring_polymer_state` per quantum particle: `cmat` (P×P orthonormal DFT matrix) and `omega_k` (§6).
5. Bead positions: every bead of every particle starts at the geometry position (`build_positions`, line 469) — **the ring starts collapsed to a point**. (The non-config legacy test path applies a deterministic linear spread, `perturb_pimd_beads`.) `q_nm = C·q_bead`.
6. Velocities (`initialize_group_velocities`, line 607): classical atoms ~ N(0, √(kT/m)); each quantum normal-mode velocity ~ N(0, √(kT·**P**/m)) (extended temperature P·T). Optional one-shot rescale of centroid and internal-mode KE to exactly 3/2·kB·(P·T) per mode block (`rescale_initial_pimd_velocities`).
7. Create the oracle adapter (spawns persistent worker pool), charge memory (`charges_by_bead(N,P)`, zero-initialised), solver guess (zero).
8. Thermostats: PILE-L state per quantum particle initialised with **temperature = T·P** (line 952); classical groups keep their config parameters.
9. One `initial_force` batch.

### 3.2 One MD step (`mixed_group_baoab_step`, driver line 975; `mixed_baoab_quantum_step`, `integrator/mixed_baoab_integrator.f90`)

```
 (a) replicate classical positions into all bead slots
 (b) F1 = oracle(all beads)                         stage "group_pre_force"
 (c) classical atoms:  v += ½dt·F̄1/m ;  x += ½dt·v        (B, A)      F̄ = Σ_j F_j/P
 (d) classical thermostat step (Berendsen / Langevin / none)           (O)
 (e) classical atoms:  x += ½dt·v ; replicate to all beads             (A)
 (f) quantum step, for every PIMD particle q:
       update bead positions from q_nm (wrap each bead into the cell)
       F2 = oracle(all beads)                       stage "quantum_pre_force"
       v_nm += ½dt · C·(F2_j/P · P_max/P_q) / m                        (B)
       exact free-ring propagation  ½dt   (A)
       PILE-L step  (full dt coefficients)                             (O)
       exact free-ring propagation  ½dt   (A)
       update bead positions from q_nm
       F3 = oracle(all beads)                       stage "quantum_post_force"
       v_nm += ½dt · C·(F3_j/P …) / m                                  (B)
 (g) replicate classical; F4 = oracle(all beads)    stage "group_post_force"
 (h) classical atoms:  v += ½dt·F̄4/m                                   (B)
```

Wrapper nesting: the ring-polymer is thermostatted in normal-mode space; the free-ring part (spring + centroid drift) is integrated *exactly* (rotation by ω_k·dt/2 in each internal mode, free drift for the centroid) — `free_ring_propagator.f90`. Forces are applied as kicks in normal-mode space after the bead→NM transform (`half_kick_from_summary`).

Force evaluations per step: **4** batches (b, f-pre, f-post, g) plus one `initial_force` at start → "80 oracle calls for P=16, one step" in the run notes (5×16). [code, doc]

### 3.3 Integrator-correctness observations

- **Redundant evaluations.** F4 is evaluated at exactly the geometry of F3 (nothing moves between them), and F1 at step n+1 is the same geometry as F4 at step n. A velocity-Verlet/BAOAB step needs 1 batch, not 4. This is a 3–4× avoidable cost on the dominant bottleneck (the oracle).
- **Mixed-configuration kick.** F2 is evaluated after the carbons have already been moved (steps c–e) but before the electrons move. The first electron half-kick therefore uses F(C_new, e_old) rather than F(C_old, e_old) (or a consistent re-ordering). The classical atoms' kicks (F1 at old/old; F4 at new/new) are consistent. The effective splitting is not symmetric/reversible for the coupled system; the electron kick error is O(dt) per half-kick → first-order global accuracy, no symplecticity. This is a plausible contributor to the heating seen in §12, though it has not been isolated.
- **Per-quantum pre/post-force batches** are re-run inside `mixed_baoab_quantum_step` once for *all* quantum particles together (it loops `q` after the shared force evaluation) — OK, but with multiple electrons each "bead j" replica contains all electrons' bead j (standard RPMD pairing).
- Project notes record that an extra classical half-kick was already removed once (2026-09-11 audit); the code I read does not contain it.

---

## 4. How PIMD connects to eReaxFF

### 4.1 Architecture

```
 ereaxff_mixed_baoab_main  ──►  bead_force_manager.evaluate_oracle_beads
                                   │   (dispatch: serial | openmp_all | bounded_pool)
                                   ▼
                           ereaxff_oracle_adapter          (src/force)
                                   │   file protocol: request.in / response.ready
                                   ▼
                         persistent_oracle_worker (one process per pool slot)   (tools/)
                                   │   Fortran CALLs
                                   ▼
                  legacy V2: ereaxff_force_oracle_init / ereaxff_force_eval
                              (+ reac_v2, poten_v2, blas, cg, shanno_v2, minres, link stubs)
```

Design rule (stated in `DESIGN.md`, respected by the code [code]): the oracle evaluates energy, forces, charges only; it owns no dynamics, thermostat, or trajectory state.

The **oracle sees a full classical system per bead.** For bead j, the `positions_A(:, :, j)` slice holds every particle: carbons (identical in all bead slots, replicated), and for each electron its j-th bead position. So the oracle is called P times per stage, each call a complete N-atom (e.g. 225 = 112 C + 113 El) eReaxFF evaluation. Electron–electron and electron–carbon interactions are therefore only between *same-index* beads (standard RPMD for distinguishable particles; no exchange).

### 4.2 Data crossing the boundary

| Direction | Quantity | Units | Notes |
|---|---|---|---|
| in | positions (3×N) | Å | binary double via worker request file; for the legacy external-process path they go through a `geo` file written with `3f10.5` (1e-5 Å quantisation) |
| in | charges_in(N) | e | bead-local charge memory |
| in | solver_guess_in(N) | solver units | ACKS2 warm-start vector (see §8) |
| out | energy | kcal/mol | × `kcal2eV = 0.043364` → eV |
| out | forces (3×N) | legacy "force_work" | ÷ `legacy_force_accel_scale = 9.6485e-3` → eV/Å (this constant is (1 eV/Å)/amu in Å/fs², i.e. forces are in amu·Å/fs²) [code] |
| out | charges(N), solver_out(N) | e | committed to bead-local memory after the batch |
| out | `status`, 14 energy terms, bond orders, nonbonded buckets | kcal/mol | extra diagnostics exported by `ereaxff_energy_terms_export.f` |

### 4.3 Unit handling findings

- Energy: `kcal2eV = 0.043364` (adapter, worker, JAX) vs `kcalmol_to_eV = 0.0433641153087705` in `core/constants.f90` — relative difference 2.7×10⁻⁶.
- Force: 9.6485e-3 vs exact 9.648533e-3 (3.4×10⁻⁶ relative). Combining the two, the force and energy conversions are mutually inconsistent at ~5×10⁻⁶ relative. Negligible physically but non-conservative at the level of a precision audit; the fix is a single shared constant.
- Time/mass: `amu_angstrom2_per_fs2_to_eV = 103.6427…` computed from SI constants [code].

### 4.4 Bead-force reduction (`bead_force_manager.f90`)

For bead j: `scaled_forces(:,:,j) = F_j / P` (lines 82, 241, 278, 336, 576) and `pimd_potential = (1/P)·Σ_j E_j`. Classical atoms receive `Σ_j scaled_forces = mean_j F_j` (mean-field force from the electron ring) — correct for classical nuclei coupled to a ring polymer. **Electron beads also receive F_j/P** — see F1 in §6.4.

For particles with fewer active beads than `max_beads`, unused bead slots are centroid-filled and forces are rescaled by `P_max/P_active` (`half_kick_from_summary`, lines 157–165).

---

## 5. How the beads are formed

- **Representation.** Per quantum particle and Cartesian component, the state is held in **normal-mode coordinates** `q_nm(3, nquantum, P)` and `v_nm`. Bead positions are derived on demand: `q_bead = Cᵀ q_nm`.
- **Transform** (`pimd/normal_modes.f90`, `build_cmat`): real orthonormal closed-ring DFT matrix. Row 1 = centroid, `1/√P`; for k = 1..⌊(P−1)/2⌋: row 2k = √(2/P)·cos(2πk n/P), row 2k+1 = √(2/P)·sin(2πk n/P), n = bead−1; for even P the last row is the Nyquist mode `(−1)^n/√P`. Centroid = `q_nm(1)/√P` = bead average. Tests assert `C Cᵀ = I`, round-trip, centroid consistency, up to P=128.
- **Initial ring.** All beads at the same point (collapsed), with thermal normal-mode velocities at P·T. The ring is expected to expand during dynamics/thermostatting; no separate ring pre-equilibration is performed. Initial bead RMS spread is therefore exactly zero. At P=16 a short run only reaches ≈0.005 Å spread after one step (`jamil_p_sweep_summary.out`), 0.05–0.15 Å in 10–100-step runs — i.e. far from equilibrium ring size for most logged runs.
- **Per-particle bead counts** are supported (`quantum_bead_counts`) with a shared storage axis `max_beads`; inactive slots hold the centroid.
- **Periodic wrapping.** Each bead is wrapped into the cell independently *only for the oracle input*; `q_nm` stays unwrapped so spring energies/forces are not corrupted by wrapping. [code]
- **Beads are not copied for classical atoms**: classical atoms exist once (`positions(:,i,1)`) and are replicated into all bead slots before every oracle call.

---

## 6. Spring constant and ring-polymer Hamiltonian

### 6.1 Frequencies

(`normal_modes.f90`, lines 75, 88; units fs⁻¹ with kB in eV/K and ħ = 0.6582119569 eV·fs)

```
ω_P  = P · kB · T / ħ
ω_k  = 2 · ω_P · sin(π k / P),   k = 0 … P−1,    ω_0 = 0
```

`T` is the **physical** group temperature. The legacy prototypes used `√P · kB T/ħ` — the new code deliberately replaced it (guarded by a regression test), see F1 for why this change needs care.

### 6.2 Spring "constant"

There is no stored spring constant. The ring-polymer springs are the diagonal normal-mode quadratic form with stiffness per mode

```
k_k = m · ω_k²          [amu·fs⁻² → ×103.6427 → eV/Å²]
E_spring = Σ_{k≥1} ½ m ω_k² q_nm,k²      (spring_energy.f90, centroid excluded)
```

which equals the bead-space form ½ m ω_P² Σ_j |q_j − q_{j+1}|² with ω_P = P kT/ħ (cyclic). Worked numbers at T = 300 K, P = 16: ω_P = 0.6284 fs⁻¹ (≈3340 cm⁻¹), ω_max = 2ω_P = 1.257 fs⁻¹ (≈6670 cm⁻¹, period ≈ 5 fs), k_max = 163.7 eV/Å² for m = 1 amu, 0.090 eV/Å² for the physical electron mass 0.00054858 amu.

### 6.3 Mass

Each quantum particle has its own mass (config `mass_amu`, required for every PIMD group). Jamil configs use **1.0 amu** for El (matching the paper's pseudo-electron mass, which is a classical-eReaxFF convention, not a physical electron mass); the paired-graphene config uses the physical 0.00054858 amu. Because ω_k depends only on P·T (not mass), the mass affects only the kinetic timescale and bead spread. The project notes this ("force-off mass checks…").

### 6.4 **F1 — Convention consistency (high-priority review item)**  [code + derived]

The standard primitive path-integral discretisation is exp(−Σ_j [m P/(2βħ²)(q_{j+1}−q_j)² + (β/P)V(q_j)]). It can be sampled by classical dynamics in two self-consistent ways:

| Convention | spring ω | potential per bead | thermostat/velocity T | 
|---|---|---|---|
| "P·T" (i-PI, Ceriotti) | `P kT/ħ` | **V(q_j)** (unscaled) | **P·T** |
| "Tuckerman" | `√P kT/ħ` | **V(q_j)/P** | **T** |

What the production path does (all [code]):
- ω_P = P kT/ħ (`normal_modes.f90:75`) → P·T convention;
- thermostat initialised at `T·P` and velocities drawn at `kT·P` (`main.f90:952, 635`) → P·T convention;
- bead forces `F_j/P` (`bead_force_manager.f90:241`, restored to `F_j/P_active` in the kick) → **Tuckerman-style potential**.

With ω_P and T_ext fixed as above, the sampled weight of the potential is exp(−V̄/(kB·P·T)) = exp(−ΣV_j/(P²kT)), i.e. the beads see a potential **P times too weak** (for the electrons). The harmonic validation (`test_03`, `harmonic_pimd_main`) passes `physical_force_scale = 1.0` with T_ext = P·T, which is the self-consistent P·T combination — so that test passes, but it does not test the combination used with eReaxFF. `DESIGN.md` itself warns that these three choices "are coupled" and `LESSONS_LEARNT.md` (2026-04-30 entries) records the harmonic path choice, but the production driver was never reconciled.

Observation that raises the stakes: the *legacy* prototypes (`√P` frequency, `1/P` forces, physical-T velocities) are exactly the self-consistent Tuckerman convention. The "convention bug fix" flipped only the frequency (and added P·T temperature) without flipping the force scaling; if my reading is right, the legacy code was consistent and the new one is not. The classical atoms are fine in either convention (they get the mean force Σ_j F_j/P at temperature T).

**Suggested resolution/test (not performed):** run the *production* `mixed_baoab_quantum_step` with a harmonic force provider through `bead_force_manager` (so the 1/P is in the loop) at P = 8/16 and compare ⟨x²⟩ to `harmonic_finite_p_bead_variance_A2`; either remove the 1/P for electron beads (keeping a mean-field 1/P only for classical atoms) or switch to the Tuckerman parameter set. Please verify the derivation against i-PI source before changing code, since I could not.

> **Status update (Part II, §II.1):** the P·T convention below is now **confirmed in the current i-PI source** (`ipi/engine/normalmodes.py:461`, `dynamics.py:235`, `ensembles.py:259`). F1 stands as a real inconsistency, and §II.3 gives the minimal fix.

---

## 7. Thermostats and how they "communicate"

### 7.1 What is implemented

| Name | Applies to | Where | Notes |
|---|---|---|---|
| `none` | any group | `thermostat_none` | zero heat |
| `langevin` | classical groups | `apply_classical_thermostat` (driver line 1115); `thermostat_langevin.f90` (unit-tested) | exact OU: v ← e^{−γdt}v + √(kT/m)·√(1−e^{−2γdt})·ξ, per Cartesian component |
| `berendsen` | classical groups | `apply_berendsen_group` (line 1190) | weak coupling: λ = √(1 + dt/τ (T₀/T − 1)); group-wide rescale; **non-canonical**, suppresses KE fluctuations |
| `pile_l` | PIMD groups | `thermostat_pile_l.f90`, used inside `mixed_baoab_quantum_step` | local Langevin per normal mode, see below |
| `pile_g` | — | `thermostat_pile_g.f90` | **stub: `error stop`** — deliberately fails rather than silently mapping to PILE-L |
| `svr` | — | listed in `docs/THERMOSTATS.md` | not in the registry; not implemented |

**PILE-L** [code]: mode friction γ_0 = `centroid_gamma_fs_inv` (config; Jamil: 0.001 fs⁻¹, i.e. τ = 1 ps); γ_k = 2λ ω_k for k ≥ 1 with **λ = 1 hard-wired** (not configurable); c1_k = e^{−γ_k dt}; v_nm,k ← c1_k v + √(k_B T_ext/m)·√(1−c1_k²)·ξ with **T_ext = P·T**. Applied once per step with full dt between two half-step free propagations (BAOAB). The same unit-tested OU kernel semantics.

**Heat ledger.** Every thermostat step accumulates `heat += KE_before − KE_after` (energy removed from the system). Classical heat and per-quantum heat are kept separately and summed into `thermostat_heat_eV`; both are written to restart files. The ledger is exact for OU/PILE-L; for Berendsen it is the KE change of the rescale.

### 7.2 Do the two baths communicate?

**Not directly, by design.** In the code there is no coupling term, no shared state, and no energy exchange routine between the classical bath and the ring-polymer bath:

1. They act on disjoint degrees of freedom (carbon Cartesian velocities vs electron normal-mode velocities).
2. They are applied in sequence inside one step: classical O (step d) → electron O inside the quantum step (f). The only communication is **through the oracle forces**, which include the C–El interaction (so the baths exchange energy only physically, through the force field).
3. They run at **different effective temperatures** — carbons at T, electron normal modes at P·T (the extended-ensemble convention) — which is consistent only if the P·T convention is used consistently (F1). The physical electron temperature is reported as `T_ext/P`.
4. They **share one RNG stream** (`rng_state` passed through both). This couples the random sequences (order-dependent reproducibility) but does not couple the physics. The generator is a Park–Miller minimal-standard LCG (16807 mod 2³¹−1, period ≈2.1×10⁹) with Box–Muller using only the cosine branch — adequate for tests, thin for production-length runs of ~10⁸ draws.
5. There is no PILE-G, hence no global/centroid-to-carbon-bath coupling and no momentum/COM removal for the mobile carbon patch.

Since the question of "how the PILE (ring-polymer) thermostat talks to the regular particle thermostat" has the honest answer *"it doesn't, except via forces and a shared RNG"*, a reviewer should confirm that this separate-bath scheme is what is intended for a mixed classical/ring-polymer system with a strongly coupled fast electron subsystem — hot-electron / energy-partition problems are documented in the run notes (carbon heating, §12).

### 7.3 Registry/doc mismatches

`thermostat_registry.f90` accepts `none, langevin, berendsen, pile_l, pile_g`; `docs/THERMOSTATS.md` lists `svr` and omits `berendsen`. Docs lag code.

---

## 8. Charge equilibration (QEq / ACKS2)

### 8.1 What the PIMD engine owns  [code]

`charge_coupling_policy.f90` + `bead_force_manager.f90`:

- `bead_local_acks2` (default production): one independent charge vector and solver guess per bead (`charges_by_bead(N,P)`, `solver_guess_by_bead(N,P)`). Each bead is an independent oracle call; results are committed to memory in deterministic bead order *after* the batch completes (no races; persistent-worker design doc states this explicitly).
- `centroid_acks2`: one oracle call at the centroid geometry (`evaluate_oracle_centroid`), with charge memory and solver guess averaged over beads as input; the single result is broadcast to all beads; `pimd_potential` = that single energy; forces are ÷P and broadcast. This is a documented *approximation* (it also means the electron beads all feel the force of the centroid geometry).
- `frozen_charges`, `no_charge_coupling`: valid names but rejected by the oracle manager.
- Diagnostics: `charge_by_bead.out`, `charge_centroid.out`, max charge difference across beads.
- Initial charges and solver guesses are zero; the first call must converge from a cold start.

### 8.2 What the legacy oracle does  [indirect]

- The ffield header reads `[ ereaxff acks2 ]`; the V2 link line includes `cg.o`, `minres.o` (Krylov solvers) and `shanno_v2.o` — consistent with an ACKS2/EEM matrix solve done in the legacy code. The `ffield` carries EEM parameters (`gammaEEM`, `chiEEM`, `etaEEM`) per element and an `El` atom type.
- The wrapper passes `charge_in(N)` and `solver_guess_in(N)` into `ereaxff_force_eval(...)` and gets `charge_out`, `solver_out`. **I could not verify how the legacy code uses `charge_in` / `solver_guess_in`** (initial guess, constraint, or ignored). Note the guess has length N, while a full ACKS2 solve has 2N+2 unknowns (charges + potentials + 2 Lagrange multipliers), so `solver_guess` is presumably a reduced vector; this needs checking in the legacy source.
- The Jamil control file in the work directories has `icharg = 8`. Per the JAX backend docs, **`icharg=8` with explicit `MOLCHARGE` constraints is the actual charge path for the Jamil 112C/113El target**, and the JAX charge kernel simply returns the fixed labelled charges; "full dynamic ACKS2 matrix construction/solve is fail-closed". `targets.md` says fixed C⁺/e⁻ pairs "can avoid ACKS2". Net charge in logs: ≈ −1.4×10⁻¹³ e (neutral).

**Implication.** For the Jamil target, charge coupling is essentially fixed per particle, so bead-local vs centroid ACKS2 — the elaborate part of the charge machinery — is exercised in plumbing and tests but is physically nearly trivial there. A dynamic-ACKS2 PIMD validation (charges that really respond to bead geometry) is not evidenced for any production-sized system.

### 8.3 Charge-related risks for review

- Cold-start solver guess with nonzero tolerance (`charge_tolerance: 1e-8` is parsed but I did not find where it is applied).
- Under `centroid_acks2`, averaging bead charge memories and solver guesses is not a theory-derived contraction.
- JAX fast backend supports only bead-local policy.

---

## 9. Energy terms

### 9.1 Terms owned by the PIMD engine  [code]

Let masses m, extended-temperature kinetic energy, bead-averaged potential V̄ = (1/P)Σ_j E_j.

| Quantity | Definition | Code |
|---|---|---|
| Quantum KE (printed) | Σ_modes ½ m v_nm² (this is the *extended-ensemble* KE; its expectation is 3/2·kB·P·T per Cartesian-triplet-per-mode) | `kinetic_energy_*` in `thermostat_base.f90` |
| Spring energy | Σ_{k≥1} ½ m ω_k² q_nm,k² | `compute_spring_energy_nm` |
| Free-ring energy | KE + spring | `compute_free_ring_energy_nm` (conserved in exact free propagation: 20 000 steps, 8.9e-13 eV drift) |
| Physical potential | V̄ = (1/P)Σ E_j | `summary%pimd_potential_eV` |
| Thermostat heat | Σ (KE_before − KE_after) | per thermostat |
| "conserved-like" | KE + V̄ + spring + heat | `nvt_conserved_like` |
| Temperature diagnostics | physical T_el = T_ext/P; T_ext; centroid vs internal split; per-mode | `write_temperature_diagnostics` |
| Bead RMS spread | √⟨|q_j − q̄|²⟩ | `write_step_outputs` |

**F5 – "conserved-like" is not a conserved quantity of the whole system.** It uses only the first (or "primary") quantum particle's KE and spring (`primary_quantum_slot`), excludes classical KE and all other electrons' KE/spring, and adds an extended-ensemble KE (∝P) to V̄. With 113 electrons only one ring contributes. The logged "conserved span" of 20 → 10 157 eV across P = 4 → 64 (`FAST_SPEED_CONVERGENCE`) shows it scales with P and size, and the run notes themselves caution it is only a regression diagnostic. It cannot certify dynamics. Also, the `electron_localization.out` writer is passed a literal `0.0_dp` for one argument; check what that field is meant to hold.

**No physical estimators beyond these**: no thermodynamic/virial energy estimator, no quantum KE estimator, no radius-of-gyration distribution, no centroid-virial pressure.

### 9.2 Terms owned by the eReaxFF oracle  [indirect: names from `ereaxff_energy_terms_export.f`; formulas from the JAX port]

The legacy COMMON block exports a 14-term vector (kcal/mol) plus `estrc` total; the wrapper maps:

| # | Name | Legacy var | Meaning |
|---|---|---|---|
| 1 | `bond` | `eb` | bond energy from corrected bond orders |
| 2 | `overcoordination` | `ea` | over/under-coordination (`ovcor`) |
| 3 | `lone_pair` | `elp` | lone-pair energy (+ C–C "C2 correction") |
| 4 | `molecule` | `emol` | molecular energy (parameter "not used") |
| 5 | `valence_penalty` | `ev + epen` | valence-angle + double-bond angle penalty |
| 6 | `angle_conjugation` | `ecoa` | valence-angle conjugation |
| 7 | `hydrogen_bond` | `ehb` | hydrogen bonds (none in Jamil ffield) |
| 8 | `torsion` | `et` | torsion |
| 9 | `torsion_conjugation` | `eco` | torsion conjugation |
| 10 | `vdw` | `ew` | shielded, tapered vdW (Morse-type) |
| 11 | `coulomb` | `ep` | shielded, tapered Coulomb with EEM charges |
| 12 | `charge` | `ech` | charge-equilibration (self/EEM) energy |
| 13 | `external_field` | `efi` | external field |
| 14 | `electron` | `eelect` | explicit electron pair term |

Formulas below are from the JAX re-derivation (`jax_force_kernel.py`, docstrings state they follow `poten.f`/`reac.f`) and reproduce the oracle to 3.7×10⁻⁷ eV/Å on sampled frames:

- **Shielded Coulomb (all particle pairs incl. El, charges from the charge solve):** E = C_c · q_i q_j / (r³ + γ_ij⁻³)^{1/3} × Tap(r), C_c = 332.0638 kcal·Å/(mol·e²), γ_ij from `gammaEEM`, Tap = 7th-order taper to the upper radius (10 Å in the ffield).
- **Shielded vdW:** E = p₂(h₁ − 2h₂), h₁ = exp[p₃(1 − r_w/p₁)], h₂ = exp[½p₃(1 − r_w/p₁)], r_w = (r^{vpar29} + γ_w)^{1/vpar29}, × Tap.
- **Explicit electron–core attraction (core C ↔ El):** E = −C_c · β_el · Z_core / r · erf(√(2α) r) × Tap, with β_el, α_g from the `El` ffield entry and Z_core the core charge.
- **Explicit electron–electron repulsion:** E = +C_c · √(β_i β_j) / r · erf(√(√(α_iα_j)) r) × Tap.
- **Electron valence sharing (links electrons into the covalent terms):** for each core atom i, eshare₂,i = Σ_{El e} exp(−vincr_i r²_ie); the corrected valence a_val,c = a_val ± eshare₂ and lone-pair count stlp_c = stlp + eshare₂ feed the lone-pair and over/under-coordination terms (`lonpar`, `ovcor`). This is how El positions change the bonded energy and forces of carbons.
- **Bond/valence/torsion terms** act among core atoms through bond orders (corrected for over-coordination); the El particles influence them only through the sharing above and the `electron_valence_sign` in the bond inputs.
- **Charge energy `ech`:** EEM/ACKS2 charge energy; for fixed MOLCHARGE it is whatever the legacy code returns — not reconstructed in JAX.

**Between particles vs between particle and electron** (summary):

| Interaction | Terms |
|---|---|
| C–C (same bead) | bond, over/under, lone pair, valence angle/penalty/conjugation, torsion/conjugation, vdW, Coulomb |
| C–El (same bead) | electron–core Gaussian attraction (`eelect`), shielded Coulomb (charges), shielded vdW bucket, valence sharing into lone-pair/overcoordination |
| El–El (same bead index) | explicit repulsion (`eelect`), shielded Coulomb/vdW buckets |
| El–El (different bead index) | none from eReaxFF; ring-polymer springs only (engine side) |

**Caveat (reviewer item).** The JAX port computes nonbonded/electron-sharing distances with **central-cell distances, not minimum image** (explicit comments: minimum image "over-counts" relative to the legacy Verlet list in the 17.2×17.1 Å cell, and breaks overcoordination-force parity). That matches the legacy oracle numerically but means the legacy code's periodic treatment of El interactions deserves a check against the physical intent of a periodic graphene sheet.

### 9.3 Ring-polymer coupling energy

The total quantity the integrator actually conserves in the NVE limit (if the scheme were consistent) is KE_ext + E_spring + V̄(R_C, {q_j}) (+ classical KE). Note the oracle returns *physical* E_j per bead; the spring energy is never part of the oracle.

---

## 10. The legacy-Fortran bridge in detail

### 10.1 Three backends (`ereaxff_oracle_adapter.f90`, 1048 lines)

1. **external_process**: for every call, create `adapter_call_NNNN/`, write `control` (axes rewritten from `CRYSTX`), `ffield`, `fort.20` (iopt), `fort.35`, a `charge_memory.in`, `geo`; run a shell wrapper under a Python timeout runner that launches `oracle_single_point_force_export` (and optionally `single_point_repeatability`); parse text output. ~51 s per 225-atom call in the notes.
2. **persistent_worker**: a single long-lived `persistent_oracle_worker` process.
3. **persistent_worker_pool** (current production): `oracle_max_workers` workers, each its own process with its own initialised legacy COMMON state, bead→worker by `call_slot`.

### 10.2 Worker protocol (`tools/persistent_oracle_worker.f90`)

- Start-up: stage `ffield`, `fort.20`, `fort.35`, `geo` (extract `CRYSTX`), rewrite `control` axes; `call ereaxff_force_oracle_init()` once; write `persistent_worker.ready`.
- Loop: poll for `persistent_request.in` (10 ms spin on `system_clock`, i.e. a CPU busy-wait); read request id, (step, bead, natom), stage label, then per-atom `index x y z charge solver_guess`; write geo/control/charge-trace into a per-call directory; `call ereaxff_force_eval(x, natom, charge_in, solver_in, energy, force, charge_out, solver_out, status)`; then `ereaxff_oracle_get_energy_terms`, `…get_bond_orders`, `…get_nonbonded_buckets`; write summary/forces/terms/orders/buckets files; publish `persistent_response.ready` atomically (`mv`); delete request. Shutdown via flag file.
- Hard limits: `natmax = 2600`, `nbondmax = 400000`, `nbucketmax = 256`.
- Adapter side waits with `wait_for_file` polling every 25 ms; **the same 60 s `persistent_worker_start_timeout_s` is used as the per-request timeout** (a 225-atom call took ~20 s per worker in the P=16 pool runs; a larger system on a slower node would hit it).

### 10.3 Dispatch (`bead_force_manager.evaluate_oracle_beads`)

`serial`, `openmp_all` (all beads in parallel; documented as poor for the external-process backend), `bounded_pool` (batches of `oracle_max_workers`). Per-call timing records and a batch log (`oracle_timing.log`) are written. Failures: `fail_fast` → `error stop`.

### 10.4 Risks

- Legacy code relies on **global COMMON state** and fixed filenames; isolation is achieved by one process per worker, one request at a time. Memory/disk: P=16 one step produced 715 MB of `adapter_call_*` directories.
- Per-call IPC is file-based text (formatted `es24.15`), not in-memory; the worker's `step_index` argument is hard-coded 0 in the request header from the adapter (stage label carries the info).
- The legacy `geo` file written by the worker carries positions at `f10.5`; harmless for the persistent path (positions go in binary) but worth noting for the external-process path.
- Build depends on absolute `V2_ROOT` and pre-built `.o` files that this machine no longer has.
- `ipi_socket/` is an 8-line stub: no i-PI client exists despite design docs.

### 10.5 Fast JAX backend (experimental, gated)

`force_provider.type: fast_jax_oracle` calls a JAX worker (`jax_live_worker.py`), evaluating all beads in one compiled batch; currently covers bond, shielded Coulomb/vdW, valence, penalty, angle conjugation, torsion/conjugation, lone-pair/overcoordination (incl. C2 correction), electron pair terms, with **fixed** charges only. Not ported: hydrogen-bond (absent from Jamil ffield), dynamic ACKS2 solve, hole particles (`Ho`). Parity: max |ΔF| 3.7×10⁻⁷ eV/Å on sampled P=16 frames (2 cases) and 6.6×10⁻⁷ eV/Å on 9 frames after fixes; reported speed ≈2.1 s per P=16 force stage vs ≈20 s for the 4-worker exact pool. The exact Fortran oracle remains the declared truth source.

---

## 11. Operational / reproducibility findings

- Tree size: 6.3 GB, 36 591 files, dominated by `.venv-fast/`, `build/` run outputs, and pycache; no `.gitignore` coverage of `projects/*/__pycache__` observed beyond root rules.
- Hard-coded paths: `/Users/yeatanir/PhD/PIMD/ereaxff_pimd/{V1,V2}` in Makefile, driver defaults (lines ~143–144), configs; two copies of the Jamil PDF (`124104_1_online.pdf`) in repo.
- Config parser is a YAML *subset* (no anchors, nested maps, arbitrary sequences); selectors support only a small grammar (one `and`).
- Restart files (text) preserve positions, normal modes, charges, solver guesses, classical velocities, heats, RNG state, per-quantum beads/masses; split-run equals continuous-run is tested (harmonic and mixed).
- ReaxFF-compatible outputs (`fort.7`, `fort.8`, `molfra.out`, `xmolout`) are written once at the end as **centroid-projected** classical views with minimum-image topology.

---

## 12. Validation evidence (as reported by the project) and what it does/doesn't show

**Tests (project-reported all PASS, `test_00`–`test_35`):** constants/units; normal-mode orthogonality/round-trip; free-ring energy conservation (20 000 steps, 8.9e-13 eV drift); harmonic finite-P variance convergence P=1…32 and BAOAB+PILE-L sampled variance within 2.1% at P=16; OU variance within 0.7%; PILE-L internal-mode temperatures within 2.3%; force scaling; charge memory; legacy-oracle smoke/repeatability/NaN rejection; P=1 reproduction through the force manager; mixed classical/PIMD group assignment; structural stress (P up to 128, 10 000-particle groups); restart continuity; geometry/user-geo/periodic loaders; oracle dispatch modes; persistent worker; JAX parity tests; classical-thermostat force-off (Berendsen/Langevin).

**What they do not cover:**
- No test samples a *known non-harmonic* or mixed (classical + ring) distribution through the production path with 1/P scaling (F1).
- No NVE conservation test of the coupled mixed system with real oracle forces.
- The `jamil_force_off` run (10 000 steps, P=16, dt=0.05) validates the thermostat/free ring only (physical T_el = 292.8 K vs 300 K; extended internal 4684 K vs 4800 K target).

**Logged production-like runs that need explaining** (`RUN_NOTES.md`, `FAST_SPEED_CONVERGENCE_2026-09-11.md`):

| Run | Observation |
|---|---|
| P=16, 1 step, dt 0.1, exact oracle | electron PILE heat +0.648 eV after one step; physical T_el 505 K |
| P=16, 100 steps, dt 0.05, Berendsen τ=100 fs | carbon T → **2139 K** in 5 fs (flagged THERMAL FAIL); after removing a duplicate kick: ~1156 K; Berendsen τ=10 fs → 1016 K |
| P=16, 100 steps, dt 0.01 | carbon T 399.7 K, "conserved span" 303.6 eV, status **PASS** |
| P=4→64, 10 steps, dt 0.05 | conserved-like span 19.9 → 10 157 eV |
| P=16 vs 32 vs 64 | bead RMS spread 0.051 / 0.058 / 0.068 Å — not converged in P |

Interpretation for a reviewer: a stable dt window of ≈0.01 fs with 1 amu pseudo-electrons whose ring modes oscillate at up to 6670 cm⁻¹ (ω_max dt = 0.063 at dt=0.05 for exact free propagation is not the limit), carbon heating that scales with timestep, and non-scaling "conserved" quantities are consistent with an integrator/force-consistency problem (F1, F3), not just a thermostat problem — the force-off thermostat tests pass, and JAX-vs-oracle force parity is ~10⁻⁷ eV/Å. The project's own conclusion ("timestep/integrator stability") points the same way. The 300 eV "PASS" shows the guard thresholds are loose.

---

## 13. Missing / not implemented (explicit list)

1. PILE-G (global SVR on centroid) and plain SVR; configurable λ for PILE.
2. A consistent, validated force-scaling/temperature convention for eReaxFF-driven rings (F1).
3. A 1-force-batch-per-step BAOAB (or OBABO) scheme for the coupled system; multiple-time-stepping.
4. Whole-system conserved-quantity diagnostic and an NVE/NVT energy-drift test with real forces.
5. Quantum estimators: thermodynamic/virial energy and quantum KE, radius of gyration, centroid-virial pressure; any exchange (bosonic/fermionic) statistics; ring-polymer contraction.
6. Barostat/NPT; COM-momentum removal; constraints; multiple-temperature-bath documentation for mixed systems.
7. Dynamic ACKS2 validated in PIMD; ACKS2 in the fast backend; charge tolerance enforcement/diagnostics.
8. i-PI socket client (stub only).
9. Hole (`Ho`) particles in the fast backend; hydrogen bonds in the fast backend.
10. Ring pre-equilibration protocol (rings start collapsed), production-quality RNG, in-memory (non-file) oracle IPC, timeout separation, shared unit constants.
11. Reproducible build that does not depend on an external, absolute-path legacy checkout (vendor/pin the legacy source or a library interface).

---

## 14. Suggested review questions

1. Is the electron-ring convention intended to be P·T (unscaled V per bead) or Tuckerman (V/P at T)? The code mixes them (F1).
2. For a mixed classical-carbon/ring-electron system, is separate-bath thermostatting (carbons at T, electron modes at P·T) acceptable, and should the carbon bath be canonical (Langevin/SVR) rather than Berendsen?
3. With fictitious 1 amu electrons, what physical observable is claimed? Static (equilibrium) distributions need a correct ensemble (F1); dynamical quantities (host hopping rates à la Jamil Fig. 3) are not meaningful for RPMD of pseudo-particles without calibration.
4. How does `icharg=8`/`MOLCHARGE` treat charges when beads sit on different atoms — are charges bead-independent by construction?
5. Is central-cell (non-minimum-image) treatment of El interactions in the legacy Verlet list physically intended for the periodic sheet?
6. What are the acceptance criteria (drift, temperature, bead-spread convergence in P) before long production runs?

---

## Addendum — base eReaxFF source found (added after initial report)

Location: `/Users/yeatanir/PhD/Literature/eReaxFF/Duin_group/Margaret/eReaxFF_src` (`reac.f`, `poten.f`, `cbka.blk`, `minres.f`, `cg.f`, `shanno.f`; README dated Aug 05 2016). The V2 oracle files (`ereaxff_force_oracle.f`, `cbka_v2.blk`, `reac_v2`, `poten_v2`, `legacy_link_stubs`) are **still not found**; this is the upstream code V2 was derived from (it also contains `cbka_pimd.blk`/`pimd_main_only.f`, an earlier PIMD prototype). Everything below is verified against this source but not against V2 itself, so any V2-specific edits (e.g. the `ereaxff_force_eval` wrapper, how `charge_in`/`solver_guess` are consumed) remain unverified.

**Charge equilibration (`reac.f: charges`, `chargesa`) — upgrades §8.2 from [indirect] to [code, base source]:**
- `icharg=8`: every atom's charge is set directly from `MOLCHARGE` (`ch(i) = vmcha/natoms-in-group`) and the routine returns. **No solver runs.** This is the Jamil production path (control file has `icharg 8`), so charges are fixed and geometry-independent there; bead-local vs centroid "ACKS2" is physically moot for that target.
- `icharg=5`: fixed charges. `icharg=7`: ACKS2 (`chargesa`). Other values select EEM variants (`1` Mortier, `3` shielded Mortier, `4` full system, `6` fragment-aware).
- `chargesa` with `El` present uses a **compact C-only ACKS2**: unknowns = 2·n_C + 1 + n_MOLCHARGE-groups (e.g. 122 for 60 C), El charges frozen at their `MOLCHARGE` values and their shielded Coulomb field folded into the C right-hand side; El–El block skipped. Matrix: EEM `2η` diagonal, shielded Coulomb `14.40/(r³+γ⁻³)^{1/3}` × taper, a "kinetic"/softness block (`bsoft`) coupling the effective-potential unknowns, and charge-conservation rows. Solved with **MINRES** (max 1000 iterations, tolerance `convg`), warm-started from `smortold` (this is the `solver_guess`, length up to `neem`). Warns only if `istop ≥ 5`.
- Consequence for PIMD: with electrons present, even true ACKS2 treats El charges as fixed, so charge response is carbon-only.
- Charge energy: `ech = 23.02·Σ(χq + ηq²)·q + 23.02·Σ q·effpot` (kcal/mol), added to the total only if `ncha2=1`.

**Energy (`reac.f` line ~5383):** `estrc = eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+epist + ncha2·ech + efi + eelect`. Note `epist` (piston) exists in the source but is not in the 14-term export (maps `ev+epen` jointly; `epist` dropped), so the exported terms sum to `estrc` only when `epist=0`.

**Electron terms (`poten.f: electron`) match the JAX port exactly:** electron–core `eelect = −332.0638·β_el·Z_core·erf(√(2α_g)·r)/r × taper`, with `Z_core = znucel`; cutoff `swb = vpar(38)` (a separate electron taper radius, `tap7th`), loop over the Verlet pair list with images (`ix,iy,iz`), skip `r<1e-10`. `elecpos` (reac.f:8779) gives the electron-sharing: for every `El`–non-`El` pair, `eshare = exp(−vincr·r²)`, then `avalc −= eshare` for C, O, H, F, N, Os and `+= eshare` for Li, Lv, P, B, Av, Al; `stlpc += eshare`. It uses `dista2` over all atoms (no cutoff), which is consistent with the JAX note that electron sharing uses central-cell distances. Sign per element is therefore element-specific, not a single `electron_valence_sign` assumption.

**Caveats that carry over to the report:** the periodic treatment differs between `electron` (Verlet list with images) and `elecpos` (`dista2`); whether `dista2` applies minimum image should be checked in `reac.f` before trusting El sharing in periodic cells.

---

# Part II — Literature review outcome, i-PI/LAMMPS source checks, and development strategy

(Added 2026-11-11. Sources checked: i-PI `master` shallow clone, HEAD `93162d2` dated 2026-10-24; the user's `Metal-Reax` repository, HEAD `735ad6a`; an externally generated literature review pasted by the user, whose citations I did not independently re-fetch.)

## II.1 Verified against i-PI source  [code]

| Item | Finding in i-PI source |
|---|---|
| Ring frequency | `omegan = temp · nbeads · kB / ħ` (`engine/normalmodes.py:461`) → ω_P = P·kT/ħ, with T the physical temperature |
| Dynamics temperature | `get_ntemp`: "P times the physical T" (`engine/motion/dynamics.py:235`); piped into the thermostat and barostat |
| Force on beads | `beads.p += forces.f · dt` and `nm.pnm += nm.fnm · dt/2` (`dynamics.py:441,580`) — **no 1/P factor** |
| Conserved quantity | `econs = vspring + kin + forces.pot (+ bias + thermostat/barostat terms)` where `forces.pot` is the **sum** over bead potentials (`ensembles.py:259`); the physical potential printed to users is divided by nbeads (`properties.py:374,392`) |

**Conclusion:** i-PI uses the "P·T" convention with unscaled bead potential. The production engine here uses ω_P = P·kT/ħ and T_ext = P·T **but** divides bead forces by P. That mixes conventions and F1 is real (the literature review reaches the same conclusion independently).

**Modern i-PI features the 2014/2018 papers do not describe, found in the current source:**
- Socket protocol with 12-byte headers (`drivers/f90/driver.f90`, `MSGLEN=12`), messages `STATUS/POSDATA/GETFORCE/FORCEREADY`.
- **Batched evaluation** (`batch_size`, "batched POSDATA/GETFORCE", `sockets.py:234,762`) and **consolidated messages** (one send per request; `select()` multiplexing of results).
- **Shared-memory transport** (`mode='shm'`) in the Fortran reference driver: bulk arrays via POSIX shared memory, only handshake on the socket.
- **MPI interface** (`ipi/interfaces/mpi.py`, `FFMPI`, tags HELLO/INIT/WORK/RESULT/EXTRA/EXIT) with driver group roots.
- Reference **Fortran client** (`drivers/f90/{driver,fsockets,fmpi}.f90`) — a direct template for wrapping a Fortran force engine.

These close the "protocol not verified" gap that the literature review reported (its i-PI sources were the 2014 and 2018 papers).

**Gap found:** I grepped i-PI for per-atom or per-species bead counts and found none; the code, docs and examples use a single `nbeads`. The only per-atom controls I saw are atom-subset masks (`activeatoms`) and spatially localized contraction (a force-evaluation device). **A mixed system with classical carbons and ring-polymer electrons is therefore not a native i-PI configuration**, as far as I could determine; confirm with the i-PI developers.

## II.2 Does the literature review help? Yes, with caveats

**Useful, and I checked the maths:**
- Equation (3) (free-ring radius of gyration, d·ħ²/(12 m kT)·(1−1/P²)) and equation (4) (harmonic finite-P variances) are correct; (4) is identical to what `harmonic_finite_p_bead_variance_A2` already computes in this repo.
- The two-convention table (A: i-PI, B: Tuckerman) is the right framing.
- Its mixed-system treatment is correct: carbon with propagation mass P·M, thermostat at P·T and force Σ_j F_j is exactly equivalent to physical mass M, thermostat at T and **mean** force. This is exactly what the current engine does for carbons (§II.3).
- The 1-force-batch-per-step O–B–A–B–O algorithm answers F3. Cayley/BCOCB integrators, PILE-G/SVR, centroid-virial estimator, and the validation ladder are all directly usable.
- Charge-handling guidance (per-bead history, solver-residual vs force-error, iEL-SCF, shadow-MD ACKS2) matches my F6 concerns.

**Caveats:**
- All numeric acceptance thresholds in its validation table are explicitly **proposed**, not published standards.
- The review's mixed-system recipes are **inferred**, not taken from a code that does it.
- It flags the "1 amu pseudo-electron mass" problem: ring stiffness depends on the mass used in the imaginary-time action, so a 1 amu fictitious mass defines a different quantum model from a real electron. This deserves an explicit modelling decision (see §II.6).
- Several 2026 references (LAMMPS `fix pimd/langevin`, shadow-ACKS2, etc.) were not re-fetched by me.

## II.3 Minimal fix for F1 (derived)

Keep what is already right and change one thing:
- **Carbons (classical):** physical mass M, thermostat at T, force = mean_j F_j. This is exactly equivalent to the i-PI convention (mass P·M, thermostat P·T, force Σ_j F_j). **Leave as is.**
- **Electron beads:** force = raw F_j (unscaled), thermostat at P·T, ω_P = P·kT/ħ. **Remove the 1/P on electron kicks** (`half_kick_from_summary`; keep P_max/P_active logic only if mixed bead counts are retained, in which case a per-particle consistency check is needed).
- **Diagnostics:** conserved quantity = P·KE_C + KE_e(ext) + E_spring + Σ_j E_j (+ heat), summed over *all* particles; report physical potential as (1/P)Σ_j E_j.
- **Regression test:** run the production `mixed_baoab_quantum_step` through `bead_force_manager` with a harmonic provider and compare to `harmonic_finite_p_bead_variance_A2` at P = 4, 8, 16, plus a coupled harmonic classical–ring pair. This test fails today if F1 is real.

## II.4 Metal-Reax assessment (`github.com/yeatanir/Metal-Reax`)

**What it is:** a validated ReaxFF implementation as a LAMMPS plugin (`pair_style reaxff/metal`, `fix qeq/reaxff/metal`), with CPU-64, CPU-32 and Apple-Metal GPU backends. It reproduces stock LAMMPS `stable_30Sep2026` to ~3e-15 relative in energy and 3e-12 kcal/mol/Å in forces (CPU-64), supports QEq, QEq-R, QTPIE and **ACKS2**, periodic images via LAMMPS ghost atoms, a 58-configuration frozen test fixture set, bitwise-repeatable FP32, and a documented numerical policy and decision register. GPL-2.0.

**Why it matters here:** it is a **properly validated, element-agnostic ReaxFF engine whose source you control**. The lost V2 oracle was a legacy wrapper around ~2016 base code; Metal-Reax is a cleaner, tested foundation.

**What eReaxFF needs that is not there (from my reading of README/FEATURE_MATRIX and the base source `reac.f`/`poten.f`):**
1. `electron` pair term (electron–core Gaussian attraction, electron–electron repulsion) with its own taper radius `vpar(38)` and `znucel/betael/alphag` parameters.
2. `elecpos` electron sharing (per-element sign table, `vincr`) feeding `avalc/stlpc` into lone-pair and over/under-coordination, plus the corresponding force terms.
3. `El`-aware charge handling: fixed electron charges (`MOLCHARGE`/`icharg 8`) and the compact carbon-only ACKS2 with electron field folded into the RHS.
4. ffield parsing of the extra `El` parameters (ADR-003 style compatibility flags).
5. Periodic treatment decision for electron terms (§9.2 caveat): Metal-Reax's ADR-004 forbids minimum image and counts each pair once via ghosts; legacy `elecpos` used all-atom `dista2` distances. Settle by checking `dista2` and comparing against the reference binary.
6. **Multi-geometry batching** (`backend.metal_batching` is "Planned", M8) — the one feature PIMD really needs.
7. Precision: all-FP32 GPU mode is noisier in NVE; use the mixed `bonded cpu64` mode or CPU-64 for PIMD validation.
8. Platform: Metal is Apple-only. The CPU-64 path is portable to Linux/HPC; a CUDA/Kokkos path is a separate project.

**Verification resources for the eReaxFF additions (all on disk):** the original source and compiled binaries in `…/Margaret/eReaxFF_src` (`ereac`, `reac.f`, `poten.f`), the stored 225-atom reference energy (−12206.99157847738 kcal/mol), the saved `adapter_call_*` force/charge outputs, and the `build/` oracle binaries.

## II.5 i-PI vs LAMMPS vs own engine — recommendation

Short answer: **build the eReaxFF force engine as a LAMMPS pair style first, and use i-PI as the path-integral driver (LAMMPS acts as the i-PI client).** Do not choose between them.

Reasons:
1. **Socket overhead does not matter at your force cost.** The literature numbers (LAMMPS in-process 3–12× faster than i-PI sockets) are for ML potentials with ms-scale force calls. eReaxFF calls take hundreds of ms to seconds even on fast backends, so socket/shm/MPI latency is negligible. i-PI is therefore not a performance liability for you.
2. **i-PI provides, validated, what the in-house engine is missing:** PILE-L/PILE-G, SVR, NHC, GLE/PIGLET, Cayley/BCOCB-class splitting, estimators (centroid-virial etc.), RPC/MTS, restart, and a regression-test culture. Writing and validating those yourself is the largest correctness risk in the current engine (F1–F5).
3. **A LAMMPS pair style gives three things at once:** (a) classical eReaxFF MD, which the Jamil benchmark needs and is the ground truth for validation against the original `ereac`; (b) i-PI PIMD through LAMMPS' i-PI client fix; (c) in-process LAMMPS PIMD (`fix pimd/*`) as an independent cross-check.
4. **Metal-Reax already gives a validated LAMMPS-integrated ReaxFF engine**; adding eReaxFF is an extension of known-validated code rather than a new engine.
5. Wrapping the Fortran oracle as a native i-PI client is possible (i-PI ships a Fortran driver) but would still need an eReaxFF engine you can build and validate, which is the same work as (4).

**Real limitation to resolve first:** both i-PI and LAMMPS PIMD treat *all* atoms of the simulated group as ring polymers (i-PI: confirmed by my search; LAMMPS: checked in §II.8, only partially group-aware). Your science needs classical carbons plus quantum electrons. Options:
- **Option 1 (simplest, physically more complete):** make carbons ring polymers too. Carbon NQEs are small but real at 300 K; cost is P× more atoms per replica, same force calls (each bead replica is one force call either way). The beads of C are cheap relative to the force call.
- **Option 2:** write an i-PI `Motion` class for mixed species (modular architecture makes this plausible but it is real development).
- **Option 3:** keep the in-house Fortran engine for mixed systems after fixing F1–F5, using the eReaxFF engine through the same client interface.

I recommend Option 1 first because it produces results fastest, then use the in-house engine as the independent check for the classical-carbon limit.

## II.6 Revised roadmap

1. **Decide the electron physical model** (mass in the imaginary-time action vs 1 amu fictitious mass), because it determines P and what observables are meaningful.
2. **Fix F1 in the current engine** (§II.3) and add the production-path harmonic test. Small change, immediate credibility.
3. **eReaxFF in Metal-Reax (CPU-64 first):** electron term, electron sharing, El charge handling, ffield parser. Validate against `ereac` on the Margaret test cases and the stored 225-atom reference.
4. **LAMMPS classical eReaxFF validation:** reproduce Jamil et al. classical results (Berendsen, 1 amu electrons) and check energy conservation in NVE.
5. **i-PI integration:** run P=1 (equals classical), then harmonic and anharmonic tests, then P-convergence of bead spread; add per-bead charge warm-start handling and use i-PI's batched or MPI interface for parallel beads.
6. **Performance:** multi-geometry batching in the engine (M8 item), then GPU or MPI scaling studies.
7. **Charge model:** decide icharg 8 (fixed) vs compact ACKS2; if ACKS2, per-bead history and residual-vs-force-error tests.

## II.7 Open items I could not settle

- LAMMPS `fix pimd/*` handling of subsets of atoms: **read, not run** (see §II.8).
- Whether i-PI has any undocumented per-species bead mechanism (a negative grep is not proof).
- Whether legacy `dista2` applies minimum image (needs a look at `reac.f`).
- Exact differences between the lost V2 sources and the base code; only reference outputs can settle it.
- Licensing compatibility (Metal-Reax GPL-2.0-only; i-PI is a separate-process client; LAMMPS GPL-2.0) — confirm before distribution.

## II.8 LAMMPS `fix pimd/*` and subsets of atoms (checked in source, LAMMPS HEAD `823eafc`, 2026-11-10)

Read from `src/REPLICA/fix_pimd_langevin.cpp` and `fix_pimd_nvt.cpp`; **not run**.

- **`fix pimd/nvt`:** no use of the group mask (`groupbit`) anywhere; it acts on all atoms. Not usable for a quantum subset.
- **`fix pimd/langevin`:** honours the group mask in many places: velocity kicks, spring forces, centroid step, thermostat, barostat sums and the primitive/centroid-virial estimators (`mask[i] & groupbit`, `group->count(igroup)`).
- **But the bead bookkeeping is not group-aware:** coordinate gathering across replicas (`collect_xc`, `MPI_Allgatherv` of all atoms, `ntotal = atom->natoms`) and the normal-mode transform (`nmpimd_transform`) loop over **all** atoms. So the partitions always exchange and transform every atom, even those outside the group.
- **Docs:** `doc/src/fix_pimd.rst` does not describe subset use. The group handling appears to exist to exclude virtual sites (the 2026 paper mentions Wannier centres), not to support classical atoms.
- **Consequence for mixed systems:** atoms outside the group would be integrated by another fix **independently in each partition, with that partition's own bead force**, so classical carbons would diverge between replicas. A correct mixed system needs an extra step that averages forces (or synchronises positions) across partitions for the classical atoms; I found no such fix in the clone I read (sparse checkout of `REPLICA`, `REAXFF`, `EXTRA-FIX`, docs only; other packages not searched).
- **Verdict:** LAMMPS does not natively support classical atoms with ring-polymer electrons in one run either. Neither LAMMPS nor i-PI removes the need for a small, well-tested modification, so the question is where to put it. A cross-partition force-averaging fix for classical atoms is small; an i-PI per-atom classical mask is larger but benefits from i-PI's thermostats and estimators.

---

## Appendix A — Key code references

| Topic | File:line |
|---|---|
| ω_P, ω_k | `src/pimd/normal_modes.f90:75,88` |
| C matrix | `src/pimd/normal_modes.f90` (`build_cmat`) |
| Exact free propagation | `src/pimd/free_ring_propagator.f90` |
| Spring / free-ring energy | `src/pimd/spring_energy.f90` |
| PILE-L coefficients | `src/thermostat/thermostat_pile_l.f90:46–48` |
| PILE-L init at T·P | `src/drivers/ereaxff_mixed_baoab_main.f90:952` |
| Quantum velocity init at kT·P | `…main.f90:635` |
| Full step | `…main.f90:975` (`mixed_group_baoab_step`) |
| Quantum BAOAB step, force kicks | `src/integrator/mixed_baoab_integrator.f90:30–170` |
| 1/P force scaling | `src/force/bead_force_manager.f90:82,241,278,336,576`; rescale `mixed_baoab_integrator.f90:164` |
| Centroid-ACKS2 branch | `bead_force_manager.f90:534–584` |
| Classical Langevin / Berendsen | `…main.f90:1115`, `1190` |
| "conserved-like" | `src/output/diagnostics.f90`; `…main.f90:719` |
| Unit conversions | `src/force/ereaxff_oracle_adapter.f90:10–11`; `src/core/units.f90` |
| Worker protocol | `tools/persistent_oracle_worker.f90` |
| 14 energy terms | `tools/ereaxff_energy_terms_export.f` |
| Electron-pair / sharing / lone-pair formulas | `projects/jamil_replication/fast_backend/jax_force_kernel.py:480–565, 1502–1645` |

## Appendix B — Representative production config (Jamil P=16)

```yaml
simulation:      {dt_fs: 0.1, steps: 1, seed: 13579, rescale_initial_pimd_velocities: true}
force_provider:  {type: ereaxff_oracle, oracle_backend_mode: persistent_worker_pool,
                  oracle_dispatch_mode: bounded_pool, oracle_max_workers: 4, repeatability_check: false}
ring_polymer:    {convention: standard_ipi, default_beads: 1}      # per-group beads override
charge_coupling: {mode: bead_local_acks2, charge_tolerance: 1.0e-8}
geometry:        {template: user_geo_structure, periodic: true, crystx: [17.22871, 17.05200, 20.0, 90, 90, 90]}
groups:
  - {name: graphene_carbons, selector: "element == C",  dynamics: classical, thermostat: berendsen, temperature_K: 300, tau_fs: 100}
  - {name: electrons,        selector: "element == El", dynamics: pimd, beads: 16, thermostat: pile_l,
     temperature_K: 300, centroid_gamma_fs_inv: 0.001, mass_amu: 1.0}
```

(112 C classical, 113 El each a 16-bead ring; 225 atoms per oracle call.)
