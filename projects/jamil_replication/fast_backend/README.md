# Fast Backend Workspace

This directory contains the acceleration path for the Jamil 112C/113El
eReaxFF/RPMD target.

The exact legacy oracle remains the reference backend.  The fast backend must
pass oracle-dataset parity before it is used for production dynamics.

## Current Contract

- `force_provider.type: ereaxff_oracle` remains the default exact path.
- `force_provider.type: fast_jax_oracle` launches the persistent JAX worker
  from the mixed BAOAB driver when `fast_backend_allow_unvalidated: true` is
  set explicitly.
- The connected fast path is currently guarded to the Jamil electron-only
  C/El, fixed-`MOLCHARGE`, orthorhombic-cell target.  It batches all beads in
  one JAX force call and writes the same `.xyz`, energy, charge, restart, and
  compatibility outputs as the oracle driver path.
- `build_oracle_dataset.py` converts existing restart/timing outputs into
  JSONL reference cases for fast-backend development.

## Dataset Builder

```sh
python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py \
  --run-dir build/jamil_replication/verify_p4_pool_steps1/out \
  --output build/jamil_replication/fast_backend/oracle_cases.jsonl
```

Prefer the adapter-work mode when exact oracle call directories are available:

```sh
python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py \
  --adapter-work build/jamil_replication/verify_p4_pool_steps1/out/adapter_work \
  --output build/jamil_replication/fast_backend/oracle_labeled_p4.jsonl
```

The generated JSONL is intentionally kept under `build/` and out of git.  Each
case records positions, bead-local charge memory, solver guesses, and run-level
metadata.  Adapter-work extraction also records exact energy, force, charge,
solver, status, stage, bead, and timing labels from the legacy oracle outputs.
For the current Jamil target it also records the legacy `icharg=8`
`MOLCHARGE` constraints used by the first JAX charge kernel.
Legacy force files are written by the V2 oracle wrapper in MD acceleration
units.  The dataset builder converts them back to true `eV/A` at the JSONL
boundary so JAX autodiff forces and oracle labels use the same units.
Repo-built oracle exporters write `single_point_energy_terms.out` in every
call directory.  The dataset builder records those per-term energy labels, or
falls back to a retained legacy `fort.73` layout when that is the only term
output available.  The term keys follow the Jamil `reac.f` output order:
`bond`, `overcoordination`, `lone_pair`, `molecule`, `valence_penalty`,
`angle_conjugation`, `hydrogen_bond`, `torsion`, `torsion_conjugation`, `vdw`,
`coulomb`, `charge`, `external_field`, and `electron`.

Repo-built oracle exporters also write `single_point_bond_orders.out` with
legacy `abo`, `bo`, `bopi`, and `bopi2` values.  The dataset builder records
those labels when present.

Repo-built oracle exporters also write `single_point_nonbonded_buckets.out`
with legacy `nonbonnoinS` vdW/Coulomb sums grouped by atom-type pair.  The
dataset builder records those buckets when present so JAX nonbonded parity can
be debugged by interaction class instead of only by total energy.

Validate the generated JSONL line-by-line:

```sh
python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py \
  --check build/jamil_replication/fast_backend/oracle_cases.jsonl
```

Do not use `python3 -m json.tool` for this file; it is JSONL, with one case per
line, not one large JSON document.

Check the optional JAX entrypoint:

```sh
python3 projects/jamil_replication/fast_backend/jax_backend.py \
  --self-test \
  --dataset build/jamil_replication/fast_backend/oracle_cases.jsonl
```

Run the first JAX-side labeled readiness gate:

```sh
.venv-fast/bin/python projects/jamil_replication/fast_backend/jax_backend.py \
  --parity-check \
  --dataset build/jamil_replication/fast_backend/oracle_labeled_p4.jsonl
```

Compare unvalidated JAX energy contributions against labeled legacy terms:

```sh
.venv-fast/bin/python projects/jamil_replication/fast_backend/compare_energy_terms.py \
  --dataset build/jamil_replication/fast_backend/oracle_labeled_p4_terms_latest.jsonl \
  --output build/jamil_replication/fast_backend/term_compare_p4_latest.json
```

Compare JAX corrected bond-order state against legacy labels:

```sh
.venv-fast/bin/python projects/jamil_replication/fast_backend/compare_bond_orders.py \
  --dataset build/jamil_replication/fast_backend/oracle_labeled_p4_bond_orders_latest.jsonl \
  --output build/jamil_replication/fast_backend/bond_order_compare_p4_latest.json
```

Compare JAX nonbonded buckets against legacy labels:

```sh
.venv-fast/bin/python projects/jamil_replication/fast_backend/compare_nonbonded_buckets.py \
  --dataset build/jamil_replication/fast_backend/oracle_labeled_p4_nonbonded_latest.jsonl \
  --output build/jamil_replication/fast_backend/nonbonded_bucket_compare_p4_latest.json
```

Compare JAX force contributions and total forces against corrected `eV/A`
legacy labels:

```sh
.venv-fast/bin/python projects/jamil_replication/fast_backend/compare_force_terms.py \
  --dataset build/jamil_replication/fast_backend/oracle_labeled_p4_nonbonded_latest.jsonl \
  --output build/jamil_replication/fast_backend/force_term_report_p4_nonbonded.json
```

## Charge Kernel Status

The current Jamil target uses legacy `icharg=8` with explicit `MOLCHARGE`
constraints.  `jax_charge_kernel.py` implements that fixed-charge path and the
backend can evaluate it in batches with:

```json
{"command": "charge_batch", "cases": [ ... ]}
```

Full dynamic ACKS2 matrix construction/solve is still intentionally
unimplemented and will fail closed until it has oracle labels that exercise
nonzero solver state.

The mixed driver can exercise the current Jamil fast path with:

```yaml
force_provider:
  type: fast_jax_oracle
  fast_backend_command: ".venv-fast/bin/python projects/jamil_replication/fast_backend/jax_live_worker.py"
  fast_backend_model: "reaxff_corrected_bond_nonbonded_valence"
  fast_backend_allow_unvalidated: true
```

Run the end-to-end smoke test with:

```sh
make test_32_fast_jax_driver_smoke
```

Run the short P=16 stability regression with:

```sh
make test_33_fast_jax_stability
```

Run the heavier stability matrix manually with:

```sh
JAMIL_FAST_STABILITY_RUN_LONG=1 \
  projects/jamil_replication/run_fast_stability_check.sh
```

As of 2026-09-10, 10-step P=16 fast runs pass finite/output/shutdown and
thermal guards at `dt_fs=0.05` and `0.1`.  A 100-step `dt_fs=0.05` run remains
thermally unacceptable because carbon temperature rises to about `2139 K` over
5 fs; use that as the next off-equilibrium parity/profiling target before
removing the unvalidated gate.

Run selected off-equilibrium frame parity against the exact oracle with:

```sh
projects/jamil_replication/run_fast_offequilibrium_parity.sh
```

The regression target checks a small near-equilibrium subset:

```sh
make test_34_fast_offequilibrium_parity
```

As of 2026-09-10, the short-run subset (`step=0,10`, `bead=1`) matches the
oracle to about `3.7e-7 eV/A`.  The broader thermally bad 100-step diagnostic
subset (`step=0,50,100`, `beads=1,8,16`) now matches to `6.6e-7 eV/A` after
matching the legacy central-cell electron-sharing convention and adding the
legacy C-C C2 correction that is accumulated in the `lone_pair` energy bucket.

## Force Kernel Status

`jax_force_kernel.py` currently provides these gated paths:

- strict oracle-label replay for protocol/parity testing;
- an experimental fixed-charge Coulomb-only JAX kernel for testing real
  autodiff, batching, and timing mechanics.
- an experimental source-derived eReaxFF shielded/tapered Coulomb slice that
  parses `gammaEEM` and taper radii from each oracle call directory's `ffield`.
- an experimental source-derived `nonbonnoinS` slice for the Jamil force field,
  combining shielded Morse vdW with shielded/tapered Coulomb.
- an experimental source-derived raw bond-order and uncorrected covalent
  bond-energy slice from `reac.f:srtbon`/`poten.f:covbon`.
- an experimental source-derived corrected bond-order and covalent bond-energy
  slice that applies the Jamil `poten.f:boncor` overcoordination and 1-3
  bond-order corrections before `covbon`.
- an experimental aggregate slice that sums corrected covalent bond energy,
  no-inner-wall vdW, shielded/tapered Coulomb, sparse valence/torsion terms,
  and electron-only `El` pair terms in one JAX force call.

The ReaxFF shielded Coulomb and vdW kernels use central-cell pair distances for
the current Jamil oracle-parity target.  That matches the legacy V2
`nonbonnoinS` bucket labels for the checked 112C/113El periodic cases; expanding
all 27 neighbor images overcounts relative to the legacy Verlet list.  The
`reaxff_shielded_coulomb` path follows the Jamil source `poten.f` Coulomb
expression and `reac.f` seventh-order taper polynomial.  The
`reaxff_nonbonded_no_inner_shielded` path follows the same `poten.f`
`nonbonnoinS` routine for the no-inner-wall, shielded-vdW case selected by the
current Jamil C/El parameters, including off-diagonal vdW pair overrides from
`ffield`.  The `reaxff_uncorrected_bond` path builds raw sigma/pi/pi-pi bond
orders and the direct covalent bond term.  The
`reaxff_corrected_bond` path applies the source-derived `boncor` bond-order
attenuation before the covalent bond-energy expression.  The
`reaxff_corrected_bond_nonbonded` path combines the implemented corrected bond,
vdW, and Coulomb slices.  The
`reaxff_corrected_bond_nonbonded_valence` path additionally reports valence,
penalty, angle-conjugation, torsion, torsion-conjugation, `lone_pair`,
`overcoordination`, `electron_core`, and `electron_electron` contributions.  It
is electron-only for the current Jamil target and fails fast on `Ho`/hole
particles.  Corrected covalent bond orders, bond energy, valence/penalty,
torsion, torsion-conjugation, lone-pair, electron-pair terms, and nonbonded
vdW/Coulomb buckets now have direct oracle labels for comparison.  On the
current P=4 nonbonded Jamil oracle dataset, corrected BO topology has zero
missing/extra bonds, scalar total-energy error is below `1e-5 eV`, and total
force error is below `4e-7 eV/A` after converting legacy oracle forces to true
`eV/A`.  Full dynamic ACKS2 charge solving is still intentionally unimplemented
and will fail closed until it has oracle labels that exercise nonzero solver
state.

Force evaluation fails closed by default:

```json
{"command": "force_batch", "cases": [ ... ]}
```

For protocol/parity tests only, replay can be enabled explicitly:

```json
{"command": "force_batch", "allow_oracle_replay": true, "cases": [ ... ]}
```

For kernel-development tests only, the Coulomb slice can be selected through
the unvalidated CLI gate:

```sh
.venv-fast/bin/python projects/jamil_replication/fast_backend/jax_backend.py \
  --allow-unvalidated
```

with a stdin request:

```json
{"command": "force_batch", "force_kernel": "fixed_charge_coulomb", "cases": [ ... ]}
```

or:

```json
{"command": "force_batch", "force_kernel": "reaxff_shielded_coulomb", "cases": [ ... ]}
```

or:

```json
{"command": "force_batch", "force_kernel": "reaxff_nonbonded_no_inner_shielded", "cases": [ ... ]}
```

or:

```json
{"command": "force_batch", "force_kernel": "reaxff_uncorrected_bond", "cases": [ ... ]}
```

or:

```json
{"command": "force_batch", "force_kernel": "reaxff_corrected_bond", "cases": [ ... ]}
```

or:

```json
{"command": "force_batch", "force_kernel": "reaxff_corrected_bond_nonbonded", "cases": [ ... ]}
```
