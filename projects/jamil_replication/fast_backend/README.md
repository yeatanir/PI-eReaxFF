# Fast Backend Workspace

This directory contains the acceleration path for the Jamil 112C/113El
eReaxFF/RPMD target.

The exact legacy oracle remains the reference backend.  The fast backend must
pass oracle-dataset parity before it is used for production dynamics.

## Current Contract

- `force_provider.type: ereaxff_oracle` remains the default exact path.
- `force_provider.type: fast_jax_oracle` is accepted by config parsing but the
  mixed driver stops immediately until the backend is validated and connected.
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

The entrypoint currently fails closed for force evaluation.  That is deliberate:
the exact oracle dataset and charge parity have to exist before a JAX/native
force kernel can be accepted as a dynamics backend.

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
vdW/Coulomb buckets now have direct oracle labels for comparison.  Full dynamic
ACKS2 charge solving, the remaining overcoordination scalar mismatch, and other
total-energy terms are still unvalidated.  These kernels must not be used for
production dynamics or total-force parity.

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
