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

`jax_force_kernel.py` currently provides two gated paths:

- strict oracle-label replay for protocol/parity testing;
- an experimental fixed-charge Coulomb-only JAX kernel for testing real
  autodiff, batching, and timing mechanics.
- an experimental source-derived eReaxFF shielded/tapered Coulomb slice that
  parses `gammaEEM` and taper radii from each oracle call directory's `ffield`.
- an experimental source-derived `nonbonnoinS` slice for the Jamil force field,
  combining shielded Morse vdW with shielded/tapered Coulomb.
- an experimental source-derived raw bond-order and uncorrected covalent
  bond-energy slice from `reac.f:srtbon`/`poten.f:covbon`.

The Coulomb kernels are orthorhombic-minimum-image only.  The
`reaxff_shielded_coulomb` path follows the Jamil source `poten.f` Coulomb
expression and `reac.f` seventh-order taper polynomial.  The
`reaxff_nonbonded_no_inner_shielded` path follows the same `poten.f`
`nonbonnoinS` routine for the no-inner-wall, shielded-vdW case selected by the
current Jamil C/El parameters.  The `reaxff_uncorrected_bond` path builds raw
sigma/pi/pi-pi bond orders and the direct covalent bond term, but does not yet
apply `boncor` overcoordination or 1-3 bond-order corrections.  These kernels
still omit valence, torsion, lone-pair, over/undercoordination energy,
ACKS2 dynamic charge solving, self-image details beyond minimum image, and
other total-energy terms.  They must not be used for production dynamics or
total-force parity.

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
