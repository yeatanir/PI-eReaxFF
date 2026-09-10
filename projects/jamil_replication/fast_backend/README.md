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
the exact oracle dataset has to exist before a JAX/native force kernel can be
accepted as a dynamics backend.
