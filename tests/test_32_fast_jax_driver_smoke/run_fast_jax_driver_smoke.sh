#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

if [[ ! -x ".venv-fast/bin/python" ]]; then
  echo "test_32_fast_jax_driver_smoke: missing .venv-fast/bin/python; run the JAX install step first" >&2
  exit 77
fi

WORK_DIR="build/test_32_fast_jax_driver_smoke"
OUT_DIR="$WORK_DIR/output"
CONFIG="$WORK_DIR/config_fast_jax.yaml"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"

python3 - <<'PY'
from pathlib import Path

src = Path("projects/jamil_replication/config_jamil_112c_113el_p1_300K_prod.yaml")
dst = Path("build/test_32_fast_jax_driver_smoke/config_fast_jax.yaml")
text = src.read_text()
text = text.replace("  type: ereaxff_oracle\n", "  type: fast_jax_oracle\n")
text = text.replace("  steps: 1\n", "  steps: 1\n")
text = text.replace("  oracle_dispatch_mode: serial\n", "  oracle_dispatch_mode: serial\n")
text = text.replace(
    "  repeatability_check: false\n",
    "  repeatability_check: false\n"
    '  fast_backend_command: ".venv-fast/bin/python projects/jamil_replication/fast_backend/jax_live_worker.py"\n'
    '  fast_backend_model: "reaxff_corrected_bond_nonbonded_valence"\n'
    "  fast_backend_allow_unvalidated: true\n",
)
dst.write_text(text)
PY

build/ereaxff_mixed_baoab_main "$CONFIG" "$OUT_DIR" > "$WORK_DIR/driver.out" 2>&1

grep -q "force provider: fast backend" "$WORK_DIR/driver.out"
test -s "$OUT_DIR/pimd_energy.out"
test -s "$OUT_DIR/centroid.xyz"
test -s "$OUT_DIR/trajectory_centroid.xyz"
test -s "$OUT_DIR/trajectory_beads.xyz"
test -f "$OUT_DIR/fast_jax_work/fast_jax_worker.log"
test -s "$OUT_DIR/fast_jax_work/response.out"
grep -q "status 0" "$OUT_DIR/fast_jax_work/response.out"

echo "test_32_fast_jax_driver_smoke: PASS"
