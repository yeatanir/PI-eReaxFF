#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE_CFG="$ROOT_DIR/projects/jamil_replication/config_jamil_112c_113el_p16_300K_prod.yaml"
RUN_ROOT="${JAMIL_FORCE_OFF_ROOT:-$ROOT_DIR/build/jamil_replication/force_off_$(date +%Y%m%d_%H%M%S)_$$}"
SUMMARY="$ROOT_DIR/projects/jamil_replication/jamil_force_off_summary.out"
P="${JAMIL_FORCE_OFF_P:-16}"
STEPS="${JAMIL_FORCE_OFF_STEPS:-100}"
DT_FS="${JAMIL_FORCE_OFF_DT_FS:-0.05}"
GAMMA="${JAMIL_FORCE_OFF_GAMMA:-0.001}"
MASS_AMU="${JAMIL_FORCE_OFF_MASS_AMU:-1.0}"

if [[ ! -f "$BASE_CFG" ]]; then
  echo "run_jamil_force_off_check: missing base config: $BASE_CFG" >&2
  exit 2
fi

mkdir -p "$RUN_ROOT"
CFG="$RUN_ROOT/config_force_off.yaml"
OUT="$RUN_ROOT/out"
mkdir -p "$OUT"

python3 - "$BASE_CFG" "$CFG" "$OUT" "$P" "$STEPS" "$DT_FS" "$GAMMA" "$MASS_AMU" <<'PY'
from __future__ import annotations

import pathlib
import re
import sys

src = pathlib.Path(sys.argv[1])
dst = pathlib.Path(sys.argv[2])
out = sys.argv[3]
p = sys.argv[4]
steps = sys.argv[5]
dt_fs = sys.argv[6]
gamma = sys.argv[7]
mass_amu = sys.argv[8]

text = src.read_text(encoding="utf-8")

def replace(pattern: str, value: str, text_in: str) -> str:
    return re.sub(pattern, value, text_in, flags=re.M)

text = replace(r"^  dt_fs:.*$", f"  dt_fs: {dt_fs}", text)
text = replace(r"^  steps:.*$", f"  steps: {steps}", text)
text = replace(r"^  output_stride:.*$", "  output_stride: 1", text)
text = replace(r"^  restart_stride:.*$", "  restart_stride: 0", text)
if "rescale_initial_pimd_velocities:" in text:
    text = replace(r"^  rescale_initial_pimd_velocities:.*$", "  rescale_initial_pimd_velocities: true", text)
else:
    text = replace(r"^  restart_stride:.*$", "  restart_stride: 0\n  rescale_initial_pimd_velocities: true", text)
text = replace(r"^  output_dir:.*$", f"  output_dir: \"{out}\"", text)
text = replace(r"^  type:.*$", "  type: dummy", text)
text = replace(r"^  default_beads:.*$", f"  default_beads: {p}", text)
text = replace(r"^    beads:.*$", f"    beads: {p}", text)
text = replace(r"^    centroid_gamma_fs_inv:.*$", f"    centroid_gamma_fs_inv: {gamma}", text)
text = replace(r"^    mass_amu:.*$", f"    mass_amu: {mass_amu}", text)
dst.write_text(text, encoding="utf-8")
PY

START_NS="$(python3 -c 'import time; print(time.time())')"
"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$CFG" "$OUT" > "$RUN_ROOT/stdout.log"
END_NS="$(python3 -c 'import time; print(time.time())')"
WALL_S="$(python3 - "$START_NS" "$END_NS" <<'PY'
import sys
print(float(sys.argv[2]) - float(sys.argv[1]))
PY
)"

python3 - "$OUT" "$SUMMARY" "$RUN_ROOT" "$WALL_S" "$P" "$STEPS" "$DT_FS" "$GAMMA" "$MASS_AMU" <<'PY'
from __future__ import annotations

import pathlib
import sys

out = pathlib.Path(sys.argv[1])
summary = pathlib.Path(sys.argv[2])
run_root = pathlib.Path(sys.argv[3])
wall_s = float(sys.argv[4])
p = int(sys.argv[5])
steps = int(sys.argv[6])
dt_fs = float(sys.argv[7])
gamma = float(sys.argv[8])
mass_amu = float(sys.argv[9])

def last_row(path: pathlib.Path) -> list[float]:
    row: list[float] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        row = [float(x) for x in line.split()]
    if row is None:
        raise SystemExit(f"missing rows in {path}")
    return row

electron = last_row(out / "electron_temperature.out")
partition = last_row(out / "electron_mode_partition.out")
spread = last_row(out / "bead_spread.out")
heat = last_row(out / "thermostat_heat.out")

lines = [
    "# Jamil Force-Off PIMD Summary",
    "",
    f"run_root={run_root}",
    f"P={p}",
    f"steps={steps}",
    f"dt_fs={dt_fs}",
    f"gamma_fs_inv={gamma}",
    f"mass_amu={mass_amu}",
    f"wall_s={wall_s:.3f}",
    f"physical_electron_temperature_K={electron[1]:.6f}",
    f"all_extended_mode_temperature_K={partition[3]:.6f}",
    f"expected_extended_mode_temperature_K={partition[5]:.6f}",
    f"centroid_extended_mode_temperature_K={partition[1]:.6f}",
    f"internal_extended_mode_temperature_K={partition[2]:.6f}",
    f"centroid_ke_fraction={partition[9]:.8f}",
    f"quantum_heat_eV={heat[3]:.12f}",
    f"bead_rms_spread_A={spread[4]:.12f}",
]
summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
PY
