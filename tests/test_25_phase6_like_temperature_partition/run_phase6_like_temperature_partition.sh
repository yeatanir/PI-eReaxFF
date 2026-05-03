#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_25_phase6_like_temperature_partition.summary"
WORK_DIR="$ROOT_DIR/build/test_25_phase6_like_temperature_partition"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

cfg="$WORK_DIR/phase6_like_config.yaml"
out="$WORK_DIR/out"
mkdir -p "$out"

python3 - "$ROOT_DIR/config/example_user_geo_structure.yaml" "$cfg" "$out" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
out_dir = sys.argv[3]
source = source.replace("  steps: 100", "  steps: 20")
source = source.replace("  output_stride: 10", "  output_stride: 5")
source = source.replace("  restart_stride: 20", "  restart_stride: 10")
source = source.replace("  output_dir: \"\"", f"  output_dir: \"{out_dir}\"")
target.write_text(source, encoding="utf-8")
PY

"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$cfg" "$out" > "$out/stdout.log"

for file in carbon_temperature.out electron_temperature.out electron_mode_temperature.out thermostat_heat.out fort.7 fort.8 molfra.out centroid.xyz restart.dat; do
  if [[ ! -s "$out/$file" ]]; then
    echo "test_25_phase6_like_temperature_partition: FAIL missing $file" | tee "$SUMMARY"
    exit 1
  fi
done

if ! grep -q "connection table" "$out/fort.7"; then
  echo "test_25_phase6_like_temperature_partition: FAIL fort.7 missing connection-table header" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "connection table" "$out/fort.8"; then
  echo "test_25_phase6_like_temperature_partition: FAIL fort.8 missing connection-table header" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "molecule list per frame" "$out/molfra.out"; then
  echo "test_25_phase6_like_temperature_partition: FAIL molfra.out missing molecule-list header" | tee "$SUMMARY"
  exit 1
fi

carbon_last=$(awk 'NF && $1 !~ /^#/ {v=$2} END {print v}' "$out/carbon_temperature.out")
electron_last=$(awk 'NF && $1 !~ /^#/ {v=$2} END {print v}' "$out/electron_temperature.out")
mode_lines=$(grep -vc '^#' "$out/electron_mode_temperature.out")
heat_last=$(awk 'NF && $1 !~ /^#/ {v=$2} END {print v}' "$out/thermostat_heat.out")

python3 - "$carbon_last" "$electron_last" "$heat_last" <<'PY'
import math
import sys

for idx, value in enumerate(sys.argv[1:], start=1):
    try:
        x = float(value)
    except Exception as exc:
        raise SystemExit(f"non-numeric diagnostic {idx}: {value}") from exc
    if not math.isfinite(x):
        raise SystemExit(f"non-finite diagnostic {idx}: {value}")
PY

{
  echo "test_25_phase6_like_temperature_partition: PASS"
  echo "  carbon_temperature_K = $carbon_last"
  echo "  electron_temperature_K = $electron_last"
  echo "  thermostat_heat_eV = $heat_last"
  echo "  electron_mode_rows = $mode_lines"
  echo "  hot_electron_status = diagnostic_only"
} | tee "$SUMMARY"
