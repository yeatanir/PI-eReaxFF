#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_23_ereaxff_mixed_baoab_geometry_matrix.summary"
WORK_DIR="$ROOT_DIR/build/test_23_ereaxff_mixed_baoab_geometry_matrix"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

positions=(
  "7.38300,5.68433,10.00100"
  "6.85000,4.95000,10.06000"
)

for pos in "${positions[@]}"; do
  cfg="$WORK_DIR/config_${pos//,/__}.yaml"
  out="$WORK_DIR/out_${pos//,/__}"
  mkdir -p "$out"
  python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$cfg" "$out" "$pos" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
out_dir = sys.argv[3]
pos = sys.argv[4]
source = source.replace("  steps: 1000", "  steps: 4")
source = source.replace("  output_stride: 10", "  output_stride: 2")
source = source.replace("  restart_stride: 100", "  restart_stride: 2")
source = source.replace("  output_dir: \"\"", f"  output_dir: \"{out_dir}\"")
source = source.replace("  electron_position_A: [7.38300, 5.68433, 10.00100]", f"  electron_position_A: [{pos}]")
target.write_text(source, encoding="utf-8")
PY
  "$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$cfg" "$out" > "$out/stdout.log"
  for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out electron_localization.out centroid.xyz fort.7 fort.8 molfra.out xmolout restart.dat; do
    if [[ ! -s "$out/$file" ]]; then
      echo "test_23_ereaxff_mixed_baoab_geometry_matrix: FAIL missing $file for pos=$pos" | tee "$SUMMARY"
      exit 1
    fi
  done
  python3 - "$out/stdout.log" "$pos" <<'PY'
import re
import sys

log_path, expected = sys.argv[1], [float(x) for x in sys.argv[2].split(",")]
line = None
with open(log_path, "r", encoding="utf-8") as handle:
    for raw in handle:
      if raw.startswith("electron_start_A:"):
        line = raw
        break
if line is None:
    raise SystemExit("missing electron_start_A line")
values = [float(x) for x in re.findall(r"[-+]?\d+\.\d+E[+-]\d+", line.upper())]
if len(values) != 3:
    raise SystemExit(f"unexpected geometry line: {line.strip()}")
for got, want in zip(values, expected):
    if abs(got - want) > 1.0e-4:
        raise SystemExit(f"geometry mismatch: got {values}, want {expected}")
PY
done

bad_template_cfg="$WORK_DIR/bad_template.yaml"
python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$bad_template_cfg" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
source = source.replace("  template: phase4_fixed_graphene", "  template: unsupported_template")
target.write_text(source, encoding="utf-8")
PY

set +e
"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$bad_template_cfg" "$WORK_DIR/bad_template_out" > "$WORK_DIR/bad_template.log" 2>&1
bad_template_status=$?
set -e
if [[ "$bad_template_status" -eq 0 ]]; then
  echo "test_23_ereaxff_mixed_baoab_geometry_matrix: FAIL unsupported geometry template unexpectedly succeeded" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "unsupported geometry template" "$WORK_DIR/bad_template.log"; then
  echo "test_23_ereaxff_mixed_baoab_geometry_matrix: FAIL expected template error not found" | tee "$SUMMARY"
  cat "$WORK_DIR/bad_template.log"
  exit 1
fi

{
  echo "test_23_ereaxff_mixed_baoab_geometry_matrix: PASS"
  echo "  geometry_positions = ${positions[*]}"
  echo "  negative_case = unsupported_geometry_template"
} | tee "$SUMMARY"
