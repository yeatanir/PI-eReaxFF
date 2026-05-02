#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_21_ereaxff_mixed_baoab_config_smoke.summary"
WORK_DIR="$ROOT_DIR/build/test_21_ereaxff_mixed_baoab_config_smoke"
CFG="$WORK_DIR/config.yaml"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$CFG" "$WORK_DIR" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
work_dir = sys.argv[3]
source = source.replace("  steps: 1000", "  steps: 4")
source = source.replace("  output_dir: \"\"", f"  output_dir: \"{work_dir}\"")
source = source.replace("  default_beads: 16", "  default_beads: 4")
source = source.replace("    beads: 16", "    beads: 4")
target.write_text(source, encoding="utf-8")
PY

"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$CFG" "$WORK_DIR" > "$WORK_DIR/stdout.log"

for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out electron_localization.out centroid.xyz fort.7 fort.8 molfra.out xmolout restart.dat; do
  if [[ ! -s "$WORK_DIR/$file" ]]; then
    echo "test_21_ereaxff_mixed_baoab_config_smoke: FAIL missing $file" | tee "$SUMMARY"
    exit 1
  fi
done

grep -q "input file: $CFG" "$WORK_DIR/stdout.log"
grep -q "config-driven legacy V1/V2 fixtures" "$WORK_DIR/stdout.log"
grep -q "P: 4" "$WORK_DIR/stdout.log"

{
  echo "test_21_ereaxff_mixed_baoab_config_smoke: PASS"
  echo "  config = config/example_mixed_classical_pimd.yaml"
  echo "  beads = 4"
  echo "  steps = 4"
} | tee "$SUMMARY"
