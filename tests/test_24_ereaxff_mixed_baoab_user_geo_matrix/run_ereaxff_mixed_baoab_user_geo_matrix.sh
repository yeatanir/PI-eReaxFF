#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_24_ereaxff_mixed_baoab_user_geo_matrix.summary"
WORK_DIR="$ROOT_DIR/build/test_24_ereaxff_mixed_baoab_user_geo_matrix"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

cfg="$WORK_DIR/user_geo_config.yaml"
out="$WORK_DIR/out"
mkdir -p "$out"

python3 - "$ROOT_DIR/config/example_user_geo_structure.yaml" "$cfg" "$out" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
out_dir = sys.argv[3]
source = source.replace("  steps: 100", "  steps: 4")
source = source.replace("  output_stride: 10", "  output_stride: 2")
source = source.replace("  restart_stride: 20", "  restart_stride: 2")
source = source.replace("  output_dir: \"\"", f"  output_dir: \"{out_dir}\"")
target.write_text(source, encoding="utf-8")
PY

"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$cfg" "$out" > "$out/stdout.log"

for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out electron_localization.out centroid.xyz fort.7 fort.8 molfra.out xmolout restart.dat; do
  if [[ ! -s "$out/$file" ]]; then
    echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL missing $file" | tee "$SUMMARY"
    exit 1
  fi
done

if ! grep -q "connection table" "$out/fort.7"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL fort.7 missing connection-table header" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "connection table" "$out/fort.8"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL fort.8 missing connection-table header" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "molecule list per frame" "$out/molfra.out"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL molfra.out missing molecule-list header" | tee "$SUMMARY"
  exit 1
fi

if ! grep -q "geometry template: user_geo_structure" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL missing geometry template line" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "geometry source:" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL missing geometry source line" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "group thermostat: graphene_carbons" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL missing graphene thermostat diagnostic" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "group thermostat: excess_electron" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL missing electron thermostat diagnostic" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "thermostat=langevin" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL graphene thermostat mode mismatch" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "thermostat=pile_l" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL electron thermostat mode mismatch" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "classical_count = 66" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL unexpected classical count" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "pimd_count = 1" "$out/stdout.log"; then
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: FAIL unexpected pimd count" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi

{
  echo "test_24_ereaxff_mixed_baoab_user_geo_matrix: PASS"
  echo "  geometry_template = user_geo_structure"
  echo "  geometry_source = $ROOT_DIR/config/example_user_geo_structure.yaml"
  echo "  particle_counts = classical_66 pimd_1"
} | tee "$SUMMARY"
