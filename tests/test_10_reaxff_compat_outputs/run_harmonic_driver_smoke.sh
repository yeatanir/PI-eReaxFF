#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="$ROOT_DIR/build/test_harmonic_driver_output"
SUMMARY="$ROOT_DIR/build/test_summaries/test_10_harmonic_driver_outputs.summary"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR" "$(dirname "$SUMMARY")"

"$ROOT_DIR/build/harmonic_pimd_main" "$OUT_DIR" > "$OUT_DIR/stdout.log"

for file in pimd_energy.out bead_spread.out thermostat_heat.out centroid.xyz; do
  if [[ ! -s "$OUT_DIR/$file" ]]; then
    echo "test_10_harmonic_driver_outputs: FAIL missing $file" | tee "$SUMMARY"
    exit 1
  fi
done

if ! grep -q "ring-polymer convention: standard_ipi" "$OUT_DIR/stdout.log"; then
  echo "test_10_harmonic_driver_outputs: FAIL missing convention diagnostic" | tee "$SUMMARY"
  exit 1
fi

energy_lines=$(grep -vc '^#' "$OUT_DIR/pimd_energy.out")
if [[ "$energy_lines" -ne 21 ]]; then
  echo "test_10_harmonic_driver_outputs: FAIL expected 21 energy frames, got $energy_lines" | tee "$SUMMARY"
  exit 1
fi

final_heat=$(awk 'NF && $1 !~ /^#/ {v=$2} END {print v}' "$OUT_DIR/thermostat_heat.out")
if [[ -z "$final_heat" || "$final_heat" == "0.00000000000000E+00" ]]; then
  echo "test_10_harmonic_driver_outputs: FAIL thermostat heat did not change" | tee "$SUMMARY"
  exit 1
fi

{
  echo "test_10_harmonic_driver_outputs: PASS"
  echo "  output_dir = $OUT_DIR"
  echo "  energy_frames = $energy_lines"
  echo "  final_thermostat_heat_eV = $final_heat"
} | tee "$SUMMARY"
