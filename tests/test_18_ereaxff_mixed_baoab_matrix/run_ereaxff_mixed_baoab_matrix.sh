#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_18_ereaxff_mixed_baoab_matrix.summary"

rm -rf "$ROOT_DIR"/build/test_ereaxff_mixed_baoab_*
mkdir -p "$SUMMARY_DIR"

for beads in 1 2 4 8; do
  OUT_DIR="$ROOT_DIR/build/test_ereaxff_mixed_baoab_p${beads}"
  mkdir -p "$OUT_DIR"
  "$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$OUT_DIR" "$beads" 4 bead_local_acks2 > "$OUT_DIR/stdout.log"
  for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out electron_localization.out centroid.xyz fort.7 fort.8 molfra.out xmolout restart.dat; do
    if [[ ! -s "$OUT_DIR/$file" ]]; then
      echo "test_18_ereaxff_mixed_baoab_matrix: FAIL missing $file for P=$beads" | tee "$SUMMARY"
      exit 1
    fi
  done
  if ! grep -q "integrator mode: mixed_baoab" "$OUT_DIR/stdout.log"; then
    echo "test_18_ereaxff_mixed_baoab_matrix: FAIL missing integrator diagnostic for P=$beads" | tee "$SUMMARY"
    exit 1
  fi
  energy_lines=$(grep -vc '^#' "$OUT_DIR/pimd_energy.out")
  if [[ "$energy_lines" -ne 5 ]]; then
    echo "test_18_ereaxff_mixed_baoab_matrix: FAIL expected 5 energy frames for P=$beads, got $energy_lines" | tee "$SUMMARY"
    exit 1
  fi
done

{
  echo "test_18_ereaxff_mixed_baoab_matrix: PASS"
  echo "  p_values = 1 2 4 8"
} | tee "$SUMMARY"
