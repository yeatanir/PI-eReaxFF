#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_22_ereaxff_mixed_baoab_config_matrix.summary"
WORK_DIR="$ROOT_DIR/build/test_22_ereaxff_mixed_baoab_config_matrix"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

policies=(bead_local_acks2 centroid_acks2)
beads_list=(1 2 4 8)

for policy in "${policies[@]}"; do
  for beads in "${beads_list[@]}"; do
    cfg="$WORK_DIR/config_${policy}_p${beads}.yaml"
    out="$WORK_DIR/out_${policy}_p${beads}"
    mkdir -p "$out"
    python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$cfg" "$out" "$beads" "$policy" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
out_dir = sys.argv[3]
beads = sys.argv[4]
policy = sys.argv[5]
source = source.replace("  steps: 1000", "  steps: 4")
source = source.replace("  output_stride: 10", "  output_stride: 2")
source = source.replace("  restart_stride: 100", "  restart_stride: 2")
source = source.replace("  output_dir: \"\"", f"  output_dir: \"{out_dir}\"")
source = source.replace("  default_beads: 16", f"  default_beads: {beads}")
source = source.replace("  mode: bead_local_acks2", f"  mode: {policy}")
source = source.replace("    beads: 16", f"    beads: {beads}")
target.write_text(source, encoding="utf-8")
PY
    "$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$cfg" "$out" > "$out/stdout.log"
    for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out electron_localization.out centroid.xyz fort.7 fort.8 molfra.out xmolout restart.dat; do
      if [[ ! -s "$out/$file" ]]; then
        echo "test_22_ereaxff_mixed_baoab_config_matrix: FAIL missing $file for policy=$policy P=$beads" | tee "$SUMMARY"
        exit 1
      fi
    done
    grep -q "input file: $cfg" "$out/stdout.log"
    grep -q "P: $beads" "$out/stdout.log"
  done
done

bad_cfg="$WORK_DIR/bad_selector.yaml"
python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$bad_cfg" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
source = source.replace("element == C", "element = C", 1)
target.write_text(source, encoding="utf-8")
PY

set +e
"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$bad_cfg" "$WORK_DIR/bad_out" > "$WORK_DIR/bad_selector.log" 2>&1
bad_status=$?
set -e
if [[ "$bad_status" -eq 0 ]]; then
  echo "test_22_ereaxff_mixed_baoab_config_matrix: FAIL invalid selector unexpectedly succeeded" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "unsupported selector" "$WORK_DIR/bad_selector.log"; then
  echo "test_22_ereaxff_mixed_baoab_config_matrix: FAIL expected selector error not found" | tee "$SUMMARY"
  cat "$WORK_DIR/bad_selector.log"
  exit 1
fi

overlap_cfg="$WORK_DIR/overlap.yaml"
python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$overlap_cfg" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
source = source.replace("  - name: other_electrons\n    selector: \"element == El and id != 61\"\n    dynamics: frozen\n    thermostat: none\n", "")
source = source.replace("  - name: selected_electrons\n    selector: \"element == El and id in [61]\"\n    dynamics: pimd\n    beads: 16\n    thermostat: pile_l\n    temperature_K: 300.0\n    centroid_gamma_fs_inv: 0.001\n", "  - name: selected_electrons\n    selector: \"all\"\n    dynamics: pimd\n    beads: 4\n    thermostat: pile_l\n    temperature_K: 300.0\n    centroid_gamma_fs_inv: 0.001\n")
target.write_text(source, encoding="utf-8")
PY

set +e
"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$overlap_cfg" "$WORK_DIR/overlap_out" > "$WORK_DIR/overlap.log" 2>&1
overlap_status=$?
set -e
if [[ "$overlap_status" -eq 0 ]]; then
  echo "test_22_ereaxff_mixed_baoab_config_matrix: FAIL overlapping selector unexpectedly succeeded" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "particle matches multiple groups" "$WORK_DIR/overlap.log"; then
  echo "test_22_ereaxff_mixed_baoab_config_matrix: FAIL expected overlap error not found" | tee "$SUMMARY"
  cat "$WORK_DIR/overlap.log"
  exit 1
fi

{
  echo "test_22_ereaxff_mixed_baoab_config_matrix: PASS"
  echo "  policies = bead_local_acks2 centroid_acks2"
  echo "  p_values = 1 2 4 8"
  echo "  negative_cases = invalid_selector overlapping_selector"
} | tee "$SUMMARY"
