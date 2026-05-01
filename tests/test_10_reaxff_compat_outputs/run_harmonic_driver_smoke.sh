#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="$ROOT_DIR/build/test_harmonic_driver_output"
CONT_DIR="$ROOT_DIR/build/test_harmonic_driver_continuous"
PART1_DIR="$ROOT_DIR/build/test_harmonic_driver_part1"
PART2_DIR="$ROOT_DIR/build/test_harmonic_driver_part2"
SUMMARY="$ROOT_DIR/build/test_summaries/test_10_harmonic_driver_outputs.summary"

rm -rf "$OUT_DIR" "$CONT_DIR" "$PART1_DIR" "$PART2_DIR"
mkdir -p "$OUT_DIR" "$CONT_DIR" "$PART1_DIR" "$PART2_DIR" "$(dirname "$SUMMARY")"

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

"$ROOT_DIR/build/harmonic_pimd_main" "$CONT_DIR" 2000 none "$CONT_DIR/restart.dat" > "$CONT_DIR/stdout.log"
"$ROOT_DIR/build/harmonic_pimd_main" "$PART1_DIR" 1000 none "$PART1_DIR/restart.dat" > "$PART1_DIR/stdout.log"
"$ROOT_DIR/build/harmonic_pimd_main" "$PART2_DIR" 2000 "$PART1_DIR/restart.dat" "$PART2_DIR/restart.dat" > "$PART2_DIR/stdout.log"

python3 - "$CONT_DIR/restart.dat" "$PART2_DIR/restart.dat" <<'PY'
import math
import sys

def read_restart(path):
    data = {}
    q = []
    v = []
    mode = None
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line == "q_nm":
                mode = "q"
                continue
            if line == "v_nm":
                mode = "v"
                continue
            parts = line.split()
            if mode == "q":
                q.append(float(parts[1]))
            elif mode == "v":
                v.append(float(parts[1]))
            else:
                data[parts[0]] = parts[1]
    data["q"] = q
    data["v"] = v
    return data

continuous = read_restart(sys.argv[1])
split = read_restart(sys.argv[2])

for key in ("step", "beads", "rng_state"):
    if continuous[key] != split[key]:
        raise SystemExit(f"restart mismatch for {key}: {continuous[key]} != {split[key]}")

for key in ("thermostat_heat_eV",):
    if abs(float(continuous[key]) - float(split[key])) > 1e-12:
        raise SystemExit(f"restart mismatch for {key}")

for label in ("q", "v"):
    if len(continuous[label]) != len(split[label]):
        raise SystemExit(f"restart length mismatch for {label}")
    max_diff = max(abs(a - b) for a, b in zip(continuous[label], split[label]))
    if max_diff > 1e-12:
        raise SystemExit(f"restart {label} max diff {max_diff}")
PY

{
  echo "test_10_harmonic_driver_outputs: PASS"
  echo "  output_dir = $OUT_DIR"
  echo "  energy_frames = $energy_lines"
  echo "  final_thermostat_heat_eV = $final_heat"
  echo "  restart_equivalence = continuous_vs_split"
} | tee "$SUMMARY"
