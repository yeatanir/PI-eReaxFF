#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_19_ereaxff_mixed_baoab_restart.summary"

CONT_DIR="$ROOT_DIR/build/test_ereaxff_mixed_baoab_continuous"
PART1_DIR="$ROOT_DIR/build/test_ereaxff_mixed_baoab_part1"
PART2_DIR="$ROOT_DIR/build/test_ereaxff_mixed_baoab_part2"

rm -rf "$CONT_DIR" "$PART1_DIR" "$PART2_DIR"
mkdir -p "$CONT_DIR" "$PART1_DIR" "$PART2_DIR" "$SUMMARY_DIR"

"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$CONT_DIR" 4 4 bead_local_acks2 none "$CONT_DIR/restart.dat" > "$CONT_DIR/stdout.log"
"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$PART1_DIR" 4 2 bead_local_acks2 none "$PART1_DIR/restart.dat" > "$PART1_DIR/stdout.log"
"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$PART2_DIR" 4 4 bead_local_acks2 "$PART1_DIR/restart.dat" "$PART2_DIR/restart.dat" > "$PART2_DIR/stdout.log"

for file in charge_by_bead.out charge_centroid.out electron_localization.out fort.7 fort.8 molfra.out xmolout restart.dat; do
  if ! cmp -s "$CONT_DIR/$file" "$PART2_DIR/$file"; then
    echo "test_19_ereaxff_mixed_baoab_restart: FAIL $file mismatch after restart" | tee "$SUMMARY"
    exit 1
  fi
done

python3 - "$CONT_DIR/pimd_energy.out" "$PART2_DIR/pimd_energy.out" "$CONT_DIR/bead_spread.out" "$PART2_DIR/bead_spread.out" <<'PY'
import sys

def last_data_line(path):
    last = None
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            last = line
    return last

for left, right in ((sys.argv[1], sys.argv[2]), (sys.argv[3], sys.argv[4])):
    if last_data_line(left) != last_data_line(right):
        raise SystemExit(f"last line mismatch: {left} != {right}")
PY

{
  echo "test_19_ereaxff_mixed_baoab_restart: PASS"
  echo "  split_run = continuous_vs_restart"
} | tee "$SUMMARY"
