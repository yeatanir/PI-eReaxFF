#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_20_ereaxff_mixed_baoab_policy_matrix.summary"

rm -rf "$ROOT_DIR"/build/test_ereaxff_mixed_baoab_policy_*
mkdir -p "$SUMMARY_DIR"

policies=(bead_local_acks2 centroid_acks2)
beads_list=(1 2 4 8)

for policy in "${policies[@]}"; do
  for beads in "${beads_list[@]}"; do
    OUT_DIR="$ROOT_DIR/build/test_ereaxff_mixed_baoab_policy_${policy}_p${beads}"
    mkdir -p "$OUT_DIR"
    "$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$OUT_DIR" "$beads" 4 "$policy" none "$OUT_DIR/restart.dat" > "$OUT_DIR/stdout.log"
    for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out electron_localization.out centroid.xyz fort.7 fort.8 molfra.out xmolout restart.dat; do
      if [[ ! -s "$OUT_DIR/$file" ]]; then
        echo "test_20_ereaxff_mixed_baoab_policy_matrix: FAIL missing $file for policy=$policy P=$beads" | tee "$SUMMARY"
        exit 1
      fi
    done
    python3 - "$OUT_DIR/pimd_energy.out" "$policy" "$beads" <<'PY'
import sys

path, policy, beads = sys.argv[1], sys.argv[2], int(sys.argv[3])
last = None
with open(path, "r", encoding="utf-8") as handle:
    for raw in handle:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        last = line.split()

if last is None:
    raise SystemExit("missing energy line")

charge_diff = float(last[5])
if policy == "centroid_acks2":
    if abs(charge_diff) > 1.0e-12:
        raise SystemExit(f"centroid charge diff nonzero: {charge_diff}")
else:
    if beads == 1:
      if abs(charge_diff) > 1.0e-12:
          raise SystemExit(f"P=1 local charge diff nonzero: {charge_diff}")
    elif charge_diff <= 0.0:
        raise SystemExit(f"local charge diff not positive: {charge_diff}")
PY
  done
done

{
  echo "test_20_ereaxff_mixed_baoab_policy_matrix: PASS"
  echo "  policies = bead_local_acks2 centroid_acks2"
  echo "  p_values = 1 2 4 8"
} | tee "$SUMMARY"
