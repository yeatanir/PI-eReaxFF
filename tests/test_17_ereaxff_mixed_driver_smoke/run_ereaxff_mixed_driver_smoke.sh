#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT1="$ROOT_DIR/build/test_ereaxff_mixed_driver_out1"
OUT2="$ROOT_DIR/build/test_ereaxff_mixed_driver_out2"
SUMMARY="$ROOT_DIR/build/test_summaries/test_17_ereaxff_mixed_driver_smoke.summary"

rm -rf "$OUT1" "$OUT2"
mkdir -p "$OUT1" "$OUT2" "$(dirname "$SUMMARY")"

"$ROOT_DIR/build/ereaxff_pimd_main" "$OUT1" 4 bead_local_acks2 none "$OUT1/restart.dat" > "$OUT1/stdout.log"
"$ROOT_DIR/build/ereaxff_pimd_main" "$OUT2" 4 bead_local_acks2 "$OUT1/restart.dat" "$OUT2/restart.dat" > "$OUT2/stdout.log"

for file in pimd_energy.out charge_by_bead.out charge_centroid.out electron_localization.out restart.dat; do
  if [[ ! -s "$OUT1/$file" ]]; then
    echo "test_17_ereaxff_mixed_driver_smoke: FAIL missing $file in out1" | tee "$SUMMARY"
    exit 1
  fi
  if [[ ! -s "$OUT2/$file" ]]; then
    echo "test_17_ereaxff_mixed_driver_smoke: FAIL missing $file in out2" | tee "$SUMMARY"
    exit 1
  fi
done

if ! grep -q "force provider: eReaxFF adapter" "$OUT1/stdout.log"; then
  echo "test_17_ereaxff_mixed_driver_smoke: FAIL missing eReaxFF diagnostic" | tee "$SUMMARY"
  exit 1
fi

if ! cmp -s "$OUT1/pimd_energy.out" "$OUT2/pimd_energy.out"; then
  echo "test_17_ereaxff_mixed_driver_smoke: FAIL energy output changed across restart" | tee "$SUMMARY"
  exit 1
fi

if ! cmp -s "$OUT1/charge_by_bead.out" "$OUT2/charge_by_bead.out"; then
  echo "test_17_ereaxff_mixed_driver_smoke: FAIL bead charge output changed across restart" | tee "$SUMMARY"
  exit 1
fi

if ! cmp -s "$OUT1/charge_centroid.out" "$OUT2/charge_centroid.out"; then
  echo "test_17_ereaxff_mixed_driver_smoke: FAIL centroid charge output changed across restart" | tee "$SUMMARY"
  exit 1
fi

if ! cmp -s "$OUT1/electron_localization.out" "$OUT2/electron_localization.out"; then
  echo "test_17_ereaxff_mixed_driver_smoke: FAIL localization output changed across restart" | tee "$SUMMARY"
  exit 1
fi

python3 - "$OUT1/restart.dat" "$OUT2/restart.dat" <<'PY'
import sys

def parse(path):
    data = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line in {"positions_A", "charges_by_bead"}:
                break
            key, value = line.split(maxsplit=1)
            data[key] = value
    return data

first = parse(sys.argv[1])
second = parse(sys.argv[2])
for key in ("step", "ndim", "nparticles", "beads", "charge_policy"):
    if first.get(key) != second.get(key):
        raise SystemExit(f"restart mismatch for {key}: {first.get(key)} != {second.get(key)}")
PY

{
  echo "test_17_ereaxff_mixed_driver_smoke: PASS"
  echo "  output_dir_1 = $OUT1"
  echo "  output_dir_2 = $OUT2"
} | tee "$SUMMARY"
