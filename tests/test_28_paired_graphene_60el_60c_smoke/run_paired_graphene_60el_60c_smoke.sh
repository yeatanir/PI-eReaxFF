#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_28_paired_graphene_60el_60c_smoke.summary"
WORK_DIR="$ROOT_DIR/build/test_28_paired_graphene_60el_60c_smoke"
CFG="$ROOT_DIR/config/example_paired_graphene_60el_60c.yaml"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

"$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$CFG" "$WORK_DIR/out" > "$WORK_DIR/stdout.log"

for file in pimd_energy.out bead_spread.out charge_by_bead.out charge_centroid.out centroid.xyz restart.dat; do
  if [[ ! -s "$WORK_DIR/out/$file" ]]; then
    echo "test_28_paired_graphene_60el_60c_smoke: FAIL missing $file" | tee "$SUMMARY"
    exit 1
  fi
done

grep -q "geometry template: user_geo_structure" "$WORK_DIR/stdout.log"
grep -q "classical_count = 60" "$WORK_DIR/stdout.log"
grep -q "pimd_count = 60" "$WORK_DIR/stdout.log"
grep -q "nparticles 120" "$WORK_DIR/out/restart.dat"
grep -q "nquantum 60" "$WORK_DIR/out/restart.dat"
grep -q "beads 2" "$WORK_DIR/out/restart.dat"

python3 - "$WORK_DIR/out/restart.dat" "$WORK_DIR/out/pimd_energy.out" <<'PY'
from __future__ import annotations

import math
import pathlib
import sys

restart = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
energy = pathlib.Path(sys.argv[2]).read_text(encoding="utf-8").splitlines()

q_rows = []
in_q = False
for line in restart:
    if line.strip() == "quantum_indices":
        in_q = True
        continue
    if in_q and line.strip() == "positions_A":
        break
    if in_q and line.strip():
        q_rows.append(line.split())

if len(q_rows) != 60:
    raise SystemExit(f"FAIL: expected 60 quantum rows, got {len(q_rows)}")
if any(int(row[2]) != 2 for row in q_rows):
    raise SystemExit("FAIL: expected all paired electrons to carry P=2")

data_rows = [line.split() for line in energy if line.strip() and not line.startswith("#")]
if not data_rows:
    raise SystemExit("FAIL: missing energy rows")
for row in data_rows:
    vals = [float(x) for x in row[1:]]
    if not all(math.isfinite(x) for x in vals):
        raise SystemExit("FAIL: non-finite energy diagnostic")
PY

{
  echo "test_28_paired_graphene_60el_60c_smoke: PASS"
  echo "  geometry = paired_graphene_60el.geo"
  echo "  particles = C_60 El_60"
  echo "  quantum_electrons = 60"
  echo "  beads = 2"
} | tee "$SUMMARY"
