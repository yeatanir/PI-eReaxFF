#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
V1_ROOT="${EREAXFF_V1_ROOT:-/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1}"
V2_ROOT="${EREAXFF_V2_ROOT:-/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2}"
WORK_DIR="$ROOT_DIR/build/test_v2_oracle_smoke"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SMOKE_BIN="$V2_ROOT/phase4_single_point_oracle"
REPEAT_BIN="$V2_ROOT/single_point_repeatability"
MAKE_GEO="$V2_ROOT/tests/regression/test_phase4_fixed_graphene/make_geo.py"
TIMEOUT_RUNNER="$V2_ROOT/tests/regression/test_force_repeatability/run_with_timeout.py"

mkdir -p "$WORK_DIR" "$SUMMARY_DIR"
SUMMARY="$SUMMARY_DIR/test_08_reaxff_oracle_smoke.summary"

fail() {
  {
    echo "test_08_reaxff_oracle_smoke: FAIL"
    echo "$1"
  } > "$SUMMARY"
  cat "$SUMMARY"
  exit 1
}

[[ -x "$SMOKE_BIN" ]] || fail "Missing executable: $SMOKE_BIN"
[[ -x "$REPEAT_BIN" ]] || fail "Missing executable: $REPEAT_BIN"
[[ -f "$V1_ROOT/tests/T4_60el/control" ]] || fail "Missing V1 control fixture."
[[ -f "$V1_ROOT/tests/T4_60el/ffield" ]] || fail "Missing V1 ffield fixture."
[[ -f "$V1_ROOT/tests/reference/iopt" ]] || fail "Missing V1 iopt fixture."
[[ -f "$MAKE_GEO" ]] || fail "Missing V2 make_geo.py fixture generator."
[[ -x "$TIMEOUT_RUNNER" ]] || fail "Missing executable timeout runner."

cp "$V1_ROOT/tests/T4_60el/control" "$WORK_DIR/control"
cp "$V1_ROOT/tests/T4_60el/ffield" "$WORK_DIR/ffield"
cp "$V1_ROOT/tests/reference/iopt" "$WORK_DIR/fort.20"
printf "0\n" > "$WORK_DIR/fort.35"

(
  cd "$WORK_DIR"
  python3 "$MAKE_GEO"
  python3 "$TIMEOUT_RUNNER" 90 run_single_point.log "$SMOKE_BIN"
  python3 "$TIMEOUT_RUNNER" 90 run_repeatability.log "$REPEAT_BIN" 3
)

python3 - "$WORK_DIR" "$SUMMARY" <<'PY'
import math
import pathlib
import sys

work = pathlib.Path(sys.argv[1])
summary = pathlib.Path(sys.argv[2])
single = work / "phase4_single_point_oracle.out"
repeat = work / "single_point_repeatability.out"

if not single.exists():
    raise SystemExit("missing phase4_single_point_oracle.out")
if not repeat.exists():
    raise SystemExit("missing single_point_repeatability.out")

single_lines = [
    line.split()
    for line in single.read_text(encoding="utf-8", errors="replace").splitlines()
    if line.strip()
]
if len(single_lines) < 2:
    raise SystemExit("single-point output is missing header or data line")

header = single_lines[0]
if header and header[0] == "#":
    header = header[1:]
header = [item.lstrip("#") for item in header]
values = single_lines[1]
if len(header) != len(values):
    raise SystemExit("single-point header/data field count mismatch")

fields = dict(zip(header, values))
required = {"status", "natom", "energy_kcal", "max_force", "charge_sum"}
missing = sorted(required - fields.keys())
if missing:
    raise SystemExit(f"missing single-point fields: {', '.join(missing)}")

status = int(float(fields["status"]))
natom = int(float(fields["natom"]))
energy = float(fields["energy_kcal"])
max_force = float(fields["max_force"])
charge_sum = float(fields["charge_sum"])

checks = [
    (status == 0, f"oracle status is {status}, expected 0"),
    (natom > 0, f"natom is {natom}, expected > 0"),
    (math.isfinite(energy), "energy is not finite"),
    (math.isfinite(max_force), "max_force is not finite"),
    (math.isfinite(charge_sum), "charge_sum is not finite"),
]
for ok, message in checks:
    if not ok:
        raise SystemExit(message)

repeat_lines = [
    line
    for line in repeat.read_text(encoding="utf-8", errors="replace").splitlines()
    if line.strip() and not line.startswith("#")
]
if len(repeat_lines) < 3:
    raise SystemExit("repeatability output has fewer than 3 evaluations")

max_force_diff = 0.0
max_charge_diff = 0.0
for line in repeat_lines:
    parts = line.split()
    if len(parts) < 9:
        raise SystemExit(f"malformed repeatability line: {line}")
    rep_status = int(float(parts[8]))
    if rep_status != 0:
        raise SystemExit(f"repeatability status is {rep_status}, expected 0")
    values = [float(item) for item in parts[1:8]]
    if not all(math.isfinite(value) for value in values):
        raise SystemExit(f"non-finite repeatability value: {line}")
    max_force_diff = max(max_force_diff, abs(float(parts[2])))
    max_charge_diff = max(max_charge_diff, abs(float(parts[3])))

summary.write_text(
    "\n".join(
        [
            "test_08_reaxff_oracle_smoke: PASS",
            f"natom = {natom}",
            f"energy_kcal = {energy:.15e}",
            f"max_force = {max_force:.15e}",
            f"charge_sum = {charge_sum:.15e}",
            f"repeatability_evaluations = {len(repeat_lines)}",
            f"max_repeat_force_diff = {max_force_diff:.15e}",
            f"max_repeat_charge_diff = {max_charge_diff:.15e}",
        ]
    )
    + "\n",
    encoding="utf-8",
)
PY

cat "$SUMMARY"
