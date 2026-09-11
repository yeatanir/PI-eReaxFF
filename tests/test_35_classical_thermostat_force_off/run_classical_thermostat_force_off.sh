#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_35_classical_thermostat_force_off.summary"
WORK_DIR="$ROOT_DIR/build/test_35_classical_thermostat_force_off"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

run_case() {
  local thermostat="$1"
  local cfg="$WORK_DIR/config_${thermostat}.yaml"
  local out="$WORK_DIR/out_${thermostat}"
  mkdir -p "$out"

  python3 - "$ROOT_DIR/config/example_mixed_classical_pimd.yaml" "$cfg" "$out" "$thermostat" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
target = pathlib.Path(sys.argv[2])
out_dir = sys.argv[3]
thermostat = sys.argv[4]

source = source.replace("  dt_fs: 0.1", "  dt_fs: 0.05")
source = source.replace("  steps: 1000", "  steps: 100")
source = source.replace("  output_stride: 10", "  output_stride: 10")
source = source.replace("  restart_stride: 100", "  restart_stride: 100")
source = source.replace("  output_dir: \"\"", f"  output_dir: \"{out_dir}\"")
source = source.replace("  type: ereaxff_oracle", "  type: dummy")
source = source.replace("  default_beads: 16", "  default_beads: 4")
source = source.replace("    beads: 16", "    beads: 4")

if thermostat == "berendsen":
    source = source.replace(
        "    thermostat: langevin\n    temperature_K: 300.0\n    gamma_fs_inv: 0.001",
        "    thermostat: berendsen\n    temperature_K: 300.0\n    tau_fs: 10.0",
        1,
    )
elif thermostat == "langevin":
    source = source.replace(
        "    gamma_fs_inv: 0.001",
        "    gamma_fs_inv: 0.01",
        1,
    )
else:
    raise SystemExit(f"unsupported thermostat {thermostat}")

target.write_text(source, encoding="utf-8")
PY

  "$ROOT_DIR/build/ereaxff_mixed_baoab_main" "$cfg" "$out" > "$out/stdout.log"

  for file in carbon_temperature.out thermostat_heat.out pimd_energy.out restart.dat; do
    if [[ ! -s "$out/$file" ]]; then
      echo "test_35_classical_thermostat_force_off: FAIL missing $file for $thermostat" | tee "$SUMMARY"
      exit 1
    fi
  done

  grep -q "force provider: dummy zero-force" "$out/stdout.log"
  grep -q "thermostat=$thermostat" "$out/stdout.log"
}

run_case berendsen
run_case langevin

python3 - "$WORK_DIR/out_berendsen" "$WORK_DIR/out_langevin" <<'PY'
import math
import pathlib
import sys

def rows(path):
    data = []
    for raw in pathlib.Path(path, "carbon_temperature.out").read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        step, temp, ke = raw.split()
        data.append((int(step), float(temp), float(ke)))
    return data

for label, path in (("berendsen", sys.argv[1]), ("langevin", sys.argv[2])):
    data = rows(path)
    if len(data) < 2:
        raise SystemExit(f"{label}: insufficient carbon temperature rows")
    temps = [row[1] for row in data]
    if not all(math.isfinite(value) for value in temps):
        raise SystemExit(f"{label}: non-finite carbon temperature")
    if max(temps) > 1200.0:
        raise SystemExit(f"{label}: force-off carbon thermostat runaway, max T={max(temps):.3f} K")
    if temps[-1] <= 0.0:
        raise SystemExit(f"{label}: non-positive final carbon temperature")

berendsen = rows(sys.argv[1])
initial = berendsen[0][1]
final = berendsen[-1][1]
if abs(final - 300.0) > abs(initial - 300.0) + 1.0e-8:
    raise SystemExit(
        f"berendsen: temperature did not relax toward target, initial={initial:.6f} final={final:.6f}"
    )

print(f"berendsen_initial_K={initial:.6f}")
print(f"berendsen_final_K={final:.6f}")
print(f"langevin_final_K={rows(sys.argv[2])[-1][1]:.6f}")
PY

{
  echo "test_35_classical_thermostat_force_off: PASS"
  echo "  cases = berendsen langevin"
  echo "  force_provider = dummy"
  echo "  dt_fs = 0.05"
  echo "  steps = 100"
} | tee "$SUMMARY"
