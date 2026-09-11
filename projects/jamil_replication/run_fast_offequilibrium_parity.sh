#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/projects/jamil_replication"
BASE_CFG="$SCRIPT_DIR/config_jamil_112c_113el_p4_pool_verify.yaml"
SOURCE_RUN="${JAMIL_FAST_PARITY_SOURCE_RUN:-$ROOT_DIR/build/jamil_replication/fast_stability_long_dt005/fast_p16_steps100_dt0p05/out}"
PARITY_ROOT="${JAMIL_FAST_PARITY_ROOT:-$ROOT_DIR/build/jamil_replication/fast_offequilibrium_parity_$(date +%Y%m%d_%H%M%S)}"
STEPS="${JAMIL_FAST_PARITY_STEPS:-0 10 50 100}"
BEADS="${JAMIL_FAST_PARITY_BEADS:-1 8 16}"
SUMMARY="$SCRIPT_DIR/fast_offequilibrium_parity_summary.out"
DATASET="$PARITY_ROOT/offequilibrium_oracle_cases.jsonl"
REPORT="$PARITY_ROOT/offequilibrium_force_report.json"
MAX_FORCE_ERROR="${JAMIL_FAST_PARITY_MAX_FORCE_ERROR:-}"

if [[ ! -d "$SOURCE_RUN" ]]; then
  echo "run_fast_offequilibrium_parity: missing source run: $SOURCE_RUN" >&2
  echo "run JAMIL_FAST_STABILITY_RUN_LONG=1 projects/jamil_replication/run_fast_stability_check.sh first, or set JAMIL_FAST_PARITY_SOURCE_RUN" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/build/oracle_single_point_force_export" ]]; then
  echo "run_fast_offequilibrium_parity: missing build/oracle_single_point_force_export; run make oracle_single_point_force_export" >&2
  exit 1
fi
if [[ ! -x "$ROOT_DIR/.venv-fast/bin/python" ]]; then
  echo "run_fast_offequilibrium_parity: missing .venv-fast/bin/python" >&2
  exit 77
fi

mkdir -p "$PARITY_ROOT"
rm -f "$SUMMARY"

python3 - "$ROOT_DIR" "$BASE_CFG" "$SOURCE_RUN" "$PARITY_ROOT/adapter_work" "$STEPS" "$BEADS" <<'PY'
from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import sys
import time

root = pathlib.Path(sys.argv[1])
base_cfg = pathlib.Path(sys.argv[2])
source_run = pathlib.Path(sys.argv[3])
adapter_work = pathlib.Path(sys.argv[4])
selected_steps = [int(value) for value in sys.argv[5].split()]
selected_beads = [int(value) for value in sys.argv[6].split()]


def config_value(text: str, key: str) -> str:
    match = re.search(rf"^  {re.escape(key)}:\s*(.*)$", text, flags=re.M)
    if not match:
        raise ValueError(f"missing force_provider.{key} in {base_cfg}")
    return match.group(1).strip().strip('"')


def resolve_path(v1_root: pathlib.Path, value: str) -> pathlib.Path:
    path = pathlib.Path(value)
    if path.is_absolute():
        return path
    return (v1_root / path).resolve()


def read_xyz_frames(path: pathlib.Path) -> dict[tuple[int, int], list[tuple[str, float, float, float]]]:
    frames: dict[tuple[int, int], list[tuple[str, float, float, float]]] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        natoms = int(lines[cursor].strip())
        header = lines[cursor + 1].split()
        meta = {}
        for token in header:
            if "=" in token:
                key, value = token.split("=", 1)
                meta[key] = int(value)
        if "step" not in meta or "bead" not in meta:
            raise ValueError(f"{path}: XYZ frame missing step/bead metadata")
        atoms = []
        for raw in lines[cursor + 2 : cursor + 2 + natoms]:
            fields = raw.split()
            atoms.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
        frames[(meta["step"], meta["bead"])] = atoms
        cursor += natoms + 2
    return frames


def write_geo(template: pathlib.Path, target: pathlib.Path, atoms: list[tuple[str, float, float, float]]) -> None:
    atom_cursor = 0
    out_lines = []
    for raw in template.read_text(encoding="utf-8").splitlines():
        if raw.startswith("HETATM"):
            if atom_cursor >= len(atoms):
                raise ValueError("too many HETATM rows in template")
            symbol, x, y, z = atoms[atom_cursor]
            prefix = raw[:30]
            suffix = raw[60:] if len(raw) > 60 else ""
            out_lines.append(f"{prefix}{x:10.5f}{y:10.5f}{z:10.5f}{suffix}")
            atom_cursor += 1
        else:
            out_lines.append(raw)
    if atom_cursor != len(atoms):
        raise ValueError("template atom count does not match XYZ frame")
    target.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def copy_control_with_axes(source: pathlib.Path, target: pathlib.Path, crystx: list[float]) -> None:
    axis = {"axis1": crystx[0], "axis2": crystx[1], "axis3": crystx[2]}
    lines = []
    seen = set()
    for raw in source.read_text(encoding="utf-8").splitlines():
        fields = raw.split()
        if len(fields) >= 2 and fields[1] in axis:
            key = fields[1]
            lines.append(f"{axis[key]:12.5f} {key}")
            seen.add(key)
        else:
            lines.append(raw)
    for key, value in axis.items():
        if key not in seen:
            lines.append(f"{value:12.5f} {key}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_crystx(template: pathlib.Path) -> list[float]:
    for raw in template.read_text(encoding="utf-8").splitlines():
        if raw.startswith("CRYSTX"):
            return [float(value) for value in raw.split()[1:7]]
    raise ValueError(f"{template}: missing CRYSTX")


cfg_text = base_cfg.read_text(encoding="utf-8")
v1_root = pathlib.Path(config_value(cfg_text, "v1_path"))
geo_template = pathlib.Path(re.search(r'^  geo_path:\s*"(.*)"$', cfg_text, flags=re.M).group(1))
ffield = resolve_path(v1_root, config_value(cfg_text, "ffield"))
control = resolve_path(v1_root, config_value(cfg_text, "control"))
iopt = resolve_path(v1_root, config_value(cfg_text, "iopt"))
crystx = parse_crystx(geo_template)
frames = read_xyz_frames(source_run / "trajectory_beads.xyz")

if adapter_work.exists():
    shutil.rmtree(adapter_work)
adapter_work.mkdir(parents=True)
timing_rows = [
    "# kind step batch job bead worker dispatch_mode stage_kind stage_label worker_state",
    "# n_jobs start_s end_s elapsed_s status workdir",
    "# generated from fast trajectory frames for off-equilibrium oracle parity",
]
call_slot = 0
for step in selected_steps:
    for bead in selected_beads:
        key = (step, bead)
        if key not in frames:
            raise ValueError(f"missing trajectory_beads.xyz frame step={step} bead={bead}")
        call_slot += 1
        call_dir = adapter_work / f"adapter_call_{call_slot:04d}"
        call_dir.mkdir()
        write_geo(geo_template, call_dir / "geo", frames[key])
        shutil.copy2(ffield, call_dir / "ffield")
        shutil.copy2(iopt, call_dir / "iopt")
        shutil.copy2(iopt, call_dir / "fort.20")
        (call_dir / "fort.35").write_text("0\n", encoding="utf-8")
        copy_control_with_axes(control, call_dir / "control", crystx)
        with (call_dir / "charge_memory.in").open("w", encoding="utf-8") as fh:
            for atom_index in range(1, len(frames[key]) + 1):
                fh.write(f"{atom_index} 0.0 0.0\n")
        start = time.time()
        with (call_dir / "run_single_point.log").open("w", encoding="utf-8") as log:
            proc = subprocess.run([str(root / "build/oracle_single_point_force_export")], cwd=call_dir, stdout=log, stderr=subprocess.STDOUT)
        elapsed = time.time() - start
        if proc.returncode != 0:
            raise SystemExit(f"oracle exporter failed for step={step} bead={bead}; see {call_dir / 'run_single_point.log'}")
        timing_rows.append(
            "CALL "
            f"{step} 1 {call_slot} {bead} 1 serial bead_local fast_frame_s{step}_b{bead} external "
            f"0.0 {elapsed:.6f} {elapsed:.6f} 0 {call_dir}"
        )

(adapter_work / "oracle_timing.log").write_text("\n".join(timing_rows) + "\n", encoding="utf-8")
print(f"wrote {call_slot} oracle call directories under {adapter_work}")
PY

python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py \
  --adapter-work "$PARITY_ROOT/adapter_work" \
  --output "$DATASET"
python3 projects/jamil_replication/fast_backend/build_oracle_dataset.py --check "$DATASET"
.venv-fast/bin/python projects/jamil_replication/fast_backend/compare_force_terms.py \
  --dataset "$DATASET" \
  --output "$REPORT"

python3 - "$DATASET" "$REPORT" "$SUMMARY" "$SOURCE_RUN" "$PARITY_ROOT" "$MAX_FORCE_ERROR" <<'PY'
from __future__ import annotations

import json
import pathlib
import sys

dataset = pathlib.Path(sys.argv[1])
report = pathlib.Path(sys.argv[2])
summary = pathlib.Path(sys.argv[3])
source_run = pathlib.Path(sys.argv[4])
parity_root = pathlib.Path(sys.argv[5])
max_force_error_arg = sys.argv[6]

cases = [json.loads(line) for line in dataset.read_text(encoding="utf-8").splitlines() if line.strip()]
data = json.loads(report.read_text(encoding="utf-8"))
rows = []
for worst in data["worst_cases"]:
    rows.append(
        {
            "stage": worst["stage"],
            "bead": worst["bead"],
            "atom": worst["atom_index_0based"] + 1,
            "component": worst["component"],
            "abs_error": worst["abs_error_eV_per_A"],
            "signed_error": worst["signed_error_eV_per_A"],
            "jax_force": worst["jax_force_eV_per_A"],
            "oracle_force": worst["oracle_force_eV_per_A"],
        }
    )

with summary.open("w", encoding="utf-8") as fh:
    fh.write("# Fast Off-Equilibrium Oracle Parity Summary\n\n")
    fh.write(f"- source_run: `{source_run}`\n")
    fh.write(f"- parity_root: `{parity_root}`\n")
    fh.write(f"- cases: `{len(cases)}`\n")
    fh.write(f"- max_abs_force_error_eV_per_A: `{data['max_abs_error_eV_per_A']:.9e}`\n\n")
    fh.write("| stage | bead | atom | component | abs_error_eV_per_A | signed_error_eV_per_A | jax_force | oracle_force |\n")
    fh.write("| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |\n")
    for row in rows:
        fh.write(
            f"| {row['stage']} | {row['bead']} | {row['atom']} | {row['component']} | "
            f"{row['abs_error']:.9e} | {row['signed_error']:.9e} | {row['jax_force']:.9e} | {row['oracle_force']:.9e} |\n"
        )

max_error = data["max_abs_error_eV_per_A"]
print(json.dumps({"cases": len(cases), "max_abs_force_error_eV_per_A": max_error, "summary": str(summary)}, sort_keys=True))
if max_force_error_arg:
    threshold = float(max_force_error_arg)
    if max_error > threshold:
        raise SystemExit(f"max force error {max_error:.9e} exceeds threshold {threshold:.9e}")
PY

echo "summary: $SUMMARY"
echo "dataset: $DATASET"
echo "report: $REPORT"
