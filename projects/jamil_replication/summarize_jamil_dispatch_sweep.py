#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import pathlib
import statistics
import sys
from dataclasses import dataclass


@dataclass
class CaseSummary:
    case_name: str
    omp_threads: int
    dispatch_mode: str
    max_workers: int
    total_wall_s: float
    n_oracle_calls: int
    mean_call_s: float
    min_call_s: float
    max_call_s: float
    p50_call_s: float
    p90_call_s: float
    p95_call_s: float
    status: str
    notes: str


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return math.nan
    values = sorted(values)
    if len(values) == 1:
        return values[0]
    rank = (pct / 100.0) * (len(values) - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return values[lo]
    frac = rank - lo
    return values[lo] * (1.0 - frac) + values[hi] * frac


def parse_timing_log(path: pathlib.Path) -> tuple[list[float], int]:
    calls: list[float] = []
    batch_count = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if parts[0] == "BATCH":
            batch_count += 1
        elif parts[0] == "CALL" and len(parts) >= 15:
            calls.append(float(parts[12]))
    return calls, batch_count


def summarize_case(case_dir: pathlib.Path) -> CaseSummary | None:
    result_path = case_dir / "result.json"
    if not result_path.exists():
        return None
    result = json.loads(result_path.read_text(encoding="utf-8"))
    timing_path = case_dir / "out" / "adapter_work" / "oracle_timing.log"
    calls: list[float] = []
    if timing_path.exists():
        calls, _ = parse_timing_log(timing_path)
    status = "PASS" if result.get("status") == 0 else "FAIL"
    notes = result.get("note", "")
    if not calls:
        notes = f"{notes} missing_timing_log".strip()
    return CaseSummary(
        case_name=case_dir.name,
        omp_threads=int(result.get("omp_threads", 0)),
        dispatch_mode=extract_dispatch_mode(case_dir / "config.yaml"),
        max_workers=extract_max_workers(case_dir / "config.yaml"),
        total_wall_s=float(result.get("wall_s", math.nan)),
        n_oracle_calls=len(calls),
        mean_call_s=statistics.fmean(calls) if calls else math.nan,
        min_call_s=min(calls) if calls else math.nan,
        max_call_s=max(calls) if calls else math.nan,
        p50_call_s=percentile(calls, 50.0),
        p90_call_s=percentile(calls, 90.0),
        p95_call_s=percentile(calls, 95.0),
        status=status,
        notes=notes,
    )


def extract_max_workers(path: pathlib.Path) -> int:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("oracle_max_workers:"):
            return int(line.split(":", 1)[1].strip())
    return 0


def extract_dispatch_mode(path: pathlib.Path) -> str:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("oracle_dispatch_mode:"):
            return line.split(":", 1)[1].strip()
    return ""


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: summarize_jamil_dispatch_sweep.py <sweep_root> <summary_out>", file=sys.stderr)
        return 2
    sweep_root = pathlib.Path(sys.argv[1])
    summary_out = pathlib.Path(sys.argv[2])
    rows: list[CaseSummary] = []
    for case_dir in sorted(p for p in sweep_root.iterdir() if p.is_dir()):
        summary = summarize_case(case_dir)
        if summary is not None:
            rows.append(summary)

    summary_out.parent.mkdir(parents=True, exist_ok=True)
    with summary_out.open("w", encoding="utf-8") as handle:
        handle.write("# Jamil Dispatch Sweep Summary\n\n")
        handle.write("| case_name | OMP_NUM_THREADS | oracle_dispatch_mode | oracle_max_workers | total_wall_s | n_oracle_calls | mean_call_s | min_call_s | max_call_s | p50_call_s | p90_call_s | p95_call_s | status | notes |\n")
        handle.write("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |\n")
        for row in rows:
            handle.write(
                f"| {row.case_name} | {row.omp_threads} | {row.dispatch_mode} | {row.max_workers} | "
                f"{row.total_wall_s:.6f} | {row.n_oracle_calls} | {row.mean_call_s:.6f} | {row.min_call_s:.6f} | "
                f"{row.max_call_s:.6f} | {row.p50_call_s:.6f} | {row.p90_call_s:.6f} | {row.p95_call_s:.6f} | "
                f"{row.status} | {row.notes} |\n"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
