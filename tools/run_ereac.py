#!/usr/bin/env python3
"""Run the rebuilt ereac in a case directory and fail loudly if the run did not really happen.

A run counts as successful only if: exit status 0, 'Normal end of MD-simulation' in the output,
no 'Energy not a number', and the expected output files are newly written (stale ones are removed first).
Start-up failures (non-zero exit, e.g. the intermittent dyld error) are retried.

Usage: run_ereac.py <ereac_binary> <case_dir> [--expect fort.73 ereac_export.out] [--tries 6]
"""
import argparse
import pathlib
import subprocess
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("case")
    ap.add_argument("--expect", nargs="*", default=["fort.73"])
    ap.add_argument("--tries", type=int, default=10)
    ap.add_argument("--timeout", type=int, default=600)
    a = ap.parse_args()
    case = pathlib.Path(a.case).resolve()
    binary = pathlib.Path(a.binary).resolve()
    for need in ("control", "fort.3", "fort.4", "fort.20", "fort.35"):
        if not (case / need).exists():
            sys.exit(f"missing input {need} in {case}")
    last = ""
    for attempt in range(1, a.tries + 1):
        for f in a.expect + ["run.out", "run.err"]:
            (case / f).unlink(missing_ok=True)
        t0 = time.time()
        with open(case / "run.out", "w") as out, open(case / "run.err", "w") as err:
            r = subprocess.run([str(binary)], cwd=case, stdout=out, stderr=err, timeout=a.timeout)
        text = (case / "run.out").read_text() + (case / "run.err").read_text()
        ok = (r.returncode == 0 and "Normal end of MD-simulation" in text and "Energy not a number" not in text
              and all((case / f).exists() and (case / f).stat().st_mtime >= t0 for f in a.expect))
        if ok:
            print(f"OK  {case.name}  attempt {attempt}  {time.time() - t0:.1f}s")
            return
        last = f"exit {r.returncode}: " + text.strip().splitlines()[-1][:160] if text.strip() else f"exit {r.returncode}"
        print(f"retry {attempt}/{a.tries}: {last}", file=sys.stderr)
        time.sleep(min(10, 1.5 * attempt))  # start-up failures cluster under machine load
    sys.exit(f"FAILED {case.name}: {last}")


if __name__ == "__main__":
    main()
