#!/usr/bin/env python3
"""Run ereac in a case dir, retrying ONLY when the process failed to start (macOS dyld error).
Any other exit status is reported, not retried (minimization runs end with 'Read past ENDFILE', exit 2).
Usage: ereac_retry.py <binary> <case_dir>   -> prints 'started try N rc R'; stdout/stderr go to run.out/run.err
"""
import pathlib, subprocess, sys, time
binary, case = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve()
for i in range(1, 61):
    r = subprocess.run([str(binary)], cwd=case, capture_output=True, text=True, timeout=7200)
    if "dyld" not in r.stderr:
        (case / "run.out").write_text(r.stdout); (case / "run.err").write_text(r.stderr)
        print(f"started try {i} rc {r.returncode}"); sys.exit(0)
    time.sleep(min(5, 1 + 0.2 * i))
sys.exit("could not start ereac in 60 tries (dyld failures)")
