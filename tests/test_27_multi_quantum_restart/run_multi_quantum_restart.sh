#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_27_multi_quantum_restart.summary"
WORK_DIR="$ROOT_DIR/build/test_27_multi_quantum_restart"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

"$ROOT_DIR/build/test_27_multi_quantum_restart_bin" "$WORK_DIR/restart.dat" > "$WORK_DIR/stdout.log"

if [[ ! -s "$WORK_DIR/restart.dat" ]]; then
  echo "test_27_multi_quantum_restart: FAIL missing restart" | tee "$SUMMARY"
  exit 1
fi
if ! grep -q "nquantum   = 2" "$WORK_DIR/stdout.log"; then
  echo "test_27_multi_quantum_restart: FAIL missing nquantum diagnostic" | tee "$SUMMARY"
  cat "$WORK_DIR/stdout.log"
  exit 1
fi

{
  echo "test_27_multi_quantum_restart: PASS"
  echo "  nquantum = 2"
  echo "  restart = $WORK_DIR/restart.dat"
} | tee "$SUMMARY"
