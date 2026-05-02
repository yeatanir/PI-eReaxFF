#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
BIN="$ROOT_DIR/build/test_12_oracle_force_manager_p1"
SUMMARY="$SUMMARY_DIR/test_12_oracle_force_manager_p1.summary"

mkdir -p "$SUMMARY_DIR"

if [[ ! -x "$BIN" ]]; then
  echo "test_12_oracle_force_manager_p1: FAIL" > "$SUMMARY"
  echo "Missing executable: $BIN" >> "$SUMMARY"
  cat "$SUMMARY"
  exit 1
fi

"$BIN" | tee "$SUMMARY"
