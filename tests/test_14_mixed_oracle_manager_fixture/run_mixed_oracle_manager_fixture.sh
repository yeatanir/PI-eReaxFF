#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
BIN="$ROOT_DIR/build/test_14_mixed_oracle_manager_fixture"
SUMMARY="$SUMMARY_DIR/test_14_mixed_oracle_manager_fixture.summary"

mkdir -p "$SUMMARY_DIR"

if [[ ! -x "$BIN" ]]; then
  echo "test_14_mixed_oracle_manager_fixture: FAIL" > "$SUMMARY"
  echo "Missing executable: $BIN" >> "$SUMMARY"
  cat "$SUMMARY"
  exit 1
fi

"$BIN" | tee "$SUMMARY"
