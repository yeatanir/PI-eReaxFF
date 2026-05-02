#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
BIN="$ROOT_DIR/build/test_13_oracle_charge_policy_modes"
SUMMARY="$SUMMARY_DIR/test_13_oracle_charge_policy_modes.summary"

mkdir -p "$SUMMARY_DIR"

if [[ ! -x "$BIN" ]]; then
  echo "test_13_oracle_charge_policy_modes: FAIL" > "$SUMMARY"
  echo "Missing executable: $BIN" >> "$SUMMARY"
  cat "$SUMMARY"
  exit 1
fi

"$BIN" | tee "$SUMMARY"
