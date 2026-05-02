#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
WORK_ROOT="$ROOT_DIR/build/test_adapter_work"
BIN="$ROOT_DIR/build/test_08_reaxff_oracle_smoke"
SUMMARY="$SUMMARY_DIR/test_08_reaxff_oracle_smoke.summary"

mkdir -p "$SUMMARY_DIR" "$WORK_ROOT"

fail() {
  {
    echo "test_08_reaxff_oracle_smoke: FAIL"
    echo "$1"
  } > "$SUMMARY"
  cat "$SUMMARY"
  exit 1
}

if [[ ! -x "$BIN" ]]; then
  fail "Missing executable: $BIN"
fi

export EREAXFF_ADAPTER_WORK_ROOT="$WORK_ROOT"

"$BIN" smoke > "$SUMMARY" 2>&1 || fail "Adapter smoke mode failed"
if ! grep -q "test_08_reaxff_oracle_smoke: PASS" "$SUMMARY"; then
  fail "Smoke output did not report PASS"
fi

set +e
"$BIN" nan > "$WORK_ROOT/test_08_nan.log" 2>&1
nan_status=$?
set -e
if [[ "$nan_status" -eq 0 ]]; then
  fail "NaN rejection mode unexpectedly succeeded"
fi
if ! grep -q "oracle_adapter_eval: non-finite position" "$WORK_ROOT/test_08_nan.log"; then
  fail "NaN rejection did not report the expected non-finite position error"
fi

echo "test_08_reaxff_oracle_smoke: PASS"
cat "$SUMMARY"
