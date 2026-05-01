#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY="$ROOT_DIR/build/test_summaries/test_11_structural_stress.summary"
mkdir -p "$(dirname "$SUMMARY")"

"$ROOT_DIR/build/test_11_structural_stress" | tee "$SUMMARY"

set +e
"$ROOT_DIR/build/test_11_invalid_group_overlap" > "$ROOT_DIR/build/test_summaries/test_11_invalid_group_overlap.summary" 2>&1
status=$?
set -e

if [[ "$status" -eq 0 ]]; then
  echo "test_11_invalid_group_overlap: FAIL - overlapping groups unexpectedly succeeded" | tee -a "$SUMMARY"
  exit 1
fi

if ! grep -q "particle matches multiple groups" "$ROOT_DIR/build/test_summaries/test_11_invalid_group_overlap.summary"; then
  echo "test_11_invalid_group_overlap: FAIL - expected overlap error not found" | tee -a "$SUMMARY"
  cat "$ROOT_DIR/build/test_summaries/test_11_invalid_group_overlap.summary"
  exit 1
fi

echo "test_11_invalid_group_overlap: PASS" | tee -a "$SUMMARY"
