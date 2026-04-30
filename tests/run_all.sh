#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

make test_01_normal_modes

mkdir -p build/test_summaries
./build/test_01_normal_modes | tee build/test_summaries/test_01_normal_modes.summary

echo "All implemented tests passed."

