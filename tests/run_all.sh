#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

make build_tests

mkdir -p build/test_summaries
./build/test_00_build | tee build/test_summaries/test_00_build.summary
./build/test_01_normal_modes | tee build/test_summaries/test_01_normal_modes.summary
./build/test_02_free_ring_polymer | tee build/test_summaries/test_02_free_ring_polymer.summary

echo "All implemented tests passed."
