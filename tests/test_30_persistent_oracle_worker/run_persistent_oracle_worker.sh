#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORK_ROOT="$ROOT_DIR/build/test_30_persistent_worker"
mkdir -p "$WORK_ROOT"

export EREAXFF_V1_ROOT="${EREAXFF_V1_ROOT:-/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1}"
export EREAXFF_V2_ROOT="${EREAXFF_V2_ROOT:-/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2}"
export EREAXFF_ADAPTER_WORK_ROOT="$WORK_ROOT"
export EREAXFF_ORACLE_SMOKE_BINARY="$ROOT_DIR/build/oracle_single_point_force_export"
export EREAXFF_PERSISTENT_WORKER_BINARY="$ROOT_DIR/build/persistent_oracle_worker"

rm -rf "$WORK_ROOT"
mkdir -p "$WORK_ROOT"

"$ROOT_DIR/build/test_30_persistent_oracle_worker_bin"
