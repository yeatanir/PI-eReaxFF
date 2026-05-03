#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUMMARY_DIR="$ROOT_DIR/build/test_summaries"
SUMMARY="$SUMMARY_DIR/test_26_crystx_periodic_loader.summary"
WORK_DIR="$ROOT_DIR/build/test_26_crystx_periodic_loader_work"
BIN="$ROOT_DIR/build/test_26_crystx_periodic_loader_bin"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR" "$SUMMARY_DIR"

geo="$WORK_DIR/crystx.geo"
out="$WORK_DIR/out"
mkdir -p "$out"

cat > "$geo" <<'EOF'
BIOGRF 200
DESCRP periodic_crystx_min_image
RUTYPE NORMAL RUN
CRYSTX    10.00000   10.00000   10.00000   90.00000   90.00000   90.00000
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
HETATM     1 C           C       0.10000   5.00000   5.00000        0 0  0.00000
HETATM     2 C           C       9.90000   5.00000   5.00000        0 0  0.00000
EOF

"$BIN" "$geo" "$out" > "$out/stdout.log"

for file in fort.7 fort.8 molfra.out xmolout; do
  if [[ ! -s "$out/$file" ]]; then
    echo "test_26_crystx_periodic_loader: FAIL missing $file" | tee "$SUMMARY"
    exit 1
  fi
done

if ! grep -q "cell_defined=T" "$out/stdout.log"; then
  echo "test_26_crystx_periodic_loader: FAIL cell not reported as defined" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "cell_periodic=T" "$out/stdout.log"; then
  echo "test_26_crystx_periodic_loader: FAIL cell not reported as periodic" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "10.00000" "$out/stdout.log"; then
  echo "test_26_crystx_periodic_loader: FAIL cell lengths not echoed" | tee "$SUMMARY"
  cat "$out/stdout.log"
  exit 1
fi
if ! grep -q "#Bonds: 1" "$out/fort.7"; then
  echo "test_26_crystx_periodic_loader: FAIL minimum-image bond not detected" | tee "$SUMMARY"
  cat "$out/fort.7"
  exit 1
fi
if ! grep -q "molecule list per frame" "$out/molfra.out"; then
  echo "test_26_crystx_periodic_loader: FAIL molecule list missing" | tee "$SUMMARY"
  exit 1
fi

{
  echo "test_26_crystx_periodic_loader: PASS"
  echo "  geo = $geo"
  echo "  periodic_bond = detected"
} | tee "$SUMMARY"
