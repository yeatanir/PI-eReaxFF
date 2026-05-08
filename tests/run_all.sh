#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

make build_tests

mkdir -p build/test_summaries
./build/test_00_build | tee build/test_summaries/test_00_build.summary
./build/test_01_normal_modes | tee build/test_summaries/test_01_normal_modes.summary
./build/test_02_free_ring_polymer | tee build/test_summaries/test_02_free_ring_polymer.summary
./build/test_03_harmonic_quantum_variance | tee build/test_summaries/test_03_harmonic_quantum_variance.summary
./build/test_04_thermostat_ou_statistics | tee build/test_summaries/test_04_thermostat_ou_statistics.summary
./build/test_05_pile_mode_temperatures | tee build/test_summaries/test_05_pile_mode_temperatures.summary

set +e
./build/test_05_pile_g_stub > build/test_summaries/test_05_pile_g_stub.summary 2>&1
stub_status=$?
set -e
if [[ "$stub_status" -eq 0 ]]; then
  echo "test_05_pile_g_stub: FAIL - PILE-G stub unexpectedly succeeded"
  exit 1
fi
if ! grep -q "PILE-G requested but not implemented" build/test_summaries/test_05_pile_g_stub.summary; then
  echo "test_05_pile_g_stub: FAIL - expected PILE-G error message not found"
  cat build/test_summaries/test_05_pile_g_stub.summary
  exit 1
fi
echo "test_05_pile_g_stub: PASS"
./build/test_06_force_scaling | tee build/test_summaries/test_06_force_scaling.summary
./build/test_07_charge_memory | tee build/test_summaries/test_07_charge_memory.summary
cat build/test_summaries/test_08_reaxff_oracle_smoke.summary
./build/test_09_mixed_classical_pimd_config | tee build/test_summaries/test_09_mixed_classical_pimd_config.summary
cat build/test_summaries/test_10_harmonic_driver_outputs.summary
cat build/test_summaries/test_11_structural_stress.summary
cat build/test_summaries/test_12_oracle_force_manager_p1.summary
cat build/test_summaries/test_13_oracle_charge_policy_modes.summary
cat build/test_summaries/test_14_mixed_oracle_manager_fixture.summary
cat build/test_summaries/test_15_mixed_oracle_manager_matrix.summary
cat build/test_summaries/test_16_mixed_oracle_manager_restart_output.summary
cat build/test_summaries/test_17_ereaxff_mixed_driver_smoke.summary
cat build/test_summaries/test_18_ereaxff_mixed_baoab_matrix.summary
cat build/test_summaries/test_19_ereaxff_mixed_baoab_restart.summary
cat build/test_summaries/test_20_ereaxff_mixed_baoab_policy_matrix.summary
cat build/test_summaries/test_21_ereaxff_mixed_baoab_config_smoke.summary
cat build/test_summaries/test_22_ereaxff_mixed_baoab_config_matrix.summary
cat build/test_summaries/test_23_ereaxff_mixed_baoab_geometry_matrix.summary
cat build/test_summaries/test_24_ereaxff_mixed_baoab_user_geo_matrix.summary
cat build/test_summaries/test_25_phase6_like_temperature_partition.summary
cat build/test_summaries/test_26_crystx_periodic_loader.summary
cat build/test_summaries/test_27_multi_quantum_restart.summary
cat build/test_summaries/test_28_paired_graphene_60el_60c_smoke.summary

echo "All implemented tests passed."
