override FC := $(if $(filter default undefined,$(origin FC)),gfortran,$(FC))
FFLAGS ?= -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -O2

BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod

CORE_SRCS := src/core/kinds.f90 src/core/constants.f90 src/core/units.f90 src/core/errors.f90 src/core/rng.f90 src/core/system_state.f90 src/core/selectors.f90 src/core/particle_groups.f90 src/core/restart.f90
PIMD_SRCS := src/pimd/normal_modes.f90 src/pimd/ring_polymer.f90 src/pimd/spring_energy.f90 src/pimd/free_ring_propagator.f90 src/pimd/harmonic_quantum_variance.f90
THERMOSTAT_SRCS := src/thermostat/thermostat_base.f90 src/thermostat/thermostat_none.f90 src/thermostat/thermostat_langevin.f90 src/thermostat/thermostat_pile_l.f90 src/thermostat/thermostat_pile_g.f90 src/thermostat/thermostat_registry.f90
FORCE_SRCS := src/force/force_provider_base.f90 src/force/dummy_force_provider.f90 src/force/harmonic_force_provider.f90 src/force/charge_coupling_policy.f90 src/force/bead_force_manager.f90 src/force/ereaxff_oracle_adapter.f90
INTEGRATOR_SRCS := src/integrator/baoab_integrator.f90
OUTPUT_SRCS := src/output/output_manager.f90 src/output/pimd_output.f90 src/output/diagnostics.f90 src/output/reaxff_compat_output.f90
HARMONIC_DRIVER_SRCS := src/drivers/harmonic_pimd_main.f90
TEST00_SRCS := tests/test_00_build/test_constants_units.f90
TEST01_SRCS := tests/test_01_normal_modes/test_normal_modes.f90
TEST02_SRCS := tests/test_02_free_ring_polymer/test_free_ring_polymer.f90
TEST03_SRCS := tests/test_03_harmonic_quantum_variance/test_harmonic_quantum_variance.f90
TEST04_SRCS := tests/test_04_thermostat_ou_statistics/test_thermostat_ou_statistics.f90
TEST05_SRCS := tests/test_05_pile_mode_temperatures/test_pile_mode_temperatures.f90
TEST05_STUB_SRCS := tests/test_05_pile_mode_temperatures/test_pile_g_stub.f90
TEST06_SRCS := tests/test_06_force_scaling/test_force_scaling.f90
TEST07_SRCS := tests/test_07_charge_memory/test_charge_memory.f90
TEST08_SRCS := tests/test_08_reaxff_oracle_smoke/test_oracle_adapter.f90
TEST09_SRCS := tests/test_09_mixed_classical_pimd_config/test_mixed_classical_pimd_config.f90
TEST12_SRCS := tests/test_12_oracle_force_manager_p1/test_oracle_force_manager_p1.f90
TEST13_SRCS := tests/test_13_oracle_charge_policy_modes/test_oracle_charge_policy_modes.f90
TEST14_SRCS := tests/test_14_mixed_oracle_manager_fixture/test_mixed_oracle_manager_fixture.f90
TEST11_SRCS := tests/test_11_structural_stress/test_structural_stress.f90
TEST11_NEG_SRCS := tests/test_11_structural_stress/test_invalid_group_overlap.f90

.PHONY: all test build_tests harmonic_pimd_main test_00_build test_01_normal_modes test_02_free_ring_polymer test_03_harmonic_quantum_variance test_04_thermostat_ou_statistics test_05_pile_mode_temperatures test_05_pile_g_stub test_06_force_scaling test_07_charge_memory test_08_reaxff_oracle_smoke test_09_mixed_classical_pimd_config test_10_driver_outputs test_11_structural_stress test_12_oracle_force_manager_p1 test_13_oracle_charge_policy_modes test_14_mixed_oracle_manager_fixture clean

all: build_tests harmonic_pimd_main

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

test_00_build: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(TEST00_SRCS) -o $(BUILD_DIR)/test_00_build

test_01_normal_modes: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST01_SRCS) -o $(BUILD_DIR)/test_01_normal_modes

test_02_free_ring_polymer: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST02_SRCS) -o $(BUILD_DIR)/test_02_free_ring_polymer

test_03_harmonic_quantum_variance: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(THERMOSTAT_SRCS) $(INTEGRATOR_SRCS) $(TEST03_SRCS) -o $(BUILD_DIR)/test_03_harmonic_quantum_variance

test_04_thermostat_ou_statistics: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(THERMOSTAT_SRCS) $(TEST04_SRCS) -o $(BUILD_DIR)/test_04_thermostat_ou_statistics

test_05_pile_mode_temperatures: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(THERMOSTAT_SRCS) $(TEST05_SRCS) -o $(BUILD_DIR)/test_05_pile_mode_temperatures

test_05_pile_g_stub: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(THERMOSTAT_SRCS) $(TEST05_STUB_SRCS) -o $(BUILD_DIR)/test_05_pile_g_stub

test_06_force_scaling: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(FORCE_SRCS) $(TEST06_SRCS) -o $(BUILD_DIR)/test_06_force_scaling

test_07_charge_memory: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(FORCE_SRCS) $(TEST07_SRCS) -o $(BUILD_DIR)/test_07_charge_memory

test_08_reaxff_oracle_smoke: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(FORCE_SRCS) $(TEST08_SRCS) -o $(BUILD_DIR)/test_08_reaxff_oracle_smoke
	tests/test_08_reaxff_oracle_smoke/run_v2_oracle_smoke.sh

test_09_mixed_classical_pimd_config: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(TEST09_SRCS) -o $(BUILD_DIR)/test_09_mixed_classical_pimd_config

test_12_oracle_force_manager_p1: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(FORCE_SRCS) $(TEST12_SRCS) -o $(BUILD_DIR)/test_12_oracle_force_manager_p1
	tests/test_12_oracle_force_manager_p1/run_oracle_force_manager_p1.sh

test_13_oracle_charge_policy_modes: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(FORCE_SRCS) $(TEST13_SRCS) -o $(BUILD_DIR)/test_13_oracle_charge_policy_modes
	tests/test_13_oracle_charge_policy_modes/run_oracle_charge_policy_modes.sh

test_14_mixed_oracle_manager_fixture: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(FORCE_SRCS) $(TEST14_SRCS) -o $(BUILD_DIR)/test_14_mixed_oracle_manager_fixture
	tests/test_14_mixed_oracle_manager_fixture/run_mixed_oracle_manager_fixture.sh

harmonic_pimd_main: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(THERMOSTAT_SRCS) $(INTEGRATOR_SRCS) $(OUTPUT_SRCS) $(HARMONIC_DRIVER_SRCS) -o $(BUILD_DIR)/harmonic_pimd_main

test_10_driver_outputs: harmonic_pimd_main
	tests/test_10_reaxff_compat_outputs/run_harmonic_driver_smoke.sh

test_11_structural_stress: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(FORCE_SRCS) $(TEST11_SRCS) -o $(BUILD_DIR)/test_11_structural_stress
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(TEST11_NEG_SRCS) -o $(BUILD_DIR)/test_11_invalid_group_overlap
	tests/test_11_structural_stress/run_structural_stress.sh

build_tests: test_00_build test_01_normal_modes test_02_free_ring_polymer test_03_harmonic_quantum_variance test_04_thermostat_ou_statistics test_05_pile_mode_temperatures test_05_pile_g_stub test_06_force_scaling test_07_charge_memory test_08_reaxff_oracle_smoke test_09_mixed_classical_pimd_config test_10_driver_outputs test_11_structural_stress test_12_oracle_force_manager_p1 test_13_oracle_charge_policy_modes test_14_mixed_oracle_manager_fixture

test:
	./tests/run_all.sh

clean:
	rm -rf $(BUILD_DIR)
