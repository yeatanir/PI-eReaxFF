override FC := $(if $(filter default undefined,$(origin FC)),gfortran,$(FC))
FFLAGS ?= -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -O2

BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod

CORE_SRCS := src/core/kinds.f90 src/core/constants.f90 src/core/units.f90 src/core/errors.f90 src/core/rng.f90
PIMD_SRCS := src/pimd/normal_modes.f90 src/pimd/ring_polymer.f90 src/pimd/spring_energy.f90 src/pimd/free_ring_propagator.f90
THERMOSTAT_SRCS := src/thermostat/thermostat_base.f90 src/thermostat/thermostat_none.f90 src/thermostat/thermostat_langevin.f90 src/thermostat/thermostat_pile_l.f90 src/thermostat/thermostat_pile_g.f90 src/thermostat/thermostat_registry.f90
FORCE_SRCS := src/force/force_provider_base.f90 src/force/dummy_force_provider.f90 src/force/harmonic_force_provider.f90 src/force/charge_coupling_policy.f90 src/force/bead_force_manager.f90 src/force/ereaxff_oracle_adapter.f90
TEST00_SRCS := tests/test_00_build/test_constants_units.f90
TEST01_SRCS := tests/test_01_normal_modes/test_normal_modes.f90
TEST02_SRCS := tests/test_02_free_ring_polymer/test_free_ring_polymer.f90
TEST04_SRCS := tests/test_04_thermostat_ou_statistics/test_thermostat_ou_statistics.f90
TEST05_SRCS := tests/test_05_pile_mode_temperatures/test_pile_mode_temperatures.f90
TEST05_STUB_SRCS := tests/test_05_pile_mode_temperatures/test_pile_g_stub.f90
TEST06_SRCS := tests/test_06_force_scaling/test_force_scaling.f90
TEST07_SRCS := tests/test_07_charge_memory/test_charge_memory.f90

.PHONY: all test build_tests test_00_build test_01_normal_modes test_02_free_ring_polymer test_04_thermostat_ou_statistics test_05_pile_mode_temperatures test_05_pile_g_stub test_06_force_scaling test_07_charge_memory clean

all: build_tests

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

test_00_build: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(TEST00_SRCS) -o $(BUILD_DIR)/test_00_build

test_01_normal_modes: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST01_SRCS) -o $(BUILD_DIR)/test_01_normal_modes

test_02_free_ring_polymer: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST02_SRCS) -o $(BUILD_DIR)/test_02_free_ring_polymer

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

build_tests: test_00_build test_01_normal_modes test_02_free_ring_polymer test_04_thermostat_ou_statistics test_05_pile_mode_temperatures test_05_pile_g_stub test_06_force_scaling test_07_charge_memory

test:
	./tests/run_all.sh

clean:
	rm -rf $(BUILD_DIR)
