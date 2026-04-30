override FC := $(if $(filter default undefined,$(origin FC)),gfortran,$(FC))
FFLAGS ?= -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -O2

BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod

CORE_SRCS := src/core/kinds.f90 src/core/constants.f90 src/core/units.f90 src/core/errors.f90 src/core/rng.f90
PIMD_SRCS := src/pimd/normal_modes.f90 src/pimd/ring_polymer.f90 src/pimd/spring_energy.f90 src/pimd/free_ring_propagator.f90
TEST00_SRCS := tests/test_00_build/test_constants_units.f90
TEST01_SRCS := tests/test_01_normal_modes/test_normal_modes.f90
TEST02_SRCS := tests/test_02_free_ring_polymer/test_free_ring_polymer.f90

.PHONY: all test build_tests test_00_build test_01_normal_modes test_02_free_ring_polymer clean

all: build_tests

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

test_00_build: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(TEST00_SRCS) -o $(BUILD_DIR)/test_00_build

test_01_normal_modes: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST01_SRCS) -o $(BUILD_DIR)/test_01_normal_modes

test_02_free_ring_polymer: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST02_SRCS) -o $(BUILD_DIR)/test_02_free_ring_polymer

build_tests: test_00_build test_01_normal_modes test_02_free_ring_polymer

test:
	./tests/run_all.sh

clean:
	rm -rf $(BUILD_DIR)
