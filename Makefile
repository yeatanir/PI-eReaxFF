override FC := $(if $(filter default undefined,$(origin FC)),gfortran,$(FC))
FFLAGS ?= -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -O2

BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod

CORE_SRCS := src/core/kinds.f90 src/core/constants.f90
PIMD_SRCS := src/pimd/normal_modes.f90
TEST01_SRCS := tests/test_01_normal_modes/test_normal_modes.f90

.PHONY: all test test_01_normal_modes clean

all: test_01_normal_modes

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

test_01_normal_modes: $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) $(CORE_SRCS) $(PIMD_SRCS) $(TEST01_SRCS) -o $(BUILD_DIR)/test_01_normal_modes

test: test_01_normal_modes
	./tests/run_all.sh

clean:
	rm -rf $(BUILD_DIR)
