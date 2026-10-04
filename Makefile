BUILD_DIR ?= build/release
BUILD_TYPE ?= Release
JOBS ?= 4

.PHONY: all configure cpu-serial benchmark-backends test clean

all: configure
	cmake --build $(BUILD_DIR) -j$(JOBS)

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

cpu-serial: configure
	cmake --build $(BUILD_DIR) --target schrodinger_helmholtz_serial -j$(JOBS)

benchmark-backends: configure
	cmake --build $(BUILD_DIR) --target sh1d_benchmark_cpu_serial sh1d_benchmark_cpu -j$(JOBS)

test: all
	ctest --test-dir $(BUILD_DIR) --output-on-failure

clean:
	cmake -E remove_directory build
