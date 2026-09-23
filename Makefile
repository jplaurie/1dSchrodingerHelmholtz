BUILD_DIR ?= build/release
BUILD_TYPE ?= Release
JOBS ?= 4

.PHONY: all configure test clean

all: configure
	cmake --build $(BUILD_DIR) -j$(JOBS)

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

test: all
	ctest --test-dir $(BUILD_DIR) --output-on-failure

clean:
	cmake -E remove_directory build

