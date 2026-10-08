cmake_minimum_required(VERSION 3.20)
project(SchrodingerHelmholtz1D VERSION 0.3.0 LANGUAGES C CXX)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Git QUIET)
set(SH1D_GIT_COMMIT "unknown")
set(SH1D_GIT_DIRTY "unknown")
if(Git_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE SH1D_GIT_COMMIT_CANDIDATE OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE SH1D_GIT_COMMIT_RESULT
    ERROR_QUIET)
  if(SH1D_GIT_COMMIT_RESULT EQUAL 0)
    set(SH1D_GIT_COMMIT "${SH1D_GIT_COMMIT_CANDIDATE}")
  endif()
  execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE SH1D_GIT_STATUS OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE SH1D_GIT_STATUS_RESULT
    ERROR_QUIET)
  if(SH1D_GIT_STATUS_RESULT EQUAL 0)
    if(SH1D_GIT_STATUS STREQUAL "")
      set(SH1D_GIT_DIRTY "false")
    else()
      set(SH1D_GIT_DIRTY "true")
    endif()
  endif()
endif()

option(SH1D_OPENMP "Enable OpenMP and threaded FFTW" ON)
option(SH1D_HDF5 "Enable optional HDF5 field snapshots and exporter" ON)

function(sh1d_enable_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

find_path(FFTW_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW_LIBRARY NAMES fftw3 REQUIRED)
find_library(FFTW_THREADS_LIBRARY NAMES fftw3_threads)
if(SH1D_HDF5)
  find_package(HDF5 COMPONENTS C QUIET)
endif()

set(SH1D_COMMON_SOURCES
    src/parameters.cpp
    src/fftw_utils.cpp
    src/hdf5_io.cpp
    src/spectral.cpp
    src/nonlinear.cpp
    src/diagnostics.cpp
    src/output.cpp
    src/output_transaction.cpp
    src/solver.cpp)
add_library(sh1d_common ${SH1D_COMMON_SOURCES})
target_include_directories(sh1d_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR})
target_link_libraries(sh1d_common PUBLIC ${FFTW_LIBRARY})
target_compile_features(sh1d_common PUBLIC cxx_std_20)
target_compile_definitions(sh1d_common PUBLIC
  SH1D_VERSION="${PROJECT_VERSION}" SH1D_GIT_COMMIT="${SH1D_GIT_COMMIT}"
  SH1D_GIT_DIRTY="${SH1D_GIT_DIRTY}")
if(SH1D_HDF5 AND HDF5_FOUND)
  target_include_directories(sh1d_common PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(sh1d_common PUBLIC ${HDF5_LIBRARIES})
  target_compile_definitions(sh1d_common PUBLIC SH1D_HAVE_HDF5=1)
elseif(SH1D_HDF5)
  message(WARNING "HDF5 not found; HDF5 field output and exporter will not be available")
endif()

if(SH1D_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(sh1d_common PUBLIC OpenMP::OpenMP_CXX)
  endif()
  if(FFTW_THREADS_LIBRARY)
    target_link_libraries(sh1d_common PUBLIC ${FFTW_THREADS_LIBRARY})
    target_compile_definitions(sh1d_common PUBLIC SH1D_HAVE_FFTW_THREADS=1)
  else()
    message(WARNING "FFTW threads library not found; FFT plans remain single-threaded")
  endif()
endif()

sh1d_enable_warnings(sh1d_common)

set(SH1D_SERIAL_COMMON sh1d_common)
if(SH1D_OPENMP AND OpenMP_CXX_FOUND)
  add_library(sh1d_common_serial ${SH1D_COMMON_SOURCES})
  target_include_directories(sh1d_common_serial
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR})
  target_link_libraries(sh1d_common_serial PUBLIC ${FFTW_LIBRARY})
  target_compile_features(sh1d_common_serial PUBLIC cxx_std_20)
  target_compile_definitions(sh1d_common_serial PUBLIC
    SH1D_VERSION="${PROJECT_VERSION}" SH1D_GIT_COMMIT="${SH1D_GIT_COMMIT}"
    SH1D_GIT_DIRTY="${SH1D_GIT_DIRTY}")
  if(SH1D_HDF5 AND HDF5_FOUND)
    target_include_directories(sh1d_common_serial PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(sh1d_common_serial PUBLIC ${HDF5_LIBRARIES})
    target_compile_definitions(sh1d_common_serial PUBLIC SH1D_HAVE_HDF5=1)
  endif()
  sh1d_enable_warnings(sh1d_common_serial)
  set(SH1D_SERIAL_COMMON sh1d_common_serial)
endif()

add_executable(schrodinger_helmholtz src/main.cpp)
target_link_libraries(schrodinger_helmholtz PRIVATE sh1d_common)
sh1d_enable_warnings(schrodinger_helmholtz)

add_executable(schrodinger_helmholtz_serial src/main.cpp)
target_link_libraries(schrodinger_helmholtz_serial PRIVATE ${SH1D_SERIAL_COMMON})
sh1d_enable_warnings(schrodinger_helmholtz_serial)

add_executable(sh1d_benchmark_cpu src/benchmark.cpp)
target_link_libraries(sh1d_benchmark_cpu PRIVATE sh1d_common)
sh1d_enable_warnings(sh1d_benchmark_cpu)

add_executable(sh1d_benchmark_cpu_serial src/benchmark.cpp)
target_link_libraries(sh1d_benchmark_cpu_serial PRIVATE ${SH1D_SERIAL_COMMON})
target_compile_definitions(sh1d_benchmark_cpu_serial PRIVATE SH1D_BENCHMARK_SERIAL=1)
sh1d_enable_warnings(sh1d_benchmark_cpu_serial)

if(SH1D_HDF5 AND HDF5_FOUND)
  add_executable(sh1d_hdf5_export src/hdf5_export.cpp)
  target_link_libraries(sh1d_hdf5_export PRIVATE sh1d_common)
  sh1d_enable_warnings(sh1d_hdf5_export)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(sh1d_parameters_tests tests/parameters.cpp)
  target_link_libraries(sh1d_parameters_tests PRIVATE sh1d_common)
  sh1d_enable_warnings(sh1d_parameters_tests)
  add_test(NAME parameter_parsing COMMAND sh1d_parameters_tests)

  add_executable(sh1d_numerics_tests tests/numerics.cpp)
  target_link_libraries(sh1d_numerics_tests PRIVATE sh1d_common)
  sh1d_enable_warnings(sh1d_numerics_tests)
  add_test(NAME nonlinear_direct_convolution COMMAND sh1d_numerics_tests)

  add_executable(sh1d_output_transaction_tests tests/output_transaction.cpp)
  target_link_libraries(sh1d_output_transaction_tests PRIVATE sh1d_common)
  sh1d_enable_warnings(sh1d_output_transaction_tests)
  add_test(NAME output_transaction_recovery COMMAND sh1d_output_transaction_tests)

  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME regression_cpu COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.py
      --executable $<TARGET_FILE:schrodinger_helmholtz>)
    set_tests_properties(regression_cpu PROPERTIES TIMEOUT 120)
    add_test(NAME temporal_convergence COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/convergence.py
      --executable $<TARGET_FILE:schrodinger_helmholtz>)
    set_tests_properties(temporal_convergence PROPERTIES TIMEOUT 120)
    if(TARGET sh1d_hdf5_export)
      add_test(NAME hdf5_output COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/hdf5_output.py
        --executable $<TARGET_FILE:schrodinger_helmholtz>
        --exporter $<TARGET_FILE:sh1d_hdf5_export>)
      set_tests_properties(hdf5_output PROPERTIES TIMEOUT 120)
    endif()
  endif()
endif()
