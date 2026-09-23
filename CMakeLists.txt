cmake_minimum_required(VERSION 3.20)
project(SchrodingerHelmholtz1D VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SH1D_OPENMP "Enable OpenMP and threaded FFTW" ON)

function(sh1d_enable_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

find_path(FFTW_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW_LIBRARY NAMES fftw3 REQUIRED)
find_library(FFTW_THREADS_LIBRARY NAMES fftw3_threads)

add_library(sh1d_common
    src/parameters.cpp
    src/fftw_utils.cpp
    src/spectral.cpp
    src/nonlinear.cpp
    src/output.cpp
    src/solver.cpp)
target_include_directories(sh1d_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${FFTW_INCLUDE_DIR})
target_link_libraries(sh1d_common PUBLIC ${FFTW_LIBRARY})
target_compile_features(sh1d_common PUBLIC cxx_std_20)

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

add_executable(schrodinger_helmholtz src/main.cpp)
target_link_libraries(schrodinger_helmholtz PRIVATE sh1d_common)
sh1d_enable_warnings(schrodinger_helmholtz)

include(CTest)
if(BUILD_TESTING)
  add_executable(sh1d_numerics_tests tests/numerics.cpp)
  target_link_libraries(sh1d_numerics_tests PRIVATE sh1d_common)
  sh1d_enable_warnings(sh1d_numerics_tests)
  add_test(NAME nonlinear_and_parameters COMMAND sh1d_numerics_tests)

  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME regression_cpu COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.py
      --executable $<TARGET_FILE:schrodinger_helmholtz>)
    set_tests_properties(regression_cpu PROPERTIES TIMEOUT 120)
  endif()
endif()
