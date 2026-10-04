# Standalone build of the pressure backend harness. No FDS sources; not part of any top-level build.
#   cmake -S Source/pressure_backend/harness -B <build> -DCMAKE_PREFIX_PATH=/workspace/amrex-install-hypre232 \
#         -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release
cmake_minimum_required(VERSION 3.25)
# Fortran is enabled only because the installed AMReX config requires it (AMReX_FORTRAN=ON); no Fortran sources.
project(PressureBackendHarness LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(AMReX REQUIRED)
foreach(_c FFT LSOLVERS MPI)
  if(NOT AMReX_${_c})
    message(FATAL_ERROR "Installed AMReX lacks AMReX_${_c}; build a new AMReX install (see README.md)")
  endif()
endforeach()

set(PB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_library(pressure_backend STATIC
  ${PB_DIR}/CommonLayer.cpp ${PB_DIR}/PressureIface.cpp ${PB_DIR}/FFTBackend.cpp ${PB_DIR}/MLMGBackend.cpp
  ${PB_DIR}/CompositeSolve.cpp ${PB_DIR}/HypreBackend.cpp ${PB_DIR}/FdsPressureLoops.cpp)
# host translations of FDS pressure-area loops: bitwise against FDS, so no FMA contraction in that file
set_source_files_properties(${PB_DIR}/FdsPressureLoops.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
target_include_directories(pressure_backend PUBLIC ${PB_DIR})
target_link_libraries(pressure_backend PUBLIC AMReX::amrex)

add_executable(pb_harness main.cpp composite_modes.cpp m2_modes.cpp)
target_link_libraries(pb_harness PRIVATE pressure_backend)

# M4: FDS-derived frozen cases (own main, reads FDS dumps; frozen/fds_*/README.md)
add_executable(pb_fds_frozen m4_modes.cpp)
target_link_libraries(pb_fds_frozen PRIVATE pressure_backend)

enable_testing()
add_subdirectory(${PB_DIR}/tests ${CMAKE_CURRENT_BINARY_DIR}/tests)
include(${PB_DIR}/tests/m4_fds_frozen.cmake)   # pb_fds_frozen_* tests
include(${PB_DIR}/tests/m5_fds_loops.cmake)   # pb_fds_loops, mutants and pb_fdsloops_* tests (host translations of pres.f90 loops)
