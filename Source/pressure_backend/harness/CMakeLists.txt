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

# HYPRE backend (HypreBackend.cpp): compiled, and the HYPRE references of the other sources enabled, only with PB_WITH_HYPRE (README.md,
# "Build options"). Default ON when the installed AMReX was built with HYPRE. This harness tests the HYPRE backend, so it needs ON; the
# build WITHOUT it (what the driver builds) is checked by the pb_nohypre_* targets below.
option(PB_WITH_HYPRE "Build the HYPRE backend (needs AMReX built with HYPRE)" ${AMReX_HYPRE_FOUND})
if(NOT PB_WITH_HYPRE)
  message(FATAL_ERROR "the harness tests the HYPRE backend: configure with an AMReX built with HYPRE (e.g. CMAKE_PREFIX_PATH=/workspace/amrex-install-hypre232)")
endif()
if(NOT AMReX_HYPRE_FOUND)
  message(FATAL_ERROR "PB_WITH_HYPRE=ON but the installed AMReX has no HYPRE")
endif()

set(PB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_library(pressure_backend STATIC
  ${PB_DIR}/CommonLayer.cpp ${PB_DIR}/PressureIface.cpp ${PB_DIR}/FFTBackend.cpp ${PB_DIR}/MLMGBackend.cpp
  ${PB_DIR}/CompositeSolve.cpp ${PB_DIR}/HypreBackend.cpp ${PB_DIR}/FdsPressureLoops.cpp)
# host translations of FDS pressure-area loops: bitwise against FDS, so no FMA contraction in that file
set_source_files_properties(${PB_DIR}/FdsPressureLoops.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
target_include_directories(pressure_backend PUBLIC ${PB_DIR})
target_compile_definitions(pressure_backend PUBLIC PB_WITH_HYPRE)   # every pressure_backend translation unit must see the same setting
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

# Without PB_WITH_HYPRE: the driver's source list (Source/driver/CMakeLists.txt) plus HypreBackend.cpp, which the driver lists when AMReX has HYPRE
# but does not define the macro, so it must compile to an empty object. No PB_WITH_HYPRE here on purpose (not linked to the pressure_backend target).
add_library(pb_nohypre_obj OBJECT
  ${PB_DIR}/CommonLayer.cpp ${PB_DIR}/PressureIface.cpp ${PB_DIR}/FFTBackend.cpp ${PB_DIR}/MLMGBackend.cpp ${PB_DIR}/CompositeSolve.cpp
  ${PB_DIR}/HypreBackend.cpp)
target_include_directories(pb_nohypre_obj PUBLIC ${PB_DIR})
target_link_libraries(pb_nohypre_obj PUBLIC AMReX::amrex)
add_executable(pb_nohypre_check nohypre_check.cpp $<TARGET_OBJECTS:pb_nohypre_obj>)
target_include_directories(pb_nohypre_check PRIVATE ${PB_DIR})
target_link_libraries(pb_nohypre_check PRIVATE AMReX::amrex)
add_test(NAME pb_nohypre_check COMMAND pb_nohypre_check)
add_test(NAME pb_nohypre_symbols COMMAND ${CMAKE_COMMAND} -DNM=${CMAKE_NM} "-DOBJS=$<TARGET_OBJECTS:pb_nohypre_obj>" -P ${PB_DIR}/tests/check_no_hypre_symbols.cmake)
