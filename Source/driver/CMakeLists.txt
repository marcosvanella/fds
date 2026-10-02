# Source/driver/CMakeLists.txt: the USE_AMREX build of FDS-AMReX (M2a, S1).
#
# This file is include()d by the top-level CMakeLists.txt (patch 0002) only when USE_AMREX=ON, after the `fds` target has been fully
# configured. It turns that target into the AMReX driver executable: the FDS Fortran sources are compiled with -DWITH_AMREX (main.f90
# then provides the callable subroutine FDS_SETUP instead of PROGRAM FDS, patch 0001) and the C++ main from this directory is added.
# With USE_AMREX=OFF this file is never read, so the FDS build is unchanged (IR-006).
#
# Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp (S2); (b) only uniform Cartesian metrics are used.

enable_language(CXX)
enable_language(C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# AMReX installation (gfortran 14 + OpenMPI 5, FFT and OMP enabled). Override with -DAMReX_ROOT=<prefix> or -DCMAKE_PREFIX_PATH.
set(AMREX_INSTALL_DIR "/workspace/amrex-install" CACHE PATH "AMReX install prefix used when USE_AMREX=ON")
find_package(AMReX 26.09 REQUIRED CONFIG PATHS ${AMREX_INSTALL_DIR} NO_DEFAULT_PATH)
find_package(MPI REQUIRED COMPONENTS CXX)

# The FDS options above are Fortran-only (-std=f2018, -fall-intrinsics, ...): restrict every existing option to the Fortran compiler.
get_target_property(_fds_opts fds COMPILE_OPTIONS)
set_property(TARGET fds PROPERTY COMPILE_OPTIONS "")
foreach(_o IN LISTS _fds_opts)
    target_compile_options(fds PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${_o}>")
endforeach()
# The version/date strings are Fortran preprocessor macros with embedded quotes: keep them out of the C++ command lines.
get_target_property(_fds_defs fds COMPILE_DEFINITIONS)
set_property(TARGET fds PROPERTY COMPILE_DEFINITIONS "")
foreach(_d IN LISTS _fds_defs)
    target_compile_definitions(fds PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${_d}>")
endforeach()
target_compile_definitions(fds PRIVATE WITH_AMREX)

target_sources(fds PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FdsAmr.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FdsSetup.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Fields.cpp
    ${CMAKE_CURRENT_LIST_DIR}/GhostExchange.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SideData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExactSum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LevelRegistry.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TimeLoop.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/CommonLayer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/PressureIface.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/FFTBackend.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/MLMGBackend.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tests/selftest_fds.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tests/kernelcheck.cpp
    ${CMAKE_CURRENT_LIST_DIR}/fds_mesh_query.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_box_shim.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_kernels.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_amrex_hooks.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_ghost_bc.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_density_split.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_clip_gather.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_step.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_alias.c
)
set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/fds_box_shim.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_mesh_query.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_kernels.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_ghost_bc.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_amrex_hooks.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_density_split.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_clip_gather.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_step.f90
    PROPERTIES COMPILE_OPTIONS "-ffree-line-length-none")
target_include_directories(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport)
target_compile_options(fds PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra>")

# C++ is the link language (C++ main; the Fortran runtime is pulled in by the gfortran libraries listed by CMake).
set_target_properties(fds PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME fds_amr)
target_link_libraries(fds PRIVATE AMReX::amrex MPI::MPI_CXX)
if (USE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_link_libraries(fds PRIVATE OpenMP::OpenMP_CXX)
endif()
message(STATUS "USE_AMREX=ON: AMReX ${AMReX_VERSION} from ${AMREX_INSTALL_DIR}; executable fds_amr")

# Driver unit tests (no FDS objects, AMReX only): index maps, ghost widths, face fill, side data, tile race. Run by tests/run_driver_tests.sh
# (1 and 4 ranks, threads 1). The FDS-linked checks are the `--selftest` mode of fds_amr (tests/selftest_fds.cpp).
add_executable(driver_unit_tests
    ${CMAKE_CURRENT_LIST_DIR}/tests/test_units.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FdsAmr.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Fields.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SideData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExactSum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/LevelRegistry.cpp)
target_include_directories(driver_unit_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport)
target_compile_options(driver_unit_tests PRIVATE -Wall -Wextra)
target_link_libraries(driver_unit_tests PRIVATE AMReX::amrex MPI::MPI_CXX)
if (USE_OPENMP)
    target_link_libraries(driver_unit_tests PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(driver_unit_tests PROPERTIES LINKER_LANGUAGE CXX)
