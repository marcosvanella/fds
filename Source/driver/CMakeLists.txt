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
    ${CMAKE_CURRENT_LIST_DIR}/CompositeOps.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TwoLevelRun.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/PressureBackendSolver.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/CommonLayer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/PressureIface.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/FFTBackend.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/MLMGBackend.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tests/selftest_fds.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tests/kernelcheck.cpp
    ${CMAKE_CURRENT_LIST_DIR}/tests/level_bind_check.cpp
    ${CMAKE_CURRENT_LIST_DIR}/fds_mesh_query.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_box_shim.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_kernels.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_box_obj.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_flux_hooks.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_wall_seam.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_amrex_hooks.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_ghost_bc.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_density_split.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_clip_gather.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_step.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_alias.c
)
set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/fds_box_shim.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_box_obj.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_flux_hooks.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_wall_seam.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_mesh_query.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_kernels.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_ghost_bc.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_amrex_hooks.f90
    ${CMAKE_CURRENT_LIST_DIR}/fds_density_split.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_clip_gather.f90 ${CMAKE_CURRENT_LIST_DIR}/fds_step.f90
    PROPERTIES COMPILE_OPTIONS "$<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>")
# Role 2's composite (multi-level) solve lives in its own file once it exists (PressureIface.cpp then references pb::solve_composite); listed only when the file is present.
if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/CompositeSolve.cpp)
    target_sources(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/CompositeSolve.cpp)
endif()
# Role 2's HYPRE backend (HypreBackend.cpp) needs an AMReX built with HYPRE; without it the link stub provides the symbols (selecting HYPRE then reports "not built").
if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/HypreBackend.cpp)
    if (AMReX_HYPRE_FOUND)
        target_sources(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend/HypreBackend.cpp)
    else()
        target_sources(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR}/HypreStub.cpp)
    endif()
endif()
target_include_directories(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport)
# Role 3 driver-level end-to-end modes (regrid_transport/DriverSources.cmake): `fds_amr <case> --rt-e2e <test>`
if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/DriverSources.cmake)
    include(${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/DriverSources.cmake)
    # Role 3's tag kernel has lines longer than gfortran's default 132 columns
    set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/rt_tag_kernels.F90 PROPERTIES COMPILE_OPTIONS "$<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>")
endif()
target_compile_options(fds PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra>")

# D-056 option B: FDS_FINE_B is defined for the driver Fortran when the tree carries patch 0007 (Source/mesh.f90 declares FINE_LEVEL_TYPE); fds_box_obj.f90 then resolves fine
# mesh numbers to FINE_LEVEL(L)%BOX(IB). Without the patch the driver is unchanged.
get_filename_component(_fds_mesh_f90 ${CMAKE_CURRENT_LIST_DIR}/../mesh.f90 ABSOLUTE)   # lexical '..': a tree whose driver directory is a symlink must read its own mesh.f90
file(READ ${_fds_mesh_f90} _fds_mesh_src)
if (_fds_mesh_src MATCHES "FINE_LEVEL_TYPE")
    target_compile_definitions(fds PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:FDS_FINE_B>")
    message(STATUS "Patch 0007 found in Source/mesh.f90: FDS_FINE_B defined")
    # Fine-level mesh objects (BUILD_FINE_BOX, fds_fine_level_create/destroy) for TimeLoop::bind_level; FDS_HAVE_FINE_LEVEL tells the C++ side that they exist.
    target_sources(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR}/fds_fine_level.f90)
    set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/fds_fine_level.f90 PROPERTIES COMPILE_OPTIONS "$<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>")
    target_compile_definitions(fds PRIVATE FDS_HAVE_FINE_LEVEL)
endif()

# D-056 option B DRAFT (patch 0007, mesh.f90): needs the patch applied to the tree; OFF by default and not used for physics. Adds draft/fds_fine_mesh_b.f90 and the
# `--fine-b-selftest` / `--fine-b-abort N` modes of fds_amr.
option(FDS_AMR_FINE_B_DRAFT "Build the fine-level mesh object draft (needs patch 0007)" OFF)
if (FDS_AMR_FINE_B_DRAFT)
    target_sources(fds PRIVATE ${CMAKE_CURRENT_LIST_DIR}/draft/fds_fine_mesh_b.f90 ${CMAKE_CURRENT_LIST_DIR}/draft/fds_fine_box_b.f90)
    set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/draft/fds_fine_mesh_b.f90 ${CMAKE_CURRENT_LIST_DIR}/draft/fds_fine_box_b.f90 PROPERTIES COMPILE_OPTIONS "$<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>")
    target_compile_definitions(fds PRIVATE FDS_FINE_B_DRAFT)
endif()

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
    ${CMAKE_CURRENT_LIST_DIR}/LevelRegistry.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RegistryTransfer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/CellTransfer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/FaceTransfer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/LevelOps.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/FluxOverrideOps.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport/FluxStageRunner.cpp)
target_include_directories(driver_unit_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../regrid_transport ${CMAKE_CURRENT_LIST_DIR}/../pressure_backend)
target_compile_options(driver_unit_tests PRIVATE -Wall -Wextra)
target_link_libraries(driver_unit_tests PRIVATE AMReX::amrex MPI::MPI_CXX)
if (USE_OPENMP)
    target_link_libraries(driver_unit_tests PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(driver_unit_tests PROPERTIES LINKER_LANGUAGE CXX)
