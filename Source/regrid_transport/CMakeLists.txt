# Role 3 library target and tests (regrid_transport).
#
# Stand-alone use (no FDS, no AMReX needed for the parser tests):
#   cmake -S Source/regrid_transport -B build-rt && cmake --build build-rt && ctest --test-dir build-rt
# Inside the FDS build the top-level file include()s this directory later (patch list in notes/); not yet included there.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.16)
    project(fds_regrid_transport_standalone LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    include(CTest)
endif()

# AMReX-independent part: &AMR parser (R0).
add_library(fds_regrid_transport_core STATIC AmrInput.cpp)
target_include_directories(fds_regrid_transport_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(fds_regrid_transport_core PRIVATE -Wall -Wextra)

# Interface target kept from the stub; consumers link this.
add_library(fds_regrid_transport INTERFACE)
target_include_directories(fds_regrid_transport INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fds_regrid_transport INTERFACE fds_regrid_transport_core)

if (BUILD_TESTING)
    add_executable(test_amr_input tests/test_amr_input.cpp)
    target_link_libraries(test_amr_input PRIVATE fds_regrid_transport_core)
    target_compile_options(test_amr_input PRIVATE -Wall -Wextra)
    add_test(NAME regrid_transport_amr_input COMMAND test_amr_input)
endif()
