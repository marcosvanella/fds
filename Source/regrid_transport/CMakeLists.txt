# Role 3 library target. Empty until R0 adds sources; not yet included by the top-level build.
add_library(fds_regrid_transport INTERFACE)
target_include_directories(fds_regrid_transport INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
