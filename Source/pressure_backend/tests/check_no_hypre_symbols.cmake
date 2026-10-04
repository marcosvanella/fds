# ctest helper: fail if any object of the pressure backend compiled without PB_WITH_HYPRE names a HYPRE symbol (defined or undefined).
# Two kinds of names are expected and ignored: the accessor PressureWorkspace::hypre_built() (part of the public API, always false
# without HYPRE) and amrex::Hypre* destructor references that AMReX's own headers emit when AMReX itself was built with HYPRE
# (an AMReX built without HYPRE has neither; the scratch link test against such an AMReX shows zero).
# Usage: cmake -DNM=<nm> "-DOBJS=<a.o;b.o;...>" -P check_no_hypre_symbols.cmake
if(NOT NM OR NOT OBJS)
  message(FATAL_ERROR "usage: -DNM=<nm> -DOBJS=<objects>")
endif()
list(REMOVE_ITEM OBJS "")
set(bad "")
foreach(o ${OBJS})
  execute_process(COMMAND sh -c "'${NM}' -C '${o}' | grep -i hypre | grep -v 'amrex::Hypre' | grep -v 'PressureWorkspace::hypre_built'"
                  OUTPUT_VARIABLE hits)
  if(hits)
    string(APPEND bad "${o}:\n${hits}\n")
  endif()
endforeach()
if(bad)
  message(FATAL_ERROR "HYPRE symbols in a build without PB_WITH_HYPRE:\n${bad}")
endif()
list(LENGTH OBJS n)
message(STATUS "no HYPRE symbol in ${n} objects")
