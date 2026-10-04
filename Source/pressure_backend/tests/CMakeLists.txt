# CTest registrations (included from harness/CMakeLists.txt). All tests run the harness through pb_test.py.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(PB_MPIEXEC NAMES mpiexec mpirun REQUIRED)
set(PB_ENV "OMP_NUM_THREADS=1;OMP_DYNAMIC=false")

function(pb_add_test name)
  add_test(NAME ${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/pb_test.py
            --harness $<TARGET_FILE:pb_harness> --mpiexec ${PB_MPIEXEC}
            --work ${CMAKE_CURRENT_BINARY_DIR}/work_${name} ${ARGN})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${PB_ENV}" TIMEOUT 600)
endfunction()

pb_add_test(pb_selector      selector)
pb_add_test(pb_exactsum      exactsum)
foreach(bc neumann periodic)
  pb_add_test(pb_fft_vs_mlmg_${bc}_64     fftmlmg --bc ${bc} --n "64 64 64")
  pb_add_test(pb_fft_vs_mlmg_${bc}_34x18x32 fftmlmg --bc ${bc} --n "34 18 32")
  pb_add_test(pb_frozen_${bc}_34x18x32    frozen  --bc ${bc} --n "34 18 32")
  pb_add_test(pb_decomp_${bc}_64          decomp  --bc ${bc} --n "64 64 64")
  pb_add_test(pb_repeat_${bc}_64          repeat  --bc ${bc} --n "64 64 64")
endforeach()
pb_add_test(pb_singular_hook singular)
# Mean removal and gauge versus the FDS ULMAT recipe (independent numpy reference), per-cell volume zero modes.
foreach(bc neumann periodic)
  pb_add_test(pb_ulmat_equiv_${bc}   ulmat      --bc ${bc} --n "8 6 5")
  pb_add_test(pb_ulmat_gauge_${bc}   ulmatgauge --bc ${bc} --n "8 6 5")
endforeach()
pb_add_test(pb_meankind_stretched   meankind   --bc neumann --n "12 10 9")

# Composite (multi-level) solves: manufactured solutions on refined-patch hierarchies; see frozen/composite-notes.md.
pb_add_test(pb_comp_selector   compsel)
pb_add_test(pb_comp_workspace  compws)
foreach(bc neumann periodic dirichlet)
  pb_add_test(pb_comp_conv_${bc}_r2     compconv --bc ${bc} --ratio 2 --ns "32 64")
  pb_add_test(pb_comp_conv_${bc}_r4     compconv --bc ${bc} --ratio 4 --ns "16 32")
  pb_add_test(pb_comp_full_${bc}_r2     compfull --bc ${bc} --ratio 2)
  pb_add_test(pb_comp_full_${bc}_r4     compfull --bc ${bc} --ratio 4)
endforeach()
foreach(bc neumann periodic)
  pb_add_test(pb_comp_decomp_${bc}      compdecomp --bc ${bc})
endforeach()
foreach(bc neumann periodic dirichlet)
  pb_add_test(pb_comp_three_levels_${bc} comp3 --bc ${bc})
  pb_add_test(pb_comp_gradient_${bc}     compgrad --bc ${bc})
endforeach()
foreach(bc neumann periodic)
  pb_add_test(pb_comp_plane2d_${bc}      compgrad --bc ${bc} --plane 1)
  pb_add_test(pb_comp_plane2d_conv_${bc} compconv --bc ${bc} --plane 1 --ns "32 64")
endforeach()
foreach(bc neumann periodic dirichlet)
  pb_add_test(pb_comp_shapes_${bc}      compshape --bc ${bc})
endforeach()
pb_add_test(pb_comp_mixed_bc compmixed)
pb_add_test(pb_comp_ns2d_16  compns2d)

# D-067: mean removal (Volume default, ScaledArithmetic parity switch) and the rho/KRES gauge on a hierarchy.
foreach(bc neumann periodic dirichlet)
  pb_add_test(pb_comp_gauge_${bc}_r2   compgauge --bc ${bc} --ratio 2)
endforeach()
pb_add_test(pb_comp_gauge_neumann_r4      compgauge --bc neumann --ratio 4)
pb_add_test(pb_comp_gauge_plane2d_neumann compgauge --bc neumann --ratio 2 --plane 1)
pb_add_test(pb_comp_gauge_plane2d_periodic compgauge --bc periodic --ratio 2 --plane 1)
foreach(bc neumann periodic)
  pb_add_test(pb_comp_gauge_decomp_${bc} compgaugedecomp --bc ${bc})
endforeach()
foreach(bc neumann periodic)
  pb_add_test(pb_meankind_uniform_${bc} meankind_uniform --bc ${bc} --n "8 6 5")
endforeach()

# FR-039 solve-trigger points (full checks versus the cheap path) and FFT plan caching in PressureWorkspace.
foreach(bc "NN,NN,NN" "PP,PP,PP")
  string(REPLACE "," "" bcn ${bc})
  pb_add_test(pb_trigger_single_${bcn} trigger1 --bc ${bc} --n "16 16 16")
endforeach()
pb_add_test(pb_comp_trigger comptrigger)
pb_add_test(pb_fftcache_neumann fftcache --bc "NN,NN,NN" --n "32 32 32")
pb_add_test(pb_fftcache_periodic fftcache --bc "PP,PP,PP" --n "32 32 32")

# Mixed open/closed faces (D-057 / M2) on a two-level hierarchy: convergence against a manufactured solution
# (N at the low and D at the high face: cos((k+1/2) pi x); D low, N high: sin((k+1/2) pi x)).
foreach(faces "ND,NN,NN" "DD,NN,ND" "NN,PP,DN" "DN,NN,DD")
  string(REPLACE "," "" fn ${faces})
  pb_add_test(pb_comp_faces_${fn} compconv --faces ${faces} --ns "32 64")
endforeach()
pb_add_test(pb_comp_faces_plane2d_DDNNDD compconv --faces "DD,DD,DD" --plane 1 --ns "32 64")
pb_add_test(pb_comp_faces_plane2d_NDDDNN compconv --faces "ND,DD,NN" --plane 1 --ns "32 64")
pb_add_test(pb_mixed_faces_sweep mixedfaces)

# Inhomogeneous boundary data folded into the right-hand side (M2 item d).
pb_add_test(pb_bcdata_exact bcdata_exact)
pb_add_test(pb_bcdata_mms   bcdata_mms)

# D-057: boundary-type mapping for one-cell directions and the sweep strings of the driver.
pb_add_test(pressure_bc_map pressure_bc_map)

# HYPRE assembled-matrix backend (frozen/hypre-notes.md): explicit request only, same discrete operator as MLMG.
# (a) FFT vs MLMG vs HYPRE agree within eps_H on the single-level cases (Neumann, periodic, Dirichlet, non-cubic, one-cell y).
foreach(bc neumann periodic dirichlet)
  pb_add_test(pb_hypre_single_${bc}_34x18x32 hypsingle --bc ${bc} --n "34 18 32")
  pb_add_test(pb_hypre_single_${bc}_64       hypsingle --bc ${bc} --n "64 64 64")
endforeach()
foreach(bc neumann periodic)
  pb_add_test(pb_hypre_single_${bc}_plane2d  hypsingle --bc ${bc} --n "16 1 16")
endforeach()
pb_add_test(pb_hypre_mixed_faces_sweep mixedfaces --be hypre)
pb_add_test(pb_hypre_bcdata_exact      bcdata_exact --be hypre)
# (c) decomposition independence (1/2/4 ranks, two box sizes) and run-to-run bitwise repeatability.
foreach(bc neumann periodic)
  pb_add_test(pb_hypre_decomp_${bc}_64 decomp --be hypre --bc ${bc} --n "64 64 64")
  pb_add_test(pb_hypre_repeat_${bc}_64 repeat --be hypre --bc ${bc} --n "64 64 64")
endforeach()
# (d) CI gate: the same frozen solves through MLMG and HYPRE must agree within eps_H (with a negative control).
pb_add_test(pb_hypre_ci_vs_mlmg hypci --n "34 18 32")
# (b) composite hierarchies: HYPRE vs MLMG, independent true residual, repeatability, rank independence.
pb_add_test(pb_hypre_comp_vs_mlmg  hypcomp)
pb_add_test(pb_hypre_comp_ns2d_16  compns2d --be hypre)
# (e) selector (NotBuilt for masks, covered cells, cylindrical, stretched, variable coefficients) is in pb_selector;
#     the composite selector checks are repeated with HYPRE as the requested backend.
pb_add_test(pb_hypre_comp_selector compsel --be hypre)
pb_add_test(pb_hypre_comp_workspace compws --be hypre)
foreach(bc "NN,NN,NN" "PP,PP,PP")
  string(REPLACE "," "" bcn ${bc})
  pb_add_test(pb_hypre_cache_${bcn} hypcache --bc ${bc} --n "32 32 32")
endforeach()

# (f) pin-aware true-residual check of singular problems (frozen/hypre-notes.md, "Residual check"): no spurious warning at 48^3..96^3,
#     negative controls still warn, MLMG/FFT checked by the same code.
foreach(bc neumann periodic)
  pb_add_test(pb_hypre_resid_${bc}   hypresid --bc ${bc} --ns "48 64 96")
  pb_add_test(pb_resid_check_${bc}   rescheck --bc ${bc} --n "48 48 48")
endforeach()
pb_add_test(pb_hypre_resid_comp hypresidcomp)
# (g) size-scaled residual limit (FR-031, A-67): limit = max(residual_tol, 10 u ||A|| ||H||/||b||) for singular and non-singular problems,
#     every backend, floor scaling with N, negative controls still warn (the synthetic formula check is in rescheck).
pb_add_test(pb_resid_limit reslimit --ns "48 96")

# (h) masked single-level MLMG branch (A-68, frozen/masked-notes.md): random solid cells, a solid slab (two components), Known cells, Dirichlet
#     faces, periodic wrap; labels/ncells/singular/removed_mean/solution/gauge against an independent dense numpy reference on 1/2/3 ranks;
#     negative control warns; FFT with a mask is NotBuilt.
pb_add_test(pb_masked masked)
