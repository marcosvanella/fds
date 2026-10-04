# R4 part 2: results of the face path, mass-weighted restriction, dynamic blob and GPU tagging checks

## Restriction of species (steering item)
- Average-down of species uses rho and rho*Z over covered cells and derives Z = (rho*Z)/rho. `LevelOps::average_down_species` is used by
  `average_down_registry` and by the coarse-fine ghost hook (`DriverAdapter.cpp`); no other restriction path touches mass fractions.
- `tests/test_species_avgdown.cpp`: bare level, registry and hook; negative control with a linear Z average changes species mass (1e-1 relative on the
  test data) and the test fails against it; RHO/ZZ/ZZS agree bitwise on covered cells with the driver's `RegistryTransfer::hierarchy_done`.

## Face path (D-060)
- `prolong_faces_level`, `average_down_faces` (all components), `max_divergence_error_local` in `FaceTransfer`; `tests/test_facetransfer.cpp`.
- All-new fine cells have the parent's divergence to about 1e-14. Next to retained fine faces (seams) max|div u - D| is of the size of the velocity
  gradient times the cell width; reported by `test_blob_registry`, not asserted. The projection of the pressure solve removes it each step.

## Dynamic blob end to end (`tests/test_blob_registry.cpp`)
- Real RegridAmrCore, LevelRegistry and RegistryTransfer; mock conservative upwind transport of rho*Z with a prescribed constant-in-time velocity;
  2-D and 3-D, ratios 2 and 4, up to 3 levels, regrid every 4 steps.
- Composite mass and species change per regrid below double resolution; no negative rho*Z, no clips; same hierarchy and data hash at 1, 2 and 4 ranks
  (`run_regrid_rank_check.sh`); AMR closer to the uniform-fine run than uniform-coarse; cells beyond the dependence cone equal the uniform-fine run bitwise.
- Negative controls: no interface flux overwrite drifts mass; linear Z restriction breaks mass; no tag buffer leaves feature cells off the finest level.
- `&AMR_REGION` must be given: the default region is empty and all tags are discarded.

## GPU tagging kernels
- `tests/tag_kernel_check.F90` built twice (host default, `RT_OFFLOAD`), each against an independent cell-by-cell reference, outputs compared bitwise
  (`regrid_transport_tagkernel_host_vs_offload`). Without an accelerator the offload build runs its target regions on the host.
- Device run: `nvfortran -fast -Mpreprocess -mp=gpu -gpu=mem:managed -DRT_OFFLOAD rt_tag_kernels.F90 tag_kernel_check.F90`; 38 cases identical to the host
  build (kernels launched on the GPU). CMake: `-DRT_OFFLOAD_FLAGS="-mp=gpu;-gpu=mem:managed"` with nvfortran.
- Using the offload kernels inside the C++ library needs an AMReX build with GPU support (device-resident TagBox/FArrayBox); the installed AMReX is CPU only.
