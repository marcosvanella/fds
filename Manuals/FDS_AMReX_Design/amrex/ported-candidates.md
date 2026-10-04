# Ported-status candidates (D-072 (c)) and the fast-math pin check

Owner: AMReX Integration Lead. Status: v0.1. `tools/ported.toml` is **unchanged: no kernel is set to ported**. This file lists what is ready, what each entry still needs, and the record of the CMake pin check.

## 1. Why no entry is set yet

D-072 (c) needs four things per kernel: a device run on record on real fields within the class tolerance; a review by the kernel owner; the V&V ulp gate for a libm kernel; and `tools/ci_checks.sh` (--strict) passing. State of the four:

| Condition | State |
|---|---|
| Device run on real fields, class tolerance | **Met** for the kernels in section 2 (run, cc 8.9 test-machine GPU; see `stage1-gpu-spike-plan.md` 7.9c, 7.10a, 7.13, 7.13a, 7.13b, 7.10d). |
| Review by the kernel owner | **Not on record.** The plan review of the stage-1 work approved work packages and read documents only ("nothing run"); it is not a kernel-result review, so it is not entered as `reviewed_by`. The wall-kernel code reviews in `solid/07` and `solid/08` are reviews on reading of the generated text, without a device run. Each entry below needs a named reviewer who has seen the device-run result. |
| V&V ulp gate (libm kernel `cfl_wall_max`) | **Not met in the gate's own terms.** The gate (`vv-runs/gpu_gate`) judges a libm kernel only by a `ULP` line printed by a device driver. The stage-1 harness result (device within 1 ulp of host on 32 real walls; 2 ulp maximum on a SYNTHETIC sweep of 20,000 walls) is a harness result, not the gate's `ULP` line. `ulp_gate_passed = true` stays unset. |
| `ci_checks.sh --strict` | **Fails** for a reason outside this work, see section 3: `s5gen_rad_wall_qin_zero` in `test/rad_kernels.golden` (K2-03 and CL-07, array dummy `W_SURF_INDEX` missing from the device-address list). |

An entry can be added in one step once a reviewer is on record and `ci_checks.sh` is green.

## 2. Candidates (kernel-map names; evidence is run on the test machine's GPU unless stated)

| Kernels | Evidence | Tolerance class |
|---|---|---|
| `vpred_us`, `vpred_vs`, `vpred_ws`, `vcorr_u`, `vcorr_v`, `vcorr_w`, `div2_pred`, `div2_corr`, `baro_p_rrho`, `baro_fvx`, `baro_fvy`, `baro_fvz`, `vflux_vort_tau` | 7.13: device = host serial = host 4 threads on every output array (`csmag_32`, four `dec2_obst` boxes); bit-equal to FDS where an FDS value exists (7.9c) | bitwise |
| `vflux_fvx`, `vflux_fvy`, `vflux_fvz` | 7.13b: real fields and real edge tables, host and device; bit-equal to FDS on `csmag_32` (33,792 of 33,792 faces each) and `dec2_obst` (272/512/272 faces per box) | bitwise |
| `cfl_max`, `vn_max`, `div_extrema` | 7.13a: device = serial = 4 threads on 155 of 155 outputs (+0 equals -0 rule, 7.10c); `div_extrema` checked independently | bitwise |
| `d_z_max`, `rho_d_interp`, `rho_d_maxloc_fix`, `del_rho_d_del_z`, `dp_div_heat` | 7.13a: same runs; `RHO_D` and `D_Z` inputs are SYNTHETIC, no FDS truth | bitwise (map shows `host` or `generated`: the map does not list a device run for these yet, so PM-08 would reject them until the device-run overlay lists the stage-1 run) |
| `rho_d_dzd`, `h_rho_d_dzd` | 7.13a neighbours | bitwise |
| WP1 non-wall and wall kernels (`7.9`, `7.10`, `7.10a`): 266 of 266 non-wall and 76 of 76 wall output arrays bit-equal across device, host serial and 4 threads (7 data sets) | 7.9, 7.10a | bitwise |
| `cfl_wall_max` | 7.10d: host serial = 4 threads = FDS-build value on 32 of 32 real walls (hot-obstruction variant of `dec2_obst`, scratch input); device within 1 ulp (4 of 32 walls differ); SYNTHETIC sweep device-vs-host maximum 2 ulp | libm, 2 ulp |

Not candidates: the six `gsfv_*` kernels (no device run), loops without a kernel or marker, and `wall_rho_d_dzdn` for its `NIC>1` branch (no case reaches it).

## 3. Fast-math pin (D-068): check on a scratch copy

The pin patch (`tools/patches/cmake-fastmath-pin.patch`) edits the generator worktree (`amrex/s4_mass/CMakeLists.txt`, `s4d/CMakeLists.txt`, `s5_gen/gpu/CMakeLists.txt`, and two new files under `amrex/cmake/`). That worktree belongs to the generator branch, so the patch was **not** applied there; it applies cleanly (`git apply --check`) and was checked on a scratch copy of the tree (run):

- `kernel_lint.py --strict` on the unpatched tree: BF-05 FAIL, BF-07 FAIL.
- Same on the patched scratch copy: BF-05 NOTE ("D-068 pin found: FORCE-OFF in amrex/cmake/AMReXCudaPin.cmake, failing guard in the same file"); BF-07 no finding (pass). The tool reported 1 remaining FAIL, which is the CL-07 finding of the radiation golden named in section 1.
- `k2_ci_check.py --strict` (patched or not): 2 FAIL, both K2-03 on `s5gen_rad_wall_qin_zero`, `W_SURF_INDEX`.
- `zone_sum_order.py --strict`: PASS. `port_kernel_map.py --strict`: self-check PASS.

Remaining step for BF-05 and BF-07: the owner of the generator worktree applies the patch there (`git apply` from the worktree root).
