# Ported-status candidates (D-072 (c)) and the fast-math pin check

Owner: AMReX Integration Lead. Status: v0.2. This file records which kernels `tools/ported.toml` sets to ported, what is blocked, and the CMake pin check.

## 1. State of `tools/ported.toml`

**22 kernels are set** (6 `[[ported]]` entries), aligned with the V&V Lead's review (`vv/ported-review.md`). `reviewed_by` says that the V&V Lead reviewed the device runs and the ulp gate; **the Architect's co-signature of the list is pending** and is not claimed. `set_by` is the Integration Lead. `port_kernel_map.py --strict` accepts all entries (self-check PASS, no PM-08 finding; the map shows 22 kernels as "ported").

| Entry | Kernels | Class | Device evidence as reviewed |
|---|---|---|---|
| 1 | `vpred_us/vs/ws`, `vcorr_u/v/w`, `div2_pred`, `div2_corr`, `vflux_vort_tau` (9) | bitwise | `csmag_32` only (real fields); `dec2_obst` is host only for these. `vpred_vs`, `vcorr_v`, `div2_corr`, `vflux_vort_tau` have no FDS truth |
| 2 | `baro_p_rrho`, `baro_fvx/fvy/fvz` (4) | bitwise | all zero (vacuous) on `csmag_32`; device evidence is the four `dec2_obst` boxes; ghost and solid inputs not physical, no FDS truth |
| 3 | `vflux_fvx/fvy/fvz` | bitwise | real fields and real edge tables; device ran `csmag_32` and `dec2_obst` box 0 only (four boxes is the host result) |
| 4 | `cfl_max`, `vn_max`, `div_extrema` | bitwise | real fields (`csmag_32`, four `dec2_obst` boxes) with `CFL_VELOCITY_NORM` 0 only for `cfl_max`; `div_extrema` only Cartesian with stored divergence off (`CARTVELDIV` all zero); other branches only by the SYNTHETIC rounds 4 to 7 hash run |
| 5 | `rho_d_dzd`, `h_rho_d_dzd` | bitwise | device coverage only: ZZ and TMP real, `RHO_D` and `D_Z` synthetic, no FDS truth; host gate PASS |
| 6 | `cfl_wall_max`, `ulp_gate_passed = true` | libm, 2 ulp | gate `ULP` lines PASS: real 32 walls max 1 ulp (only 19 distinct real values), SYNTHETIC 20,000 walls max 2 ulp; the 2 ulp bound is measured for this toolchain, not a guarantee |


**Held out, not in `ported.toml`:**
- `d_z_max`, `del_rho_d_del_z`: the review accepts them as device coverage (synthetic `RHO_D` and `D_Z`) on the stated condition that the gate's host tier shows a result for both. The latest quick run has no result line for either, so the condition is unmet and they are held out until it is.
- `rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat`: the review asks for (1) bitwise host tests in the gate (they are NO-TEST; the map shows `generated`), and (2) a real-field run: dump `RHO_D` and the `D_Z` table, rerun `WP3_DIFF` on the test machine. Cases needed: for `rho_d_interp` the real table with real `TMP`, including the end-of-table clamp; for `rho_d_maxloc_fix` the branch with more than two species (the maximum run has `NS = 2`) and a tie in `MAXV`; for `dp_div_heat` a non-synthetic `H_RHO_D_DZD*`.
- WP1 non-wall and wall kernels (266 and 76 arrays): held until the kernels are named with their truth status and the run is repeated with the device outputs kept (the earlier device outputs were removed, so the reviewer could not re-compare them).

**Open conditions of D-072 (c):**
- `tools/ci_checks.sh --strict` must be green as ruled. It is not: `k2_ci_check` (K2-03) and `kernel_lint` (CL-07) fail on `s5gen_rad_wall_qin_zero` (`W_SURF_INDEX` missing from the device-address list). Until that is fixed the entries stand for the map only.
- The ruling asks for the kernel owner's review. The V&V Lead's review is the independent review of the device results; the Architect has to accept it in place of the owner's review (or the owners' reviews are added).

## 2. Evidence (kernel-map names; run on the cc 8.9 test-machine GPU unless stated)

| Kernels | Evidence |
|---|---|
| the 9 kernels of entry 1 | 7.13: device = host serial = host 4 threads on all 54 output files of the `csmag_32` WP3 stages (gate-tool recheck: 0 failing elements); bit-equal to FDS where a value exists (7.9c) |
| the 4 baroclinic kernels | 7.13: `dec2_obst`, four boxes, device = serial = 4 threads on 8 arrays per box |
| `vflux_fvx/fvy/fvz` | 7.13b: `csmag_32` 33,792 of 33,792 faces per direction and `dec2_obst` box 0 (272 / 512 / 272), device = serial = 4 threads = FDS; negative control with neutralised edge tables fails identically |
| `cfl_max`, `vn_max`, `div_extrema` | 7.13a: 155 of 155 outputs; scalars and locations equal in five cases; `div_extrema` checked independently |
| `rho_d_dzd`, `h_rho_d_dzd` | 7.13a neighbours; `RHO_D` synthetic |
| `cfl_wall_max` | 7.10d and the gate `ULP` lines (`vv-runs/gpu_gate/device_logs/stage1_cflw`): host serial = 4 threads = FDS-build value on 32 of 32 real walls (hot-obstruction variant of `dec2_obst`, scratch input); device within 1 ulp (4 of 32 differ); SYNTHETIC sweep max 2 ulp |

Not set: the six `gsfv_*` kernels (no device run), loops without a kernel or marker, and the held-out kernels above.

## 3. Fast-math pin (D-068): check on a scratch copy

The pin patch (`tools/patches/cmake-fastmath-pin.patch`) edits the generator worktree (`amrex/s4_mass/CMakeLists.txt`, `s4d/CMakeLists.txt`, `s5_gen/gpu/CMakeLists.txt`, and two new files under `amrex/cmake/`). That worktree belongs to the generator branch, so the patch was **not** applied there; it applies cleanly (`git apply --check`) and was checked on a scratch copy of the tree (run):

- `kernel_lint.py --strict` on the unpatched tree: BF-05 FAIL, BF-07 FAIL.
- Same on the patched scratch copy: BF-05 NOTE ("D-068 pin found: FORCE-OFF in amrex/cmake/AMReXCudaPin.cmake, failing guard in the same file"); BF-07 no finding (pass). The tool reported 1 remaining FAIL, which is the CL-07 finding of the radiation golden named in section 1.
- `k2_ci_check.py --strict` (patched or not): 2 FAIL, both K2-03 on `s5gen_rad_wall_qin_zero`, `W_SURF_INDEX`.
- `zone_sum_order.py --strict`: PASS. `port_kernel_map.py --strict`: self-check PASS.

Remaining step for BF-05 and BF-07: the owner of the generator worktree applies the patch there (`git apply` from the worktree root).
