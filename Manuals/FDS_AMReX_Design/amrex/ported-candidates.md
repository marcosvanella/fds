# Ported-status candidates (D-072 (c)) and the fast-math pin check

Owner: AMReX Integration Lead. Status: v0.2. This file records which kernels `tools/ported.toml` sets to ported, what is blocked, and the CMake pin check.

## 1. State of `tools/ported.toml`

**51 kernels are set** (9 `[[ported]]` entries), aligned with the V&V Lead's review (`vv/ported-review.md`, sections 3 to 10), plus D-083 `s5clip_*`. `reviewed_by` says that the V&V Lead reviewed the device runs and the ulp gate; the Architect co-signed the list as it stands (D-077), with the V&V Lead as reviewer of record in place of the kernel owner (D-072 (c)). `set_by` is the Integration Lead. `port_kernel_map.py --strict` accepts all entries (self-check PASS, no PM-08 finding).

| Entry | Kernels | Class | Device evidence as reviewed |
|---|---|---|---|
| 1 | `vpred_us/vs/ws`, `vcorr_u/v/w`, `div2_pred`, `div2_corr`, `vflux_vort_tau` (9) | bitwise | `csmag_32` only (real fields); `dec2_obst` is host only for these. `vpred_vs`, `vcorr_v`, `div2_corr`, `vflux_vort_tau` have no FDS truth |
| 2 | `baro_p_rrho`, `baro_fvx/fvy/fvz` (4) | bitwise | all zero (vacuous) on `csmag_32`; device evidence is the four `dec2_obst` boxes; ghost and solid inputs not physical, no FDS truth |
| 3 | `vflux_fvx/fvy/fvz` | bitwise | real fields and real edge tables; device ran `csmag_32` and `dec2_obst` box 0 only (four boxes is the host result) |
| 4 | `cfl_max`, `vn_max`, `div_extrema` | bitwise | real fields (`csmag_32`, four `dec2_obst` boxes) with `CFL_VELOCITY_NORM` 0 only for `cfl_max`; `div_extrema` only Cartesian with stored divergence off (`CARTVELDIV` all zero); other branches only by the SYNTHETIC rounds 4 to 7 hash run |
| 5 | `rho_d_dzd`, `h_rho_d_dzd` | bitwise | device coverage only: ZZ and TMP real, `RHO_D` and `D_Z` synthetic, no FDS truth; host gate PASS |
| 6 | `cfl_wall_max`, `ulp_gate_passed = true` | libm, 2 ulp | gate `ULP` lines PASS: real 32 walls max 1 ulp (only 19 distinct real values), SYNTHETIC 20,000 walls max 2 ulp; the 2 ulp bound is measured for this toolchain, not a guarantee |
| 7 | `gsfv_central`, `gsfv_godunov`, `gsfv_superbee`, `gsfv_minmod`, `gsfv_charm`, `gsfv_mp5` (6) | bitwise | two device runs on record: stage-1 harness on 1,296 REAL `GET_SCALAR_FACE_VALUE` calls (`shunn3_32`, `dec2_obst`, `csmag_32`; stage1 9715a9e, plan 7.13c, docs 22d7552): device = serial = 4 threads, byte-identical in 7,776 of 7,776 kernel runs, bit-equal to FDS on 1,920,768 of 1,920,768 box cells; and the gate device driver (SYNTHETIC inputs, hash form). Limits: negative `A` with a varying scalar only in `dec2_obst`; `csmag_32` scalar uniform; the one-face callee `GET_SCALAR_FACE_VALUE_PT` is not covered; the caller must pass the address of element (0,0,0). Reviewed and accepted in `vv/ported-review.md` section 9 |
| 8 | `s5clip_density_gather/apply`, `s5clip_species_gather/apply`, `s5clip_species_one`, `s5clip_renorm` (6) | bitwise | D-083; see `tools/ported.toml` |
| 9 | 17 `rad_*` (`rad_extcoe` … `rad_wall_qin_zero`) | bitwise | `mkdev_rad` SYNTHETIC meshes at `25aad430b9`; V&V co-sign D-072 (c) ACCEPT-WITH-LIMITS (`ported-review.md` §10) |

**Held out, not in `ported.toml`:**
- `d_z_max`, `del_rho_d_del_z`: the review accepts them as device coverage (synthetic `RHO_D` and `D_Z`) on the stated condition that the gate's host tier shows a result for both. The latest quick run has no result line for either, so the condition is unmet and they are held out until it is.
- `rho_d_interp`, `rho_d_maxloc_fix`, `dp_div_heat`: the review asks for (1) bitwise host tests in the gate (they are NO-TEST; the map shows `generated`), and (2) a real-field run: dump `RHO_D` and the `D_Z` table, rerun `WP3_DIFF` on the test machine. Cases needed: for `rho_d_interp` the real table with real `TMP`, including the end-of-table clamp; for `rho_d_maxloc_fix` the branch with more than two species (the maximum run has `NS = 2`) and a tie in `MAXV`; for `dp_div_heat` a non-synthetic `H_RHO_D_DZD*`.
- WP1 non-wall and wall kernels (266 and 76 arrays): held until the kernels are named with their truth status and the run is repeated with the device outputs kept (the earlier device outputs were removed, so the reviewer could not re-compare them).

**Open conditions of D-072 (c):**
- `tools/ci_checks.sh --strict` is green on `src-s5gen` at `25aad430b9` (ZS-05 `s5rad_wall1` approved waiver; prior `rad_wall_qin_zero` / `W_SURF_INDEX` FAIL cleared). Stale OPEN CONDITION removed from the live `ported.toml` header.
- The reviewer substitution (V&V Lead in place of the kernel owner) is settled by D-077; nothing open.
- 17 `rad_*`: ACCEPT-WITH-LIMITS — SYNTHETIC `mkdev_rad` hash only; real-field / FDS-truth overlay remains open follow-up (not claimed).

## 2. Evidence (kernel-map names; run on the cc 8.9 test-machine GPU unless stated)

| Kernels | Evidence |
|---|---|
| the 9 kernels of entry 1 | 7.13: device = host serial = host 4 threads on all 54 output files of the `csmag_32` WP3 stages (gate-tool recheck: 0 failing elements); bit-equal to FDS where a value exists (7.9c) |
| the 4 baroclinic kernels | 7.13: `dec2_obst`, four boxes, device = serial = 4 threads on 8 arrays per box |
| `vflux_fvx/fvy/fvz` | 7.13b: `csmag_32` 33,792 of 33,792 faces per direction and `dec2_obst` box 0 (272 / 512 / 272), device = serial = 4 threads = FDS; negative control with neutralised edge tables fails identically |
| `cfl_max`, `vn_max`, `div_extrema` | 7.13a: 155 of 155 outputs; scalars and locations equal in five cases; `div_extrema` checked independently |
| `rho_d_dzd`, `h_rho_d_dzd` | 7.13a neighbours; `RHO_D` synthetic |
| the six `gsfv_*` kernels | 7.13c and the gate device driver, see entry 7 above |
| `cfl_wall_max` | 7.10d and the gate `ULP` lines (`vv-runs/gpu_gate/device_logs/stage1_cflw`): host serial = 4 threads = FDS-build value on 32 of 32 real walls (hot-obstruction variant of `dec2_obst`, scratch input); device within 1 ulp (4 of 32 differ); SYNTHETIC sweep max 2 ulp |

Not set: loops without a kernel or marker, and the held-out kernels above. (An earlier version of this file listed the six `gsfv_*` kernels as "no device run"; that was wrong, two device runs are on record, see entry 7 and `tools/device_runs.toml`.)

## 3. Fast-math pin (D-068): applied on the generator worktree

The pin patch (`tools/patches/cmake-fastmath-pin.patch`) is **applied** on the generator worktree (`s5-gen`). Commit `3546f746a7` (author a dedicated project identity): new `amrex/cmake/AMReXCudaPin.cmake` and `amrex/cmake/amrex-cuda-gpu.cmake`; `include` + `amrex_cuda_fastmath_check()` in `amrex/s4_mass/CMakeLists.txt`, `amrex/s4_mass/s4d/CMakeLists.txt`, and `amrex/s4_mass/s5_gen/gpu/CMakeLists.txt`. Reverse-apply of the patch is clean against that tree. Local only; not pushed.

Scratch-copy check that preceded the apply (still valid as evidence of the pin effect):

- `kernel_lint.py --strict` on the unpatched tree: BF-05 FAIL, BF-07 FAIL.
- Same on the patched tree: BF-05 NOTE ("D-068 pin found: FORCE-OFF in amrex/cmake/AMReXCudaPin.cmake, failing guard in the same file"); BF-07 no finding (pass).
- `zone_sum_order.py --strict`: PASS. `port_kernel_map.py --strict`: self-check PASS.

GNU verify **PASS** on s5-gen tip `dcada34a2e` (pin `3546f746a7` ancestor): BF-05 NOTE + BF-07 pass; configure ON→FATAL_ERROR. Evidence: `vv-runs/d068-fastmath-pin/RESULT.md`. A separate K2-03 / CL-07 finding on an older `s5gen_rad_wall_qin_zero` golden is not part of the pin; a fresh `rad_k2.F90` from `markers/rad_kernels.toml` at this tip passes `k2_ci_check` and `kernel_lint` (device-address list includes `W_SURF_INDEX`).
