# 02 — Species & Combustion sign-offs: patch 0003, the clip gather prototype, MP5 pad patches 0001/0002

Owner: AMR Species & Combustion Lead · Status: review record (DRAFT until the owner accepts). Reviewed against FireX `36975d765f`.

## Patch 0003 (`CHECK_MASS_DENSITY` species loop, two-pass split): ACCEPT
Checked by reading the patch against `Source/mass.f90:868-937` and by a dry-run apply on the FireX tree (applies, hunks offset by 6 lines).
1. **Scatter order.** Pass 2 visits cells K, J, I and applies the seven statements in the original order, so every `DELTA_RHO_ZZ` element receives its additions in the original order. Pass 1 reads only `RHO_ZZ` and `RHOP`, which neither the original nor the patched scatter modifies, so the stored amounts equal the ones the original computed inline.
2. **Rounding.** The original `DELTA - CONST*MASS_N(m)/VC(m)` evaluates as `DELTA - ((CONST*MASS_N(m))/VC(m))`. The patch stores `(CONST*MASS_N(m))/VC(m)` and subtracts it, so the value and the rounding are the same. The division sits between the multiplication and the addition, so no fused multiply-add can form there. This also means the no-FMA flag I asked for earlier is a precaution for the gather, not a requirement for this expression.
3. **Flag.** `SCATTER_FLAG` is reset before the `SOLID` cycle, set to 1 at the same place the old `CLIP_RHO_ZZ(N)=.TRUE.` was, and to 2 only when the scatter happens. `ANY(SCATTER_FLAG>0)` therefore equals the old flag, including the `SUM_MASS_N<=TWO_EPSILON_EB` cells that clip without scattering.
4. **Threading.** The routine says "Do not apply OpenMP to this routine", so the local allocatables are not shared between threads.
Requests (not blocking):
- The two arrays are allocated and freed on every call, about 60 bytes per cell, before knowing whether anything clips. On large meshes this is a per-step cost for the common no-clip case. Prefer a mesh-persistent scratch array, or allocate after a cheap pre-test for any out-of-range cell.
- The density loop (`mass.f90:799-849`) has the same structure. Do the same split in a second patch so the two loops stay consistent.
- Use `MAXVAL(SCATTER_FLAG)>0` instead of `ANY(SCATTER_FLAG>0)` only if a compiler builds a temporary for the latter; gfortran does not.

## Clip gather prototype (`prototypes/p1_mass_shim/clip_gather.f90`, D-031): ACCEPT as the device path for S1
- Gather order is right: sources are visited at (k-1), (j-1), (i-1), self, (i+1), (j+1), (k+1), and the direction index `D` of the target as seen from the source is 3, 2, 1, 0, -1, -2, -3. This reproduces the serial K, J, I scatter order into each element, starting from 0.
- The face mask indices 1..6 map to `MASS_N(-1), (1), (-2), (2), (-3), (3)` as in the original `WALL_INDEX` tests. `CLIP_TERMS` evaluates the same expressions in the same order as `mass.f90:822-839` and `:887-915`, with the species `SOLID`-first and density range-first test order preserved.
- Conditions for the device version: (a) replace `SUM(MASS_N)` by an explicit left-to-right sum over -3..3 with parentheses, because the K2 compiler may reassociate an intrinsic `SUM` (the same rule as the K2 sum rule); (b) the clipped-cell counters (`NCLIP`, `NLO`, `NHI`) must be integer reductions, which are order independent; (c) each cell recomputes `CLIP_TERMS` up to 7 times, which is acceptable because the early return for in-range cells is cheap, but it should be timed on a case where clipping fires in many cells.

## MP5 `Z_TEMP` pad (patches 0001 and 0002): ACCEPT with one physics note
- Confirmed: with `FLUX_LIMITER='MP5'` the upwind branch for `A>0` reads `U(I+IP2)` and the branch for `A<0` reads `U(I+IM1)` (`func.f90:1443-1448`). In the one-face wall calls the missing fourth element is exactly that one, so the result depends on stale scratch memory. The other limiters read only the three assigned elements in the branches the wall guards allow. All hunks of 0001 and 0002 apply to the FireX tree by dry run.
- Note: padding with `0._EB` (as `mass.f90` already does) gives MP5 a fake far-side value of zero for the first off-wall face. The monotonicity limiter bounds the effect, but a nearest-value copy would be a more physical pad. This is low priority because MP5 is not a default; for the device version match the 0 pad so that device and host agree.

## Not signed off
Nothing outstanding in the species and combustion families S1-S4.
