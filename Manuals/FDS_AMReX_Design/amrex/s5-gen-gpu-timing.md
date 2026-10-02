# Generated K2 kernels on the GPU: bitwise status and timing (34 kernels)

Branch `s5-gen`, HEAD `6c87023310` (includes the callee-directive batch `0e8b78251b`). All numbers below are **run** results on a cc 8.9 test-machine GPU (8 GB) with nvfortran 26.9 unless a line says **compiled** (compile-time information only). Synthetic inputs (random fields with realistic ranges, synthetic wall mesh); no FDS restart data.

## Setup

- Regeneration: `s5gen.py` (venv with fparser) prints `wrote ... (34 kernels, 4 callees)`; `git status --porcelain -- amrex` is empty afterwards (reproducible). Generated file, header and report unchanged.
- Flags (compiled): Fortran `-O2 -mp=gpu -gpu=cc89,nofma -Minline -Minfo=mp`; nvcc `--fmad=false`; AMReX CUDA build without fast math (existing install, not rebuilt). The generated-kernel harness additionally uses `-gpu=mem:managed` (see below).
- Two builds of the generated file: `-DS5_CALLEE_DPD` (default, callee-containing loops use `target teams distribute parallel do`) and `-DS5_CALLEE_BIND` (`target teams loop ... bind(teams,parallel)`).
- Two harnesses:
  1. `s5_perf` (committed S4d-style C++ harness, device-memory arrays): the 7 original nests, K1 (AMReX ParallelFor), hand-made K2 (`perf.gen=0`) and generated K2 (`perf.gen=1`). 10 warm-ups, 50 timed reps, device sync before and after each call, wall clock; median/p90.
  2. Fortran harness (new, not in git): the committed round-2/3 bitwise driver turned into a large-size driver by a script (`mkperf.py`): 27 kernels (18 round 2, 9 round 3), species 1, NL=0, no internal-aux walls. Arrays are CUDA managed memory (`-gpu=mem:managed`) so the host reference and the `is_device_ptr` kernels share them; each kernel is warmed up (10 calls, which also migrates its arrays to the device) and then timed 60 times with device sync before and after each call. The host never touches the arrays during the timed loop. No per-rep restore: in-place kernels see their own output again (the same applies to the S4d protocol for in-place nests).
- Bitwise method for the 27 kernels: the reference is the upstream loop text compiled with gfortran `-O2 -ffp-contract=off` (host, serial). Both programs (gfortran host, nvfortran GPU) run the same setup and print a 128-bit rotate-xor hash of every array and wall table after each kernel; a kernel is EQUAL if the GPU hash equals the gfortran reference hash and the input hashes match. This is exact-bit equality of every element (a +0/-0 difference would be a DIFF). It is a hash, not an element-by-element compare, at 128 and 192; the committed driver compares element by element at the small sizes.

## 1. Bitwise status (run)

| harness | build | size | result |
|---|---|---|---|
| `s5_perf`, 7 nests, gen0 and gen1 | DPD and BIND | 96x40x56, 128^3, 192^3 | all `BITWISE K2 vs CPU <nest>` lines EQUAL; `BITWISE SUMMARY ... ALL EQUAL (value), +0/-0-only elements 0`; `S4d exit 0` in all 12 runs |
| Fortran harness, 27 kernels | DPD | 96x40x56, 128^3, 192^3 | 27/27 EQUAL, 0 DIFF, 0 input mismatches |
| Fortran harness, 27 kernels | BIND | 96x40x56, 128^3, 192^3 | 27/27 EQUAL, 0 DIFF, 0 input mismatches |
| committed driver (408 cases) on gfortran host | - | 1x1x1 to 24x20x16 | 408/408 pass, 0 vacuous |

The 7 original nests are 7 of the 34 kernels, so all 34 kernels are checked at 128^3 and 192^3. 192^3 fits in the 8 GB card (about 2.6 GB of device memory in use during the run).

Host-compiler observation (run, not a GPU result): the committed bitwise driver compiled entirely with nvfortran (host reference too) reports 33 failing cases of 408 for three kernels (`zzs_pred`, `zz_corr`, `wall_strain_rate`), while the same driver with gfortran passes 408/408. The failing side is the nvfortran-compiled upstream reference text: it contains unparenthesized sums, which nvfortran reassociates (allowed by the standard); the generated kernels carry explicit left-to-right parentheses and match the gfortran reference. The nvfortran host reference differs from the gfortran reference in exactly these three kernels at all sizes. Consequence: the CPU reference of a GPU comparison must be built with a compiler that keeps the written order (gfortran), as in S4d.

## 2. Timings (run), median/p90 in microseconds

Ratios are medians. "gen/hand" compares the generated and the hand-made K2 of the same build; "B/D" is BIND over DPD.

## Seven original nests (s5_perf harness: device-memory arrays, K1 AMReX, hand-made K2, generated K2); median/p90 us
| n | nest | K1 | K2 hand | gen DPD | gen BIND | gen DPD/hand | gen BIND/hand | K2 hand/K1 |
|---|---|---|---|---|---|---|---|---|
| 96 | conductivity | 58.0/58.2 | 55.4/55.7 | 53.4/53.8 | 53.0/53.3 | 0.96 | 0.96 | 0.95 |
| 96 | cp_rhg | 38.1/38.4 | 31.9/32.1 | 31.9/32.1 | 32.8/33.0 | 1.00 | 1.03 | 0.84 |
| 96 | dp_kdtd | 22.9/23.1 | 21.7/22.0 | 21.4/21.7 | 21.5/21.8 | 0.99 | 0.99 | 0.95 |
| 96 | dp_species | 44.9/45.2 | 43.9/44.1 | 43.9/44.1 | 43.4/43.6 | 1.00 | 0.99 | 0.98 |
| 96 | h_rho_d_dzd | 47.7/47.8 | 47.4/47.6 | 47.4/47.6 | 46.8/47.0 | 1.00 | 0.99 | 0.99 |
| 96 | kdtd | 29.7/29.9 | 28.8/29.0 | 28.6/28.8 | 28.6/28.8 | 0.99 | 0.99 | 0.97 |
| 96 | rho_d_dzd | 29.6/29.8 | 29.0/29.3 | 28.5/28.7 | 28.6/28.8 | 0.98 | 0.99 | 0.98 |
| 128 | conductivity | 610.4/613.2 | 582.5/584.2 | 566.9/568.6 | 567.2/567.7 | 0.97 | 0.97 | 0.95 |
| 128 | cp_rhg | 509.1/510.1 | 499.6/500.8 | 506.1/508.9 | 498.5/500.1 | 1.01 | 1.00 | 0.98 |
| 128 | dp_kdtd | 531.3/531.9 | 532.8/534.6 | 532.2/533.7 | 532.2/534.2 | 1.00 | 1.00 | 1.00 |
| 128 | dp_species | 643.3/643.9 | 643.0/644.0 | 643.2/644.0 | 642.7/643.2 | 1.00 | 1.00 | 1.00 |
| 128 | h_rho_d_dzd | 518.9/521.5 | 519.5/521.6 | 519.8/521.8 | 525.4/528.0 | 1.00 | 1.01 | 1.00 |
| 128 | kdtd | 385.5/388.1 | 373.7/376.3 | 381.9/384.6 | 382.9/385.9 | 1.02 | 1.02 | 0.97 |
| 128 | rho_d_dzd | 380.4/384.6 | 384.8/388.1 | 390.2/394.0 | 390.7/394.4 | 1.01 | 1.02 | 1.01 |
| 192 | conductivity | 2069.5/2070.8 | 1843.4/1845.7 | 1830.1/1832.7 | 1860.3/1861.3 | 0.99 | 1.01 | 0.89 |
| 192 | cp_rhg | 1659.4/1663.4 | 1655.5/1659.1 | 1657.6/1661.2 | 1657.2/1662.1 | 1.00 | 1.00 | 1.00 |
| 192 | dp_kdtd | 1714.7/1721.1 | 1713.3/1717.5 | 1712.6/1720.6 | 1715.1/1723.0 | 1.00 | 1.00 | 1.00 |
| 192 | dp_species | 2060.6/2062.5 | 2058.6/2059.9 | 2058.7/2060.6 | 2058.0/2060.0 | 1.00 | 1.00 | 1.00 |
| 192 | h_rho_d_dzd | 1711.7/1714.8 | 1709.2/1712.2 | 1745.2/1747.5 | 1715.7/1718.4 | 1.02 | 1.00 | 1.00 |
| 192 | kdtd | 1242.0/1245.1 | 1242.6/1246.1 | 1273.2/1280.6 | 1271.9/1277.1 | 1.02 | 1.02 | 1.00 |
| 192 | rho_d_dzd | 1233.9/1236.6 | 1243.9/1246.3 | 1278.0/1284.5 | 1276.9/1282.7 | 1.03 | 1.03 | 1.01 |

## 27 round-2/3 kernels (Fortran harness, managed memory resident on device); median/p90 us, BIND/DPD = median ratio
| kernel | n128 DPD | n128 BIND | B/D | n192 DPD | n192 BIND | B/D | s96x40x56 DPD | s96x40x56 BIND | B/D |
|---|---|---|---|---|---|---|---|---|---|
| zzs_pred | 3092.1/3094.0 | 3092.1/3094.0 | 1.00 | 10324.0/10330.0 | 10324.9/10330.0 | 1.00 | 240.8/410.1 | 241.0/430.8 | 1.00 |
| adv_flux_store | 2678.9/2687.9 | 2675.0/2683.9 | 1.00 | 8926.2/8949.0 | 8927.8/8943.8 | 1.00 | 197.9/199.7 | 198.8/201.0 | 1.00 |
| rhos_sum | 434.8/435.9 | 434.2/435.1 | 1.00 | 1431.0/1433.2 | 1431.0/1432.9 | 1.00 | 12.9/13.1 | 12.9/13.1 | 1.00 |
| rsum_pred | 687.9/691.2 | 689.0/691.9 | 1.00 | 2286.9/2290.0 | 2285.9/2289.1 | 1.00 | 53.0/53.9 | 53.0/53.2 | 1.00 |
| zz_corr | 3506.0/3507.9 | 3505.9/3508.1 | 1.00 | 11715.9/13412.0 | 11714.0/11721.2 | 1.00 | 257.9/711.2 | 256.1/441.1 | 0.99 |
| adv_flux_avg | 3389.1/3390.7 | 3875.9/3877.1 | 1.14 | 11298.9/11301.0 | 11299.1/11301.1 | 1.00 | 283.0/284.0 | 283.0/285.2 | 1.00 |
| rho_sum | 433.0/434.0 | 433.9/434.9 | 1.00 | 1429.8/1431.9 | 1430.0/1431.0 | 1.00 | 12.9/13.1 | 12.8/13.1 | 0.99 |
| flux_mw_fix | 2616.0/2923.0 | 2921.9/2923.0 | 1.12 | 8713.0/8714.9 | 8712.1/8715.9 | 1.00 | 377.9/510.9 | 332.1/797.1 | 0.88 |
| flux_mw_fix_zz | 2890.1/2924.0 | 2923.0/2924.9 | 1.01 | 8714.2/8717.0 | 8715.1/9731.0 | 1.00 | 373.1/799.9 | 365.0/390.0 | 0.98 |
| mu_dns | 510.0/511.9 | 525.0/525.9 | 1.03 | 1668.0/1669.0 | 1729.0/1730.0 | 1.04 | 69.9/304.0 | 58.0/173.0 | 0.83 |
| up_deardorff | 456.1/459.9 | 459.9/463.9 | 1.01 | 1502.0/1509.9 | 1503.0/1509.2 | 1.00 | 15.0/16.0 | 15.0/16.0 | 1.00 |
| kres | 288.0/288.1 | 288.0/289.0 | 1.00 | 946.0/947.0 | 946.1/947.0 | 1.00 | 25.0/25.1 | 25.0/25.1 | 1.00 |
| wall_uvw_interp | 7.0/7.2 | 7.0/7.2 | 1.00 | 9.1/10.0 | 9.1/10.0 | 1.00 | 5.9/6.0 | 6.0/6.2 | 1.02 |
| wall_uvw_interp_corr | 7.0/7.1 | 6.9/7.2 | 0.99 | 9.0/10.0 | 9.1/10.0 | 1.01 | 5.9/6.2 | 5.9/6.2 | 1.00 |
| wall_up_ghost | 15.1/16.2 | 15.0/16.0 | 0.99 | 52.9/54.1 | 52.0/53.0 | 0.98 | 7.0/7.9 | 7.0/7.8 | 1.00 |
| wall_kp_ghost | 10.0/10.1 | 10.0/10.1 | 1.00 | 15.0/16.0 | 15.1/16.0 | 1.01 | 6.0/6.2 | 6.0/6.2 | 1.00 |
| wall_rho_d_dzdn | 10.1/11.2 | 10.1/11.0 | 1.00 | 16.0/16.9 | 16.0/17.0 | 1.00 | 6.9/7.2 | 6.9/7.1 | 1.00 |
| wall_un_store | 7.9/8.1 | 7.9/8.1 | 1.00 | 10.9/11.2 | 10.9/11.0 | 1.00 | 5.9/6.0 | 5.9/6.0 | 1.00 |
| wall_coriolis_ghost | 15.0/15.1 | 15.0/15.1 | 1.00 | 57.0/58.0 | 57.2/58.2 | 1.00 | 6.9/7.2 | 6.9/7.2 | 1.00 |
| wall_hs_bt | 8.1/8.8 | 8.1/9.1 | 1.00 | 11.9/12.2 | 11.9/12.2 | 1.00 | 5.9/6.0 | 5.9/6.2 | 1.00 |
| wall_us_pred | 10.0/11.0 | 10.0/11.0 | 1.00 | 16.9/17.2 | 16.2/17.1 | 0.96 | 6.9/7.2 | 6.9/7.1 | 1.00 |
| wall_strain_rate | 56.0/57.0 | 56.0/57.0 | 1.00 | 118.0/119.0 | 118.1/119.0 | 1.00 | 37.0/206.0 | 36.9/89.9 | 1.00 |
| wall_b2_work1 | 5.9/6.2 | 6.0/6.2 | 1.02 | 7.0/7.2 | 7.0/7.9 | 1.00 | 5.0/5.0 | 5.0/5.0 | 1.00 |
| rho_z_p_mass | 235.1/238.0 | 234.8/237.9 | 1.00 | 787.9/792.9 | 787.0/792.0 | 1.00 | 9.0/10.0 | 9.0/10.0 | 1.00 |
| rho_z_p_divg | 235.1/237.9 | 234.8/237.9 | 1.00 | 766.0/789.0 | 785.8/791.1 | 1.03 | 9.1/10.0 | 9.1/10.0 | 1.00 |
| rho_zz_clip_assign | 285.1/286.1 | 285.0/285.8 | 1.00 | 940.1/941.0 | 940.1/941.1 | 1.00 | 10.0/10.9 | 10.0/11.0 | 1.00 |
| delta_rho_zz_zero | 37.2/38.1 | 37.2/38.2 | 1.00 | 234.8/236.0 | 234.1/236.0 | 1.00 | 8.1/8.2 | 8.1/8.1 | 1.00 |
| sum of medians | 21927 | 22768 | 1.04 | 72008 | 72088 | 1.00 | 2055 | 1988 | 0.97 |

Notes:
- Generated over hand-made, 7 nests: 0.96 to 1.03 in every case (no value above 1.1). Hand-made K2 over K1: 0.84 to 1.01 (K2 faster at 96x40x56 and in `conductivity`, equal elsewhere).
- Only kernels with a callee or a local array differ between the two builds. All other kernels compile to the same code in both builds, so their B/D differences (adv_flux_avg 1.14, flux_mw_fix 1.12 at 128^3) are timing noise: the DPD runs show a bimodal median (about 2.62 ms vs 2.92 ms), and the 192^3 and the other pairs are 1.00.
- Real switch effect: `mu_dns` (calls `GET_VISCOSITY`, local array `ZZ_GET`): BIND is 1.03 to 1.04 slower at 128^3 and 192^3 (compiled: 96 registers and local memory with BIND, 40 registers with DPD). The 7 original nests differ by at most 3% between the two builds.
- No outlier above ratio 2. Nothing slower than 1.1 against hand-made (the 1.12/1.14 entries above are noise between identical code).
- The small size has p90 values several times the median for a few kernels (`zzs_pred`, `zz_corr`, `wall_strain_rate`, `mu_dns`): launch-latency-bound kernels of 5 to 250 us, so a single slow rep moves p90.
- The wall-table kernels cost 5 to 120 us at 128^3 and 192^3 (about 98k and 221k external walls, plus 40 internal): they are launch-latency-bound and negligible next to the cell kernels (0.2 to 12 ms).

## 3. Registers (compiled, `cuobjdump --dump-resource-usage`, no stack or local memory in any kernel except as stated)

| kernel | DPD | BIND |
|---|---|---|
| zzs_pred, zz_corr, flux_mw_fix | 64 | 64 |
| adv_flux_store, adv_flux_avg | 40 | 40 |
| mu_dns | 40 | 96 (local memory used for `zz_get`) |
| wall_strain_rate | 52 | 52 |
| wall_up_ghost / wall_coriolis_ghost | 24 / 22 | 24 / 22 |
| rho_zz_clip_assign | 37 | 37 |

## 4. The 9 round-3 wall-table kernels

Kernels: `wall_coriolis_ghost`, `wall_hs_bt`, `wall_us_pred`, `wall_strain_rate`, `wall_b2_work1` (whole wall loops) and `rho_z_p_mass`, `rho_z_p_divg`, `rho_zz_clip_assign`, `delta_rho_zz_zero` (cell sub-nests).

- Bitwise (run): 9/9 EQUAL at 96x40x56, 128^3 and 192^3, in both builds.
- Timing (run, median at 128^3 / 192^3, us): wall_coriolis_ghost 15.0/57.0, wall_hs_bt 8.1/11.9, wall_us_pred 10.0/16.9, wall_strain_rate 56.0/118.0, wall_b2_work1 5.9/7.0, rho_z_p_mass 235/788, rho_z_p_divg 235/766, rho_zz_clip_assign 285/940, delta_rho_zz_zero 37/235.
- `-Minfo=mp` (compiled): the wall loops use `target teams distribute parallel do` in both builds ("Generating nvkernel ... GPU kernel", no schedule line, so no `threads(...)` line is printed); the cell sub-nests use `target teams loop` and print `Loop parallelized across teams, threads(128) collapse(3)`. No warning, error or ICE in either build. Implicit privates and "Loop run sequentially" lines appear only for the per-cell species loops (dot1/dot2 reductions run in order inside one thread).
- Table arguments (compiled and run): the flat tables (`BC_IIG`, `BC_IOR`, `W_BOUNDARY_TYPE`, `B1_*`, `B2_WORK1`, ...) are explicit-shape arrays `(NWE+NWI)` or `(NWE)`; no derived types reach the device. No problem found.
- Unique writes: `wall_kp_ghost`, `wall_coriolis_ghost`, `wall_hs_bt` carry the `unique` marker (each wall writes a distinct element); `wall_us_pred` and `wall_strain_rate` carry `idempotent` (several walls may reach one element and store identical bits). No `omp atomic` is generated. The synthetic mesh used here is race-free by construction (internal walls on distinct faces, thin pairs, edge-line walls), so the runs show no race; a real mesh with corner or edge cells reached by two external walls is covered only by the generator's static `idempotent`/`unique` check, not by this run.

## 5. nvfortran issues

- No new nvfortran defect found: all 34 kernels compile and run bitwise equal in both switch builds; no reproducer needed.
- Reassociation of unparenthesized sums in the host reference (section 1): a standard-conforming behaviour, same cause as the S4d finding, not a defect. The generator's explicit parentheses are what make the GPU kernels match.
- Test-driver portability (not a compiler defect): the committed bitwise driver draws random numbers in `if (urand() < 0.5 .and. nb + 2 <= NWI)` and `(urand() < 0.15 .and. nx*ny*nz > 1)`. nvfortran may skip the call when the second operand is false, so the random stream (and hence the synthetic mesh) differs from gfortran for some sizes. The large-size driver splits these statements. This does not affect the committed driver, which uses one compiler per run.

## 6. Recommendation

Keep `S5_CALLEE_DPD` as the default. Across 34 kernels and two sizes it is never slower than `S5_CALLEE_BIND` by more than noise and it is faster for `mu_dns` (BIND: 96 registers, local memory, 3 to 4% slower); BIND gives no gain on any kernel, and it needs an NVHPC-only clause. BIND stays available as a switch for future callee loops that DPD cannot express.

## 7. Reproduce

Scripts and logs are in `prototypes/s4_cuda_mass/s5_gen_34/` (not in git): `build_perf.sh` (7-nest harness, both builds), `build_f.sh` (Fortran harness, both builds, plus the gfortran host program), `run_all.sh` and `run_size.sh` (runs), `mkperf.py` (driver generator), `cmp_hash.py`, `an34.py` (tables). The committed code was not modified.
