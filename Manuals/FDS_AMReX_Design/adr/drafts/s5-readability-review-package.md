# S5 readability review package: K1 (restricted C++) vs K2 (Fortran + OpenMP `target`)

Status: sent for review (ADR-001 v0.5). Reviewers (owner decision D-048): the project owner, the AMR Species & Combustion Lead, the FDS Legacy Mapper. The owner signs off the outcome. Reading time: about 20 minutes. Nothing here needs a GPU or a build.

## 1. What you are choosing
Every FDS time-step kernel must run on an NVIDIA GPU (D-027). C++ stays limited to the driver and AMReX glue (D-043). The physics kernels can be written in one of two styles, and the choice is yours on readability and maintenance for Fortran-only developers:

- **K1, restricted "Fortran-style C++".** One `amrex::ParallelFor` per kernel over `(i,j,k)`, FDS variable names, plain loops and `if`, no templates or classes. Built with the CUDA compiler AMReX already supports.
- **K2, Fortran with simple OpenMP `target`.** The kernel stays Fortran. One `!$omp target teams loop collapse(3)` directive per loop nest, built with `nvfortran`. A thin C++ wrapper passes each array's pointer and bounds.

Both were written for the same mass kernel (species update, sum, clip, post-clip: 12 kernels) and both give **identical results**: bitwise equal to the FDS/P1 reference on the CPU (60 comparisons, 69 checks at 1, 2 and 4 ranks) and bitwise equal on an NVIDIA GPU (29 of 29 comparisons each). Source: branch `s4-cuda-mass`, `amrex/s4_mass/` at `b51f4361b3`; numbers in `docs/amrex/s4-cuda-mass-findings.md`.

## 2. Facts that do not depend on taste
| | K1 | K2 |
|---|---|---|
| Code size for the 12 kernels | 279 lines | 372 lines Fortran + 96 lines C++ wrappers |
| Port from the FDS Fortran | line-by-line, but new syntax (lambda, `Array4`, 0-based component index) | mechanical from K1; keeps Fortran syntax, adds bounds plumbing and directive clauses |
| Failures found only on the GPU | none | two (below) |
| Time per step, one sample (SUPERBEE, 32³, 8 steps) | loop 6.1 ms | loop 8.3 ms (not tuned; 144 host sync points per case) |
| Toolchain | AMReX's supported CUDA path | `nvfortran` plus `nvc++`; AMReX does not officially support pairing them (not tested upstream) |
| Compiler quirks met | none | `nvfortran` 26.9 rejects `has_device_addr` (we use `is_device_ptr`); `c_f_pointer` in a target region does not compile; needs `-Minline` and a special link step |

The two K2-only failures, both fixed with no change to the numerics:
1. **A private copy of a loop-invariant value** (`QMIN = QMIN_IN` at the top of the loop, listed in `private`) made every GPU thread use a wrong value: every cell was flagged as clipped and total density grew by 306 %. The CPU builds were correct, so the CPU tests could not have caught it.
2. **An unparenthesised sum** was reordered by the compiler at `-O2`, giving last-bit differences from FDS (results stay within tolerance but are no longer bitwise equal).

They led to two K2 coding rules (ADR-001 "K2 coding rules"): no private copy of loop-invariant values, and parenthesise every sum whose order matters. Both need the GPU (or the `-Minfo` output) to check. K1 has no equivalent rule, because CUDA C++ follows the source order strictly with the flags in ADR-001 "GPU build flags".

## 3. The code
Read the FDS original first, then each style. All three compute the same species update for one cell, `RHS` and `ZZS` (predictor) or `ZZ` (corrector).

### 3.1 FDS original (`mass.f90:440-455`, predictor)
```fortran

   !$OMP DO PRIVATE(N,I,J,K,RHS)
   DO N=1,N_TOTAL_SCALARS
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
               RHS = - DEL_RHO_D_DEL_Z__0(I,J,K,N) &
                   + (FX(I,J,K,N)*UU(I,J,K)*R(I) - FX(I-1,J,K,N)*UU(I-1,J,K)*R(I-1))*RDX(I)*RRN(I) &
                   + (FY(I,J,K,N)*VV(I,J,K)      - FY(I,J-1,K,N)*VV(I,J-1,K)       )*RDY(J)        &
                   + (FZ(I,J,K,N)*WW(I,J,K)      - FZ(I,J,K-1,N)*WW(I,J,K-1)       )*RDZ(K)
               ZZS(I,J,K,N) = RHO(I,J,K)*ZZ(I,J,K,N) - DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO
```

### 3.2 K1 (`s4_mass_k1.H`, species update)
```cpp
inline void K1_DENSITY_UPDATE (const Box& bx, bool PREDICTOR, Real DT, Real RDX, Real RDY, Real RDZ,
                               Array4<Real const> const& RHO, Array4<Real> const& ZZ,
                               Array4<Real> const& RHOS, Array4<Real> const& ZZS,
                               Array4<Real const> const& FX, Array4<Real const> const& FY, Array4<Real const> const& FZ,
                               Array4<Real const> const& UU, Array4<Real const> const& VV, Array4<Real const> const& WW,
                               Array4<Real const> const& DEL_RHO_D_DEL_Z, Array4<int const> const& MASK,
                               MassConstants const& C)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (MASK(i,j,k,0) != 0) return;                                      // CELL%SOLID
        for (int N = 0; N < C.N_TRACKED_SPECIES; ++N) {
            const Real RHS = - DEL_RHO_D_DEL_Z(i,j,k,N)
                             + (FX(i+1,j,k,N+1)*UU(i+1,j,k) - FX(i,j,k,N+1)*UU(i,j,k))*RDX
                             + (FY(i,j+1,k,N+1)*VV(i,j+1,k) - FY(i,j,k,N+1)*VV(i,j,k))*RDY
                             + (FZ(i,j,k+1,N+1)*WW(i,j,k+1) - FZ(i,j,k,N+1)*WW(i,j,k))*RDZ;
            if (PREDICTOR) {
                ZZS(i,j,k,N) = RHO(i,j,k)*ZZ(i,j,k,N) - DT*RHS;
            } else {
                ZZ(i,j,k,N) = 0.5*( RHO(i,j,k)*ZZ(i,j,k,N) + RHOS(i,j,k)*ZZS(i,j,k,N) - DT*RHS );
            }
        }
    });
}

```

### 3.3 K2 (`s4_mass_k2.F90`, same kernel; first the argument list, then the loop)
```fortran
subroutine s4k2_density_update(lo, hi, PREDICTOR, DT, RDX, RDY, RDZ, RHO, rl, rh, ZZ, zl, zh, RHOS, sl, sh, ZZS, ssl, ssh, &
                               FX, fxl, fxh, FY, fyl, fyh, FZ, fzl, fzh, nf, UU, uul, uuh, VV, vvl, vvh, WW, wwl, wwh, &
                               DEL, dl, dh, MASK, ml, mh, nm, NS) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), PREDICTOR, rl(3), rh(3), zl(3), zh(3), sl(3), sh(3), ssl(3), ssh(3), &
                              fxl(3), fxh(3), fyl(3), fyh(3), fzl(3), fzh(3), nf, uul(3), uuh(3), vvl(3), vvh(3), &
                              wwl(3), wwh(3), dl(3), dh(3), ml(3), mh(3), nm, NS
real(c_double), intent(in), value :: DT, RDX, RDY, RDZ
real(c_double), intent(inout) :: RHO(rl(1):rh(1), rl(2):rh(2), rl(3):rh(3))
real(c_double), intent(inout) :: ZZ(zl(1):zh(1), zl(2):zh(2), zl(3):zh(3), 0:NS-1)
real(c_double), intent(inout) :: RHOS(sl(1):sh(1), sl(2):sh(2), sl(3):sh(3))
real(c_double), intent(inout) :: ZZS(ssl(1):ssh(1), ssl(2):ssh(2), ssl(3):ssh(3), 0:NS-1)
real(c_double), intent(in) :: FX(fxl(1):fxh(1), fxl(2):fxh(2), fxl(3):fxh(3), 0:nf-1)
real(c_double), intent(in) :: FY(fyl(1):fyh(1), fyl(2):fyh(2), fyl(3):fyh(3), 0:nf-1)
real(c_double), intent(in) :: FZ(fzl(1):fzh(1), fzl(2):fzh(2), fzl(3):fzh(3), 0:nf-1)
real(c_double), intent(in) :: UU(uul(1):uuh(1), uul(2):uuh(2), uul(3):uuh(3))
real(c_double), intent(in) :: VV(vvl(1):vvh(1), vvl(2):vvh(2), vvl(3):vvh(3))
real(c_double), intent(in) :: WW(wwl(1):wwh(1), wwl(2):wwh(2), wwl(3):wwh(3))
real(c_double), intent(in) :: DEL(dl(1):dh(1), dl(2):dh(2), dl(3):dh(3), 0:NS-1)
integer(c_int), intent(in) :: MASK(ml(1):mh(1), ml(2):mh(2), ml(3):mh(3), 0:nm-1)
integer :: i, j, k, N, i0, i1, j0, j1, k0, k1
real(eb) :: RHS, S
i0 = lo(1); i1 = hi(1); j0 = lo(2); j1 = hi(2); k0 = lo(3); k1 = hi(3)
#if defined(__NVCOMPILER)
! nvfortran 26.9 rejects has_device_addr (syntax error); ADR-001 fallback is_device_ptr on the
! explicit-shape dummies (OpenMP 5.1: non-c_ptr is_device_ptr items are treated as has_device_addr)
!$omp target teams loop collapse(3) is_device_ptr(RHO, ZZ, RHOS, ZZS, FX, FY, FZ, UU, VV, WW, DEL, MASK) private(N, RHS, S)
#else
!$omp target teams loop collapse(3) has_device_addr(RHO, ZZ, RHOS, ZZS, FX, FY, FZ, UU, VV, WW, DEL, MASK) private(N, RHS, S)
#endif
do k = k0, k1
   do j = j0, j1
      do i = i0, i1
         if (MASK(i,j,k,0) /= 0) cycle
         do N = 0, NS-1
            ! S4b: explicit parentheses = the left-to-right order of FDS/K1. Fortran lets a compiler reassociate an
            ! unparenthesised sum; nvfortran -O2 does (last-bit differences vs gfortran/K1), parentheses must be honoured.
            RHS = ((- DEL(i,j,k,N) &
                  + (FX(i+1,j,k,N+1)*UU(i+1,j,k) - FX(i,j,k,N+1)*UU(i,j,k))*RDX) &
                  + (FY(i,j+1,k,N+1)*VV(i,j+1,k) - FY(i,j,k,N+1)*VV(i,j,k))*RDY) &
                  + (FZ(i,j,k+1,N+1)*WW(i,j,k+1) - FZ(i,j,k,N+1)*WW(i,j,k))*RDZ
            if (PREDICTOR /= 0) then
               ZZS(i,j,k,N) = RHO(i,j,k)*ZZ(i,j,k,N) - DT*RHS
            else
               ZZ(i,j,k,N) = 0.5_eb*( (RHO(i,j,k)*ZZ(i,j,k,N) + RHOS(i,j,k)*ZZS(i,j,k,N)) - DT*RHS )
            endif
         enddo
         S = 0._eb
         if (PREDICTOR /= 0) then
            do N = 0, NS-1
               S = S + ZZS(i,j,k,N)
            enddo
            RHOS(i,j,k) = S
         else
            do N = 0, NS-1
               S = S + ZZ(i,j,k,N)
            enddo
            RHO(i,j,k) = S
         endif
      enddo
   enddo
enddo
end subroutine s4k2_density_update

! ZZ = RHO_ZZ/RHOP; RSUM = R0*DOT_PRODUCT(MWR_Z,ZZ); TMP = PBAR/(RSUM*RHOP) (mass.f90:526-578 / 708-760)
```
The C++ side of K2 (`s4_mass_k2.H`) that supplies the pointers and bounds:
```cpp
inline void k2_density_update (const amrex::Box& bx, bool PREDICTOR, double DT, double RDX, double RDY, double RDZ,
                               amrex::FArrayBox& RHO, amrex::FArrayBox& ZZ, amrex::FArrayBox& RHOS, amrex::FArrayBox& ZZS,
                               const amrex::FArrayBox& FX, const amrex::FArrayBox& FY, const amrex::FArrayBox& FZ,
                               const amrex::FArrayBox& UU, const amrex::FArrayBox& VV, const amrex::FArrayBox& WW,
                               const amrex::FArrayBox& DEL, const amrex::IArrayBox& MASK, s4k1::MassConstants const& C) {
    auto b = s4k2w::b3(bx), r = S4F(RHO), z = S4F(ZZ), s = S4F(RHOS), zs = S4F(ZZS), fx = S4F(FX), fy = S4F(FY), fz = S4F(FZ),
         u = S4F(UU), v = S4F(VV), w = S4F(WW), d = S4F(DEL), m = S4F(MASK);
    int pr = PREDICTOR ? 1 : 0, nf = FX.nComp(), nm = MASK.nComp(), ns = C.N_TRACKED_SPECIES;
    s4k2_density_update(b.lo, b.hi, &pr, DT, RDX, RDY, RDZ, RHO.dataPtr(), r.lo, r.hi, ZZ.dataPtr(), z.lo, z.hi, RHOS.dataPtr(), s.lo, s.hi,
                        ZZS.dataPtr(), zs.lo, zs.hi, FX.dataPtr(), fx.lo, fx.hi, FY.dataPtr(), fy.lo, fy.hi, FZ.dataPtr(), fz.lo, fz.hi, &nf,
                        UU.dataPtr(), u.lo, u.hi, VV.dataPtr(), v.lo, v.hi, WW.dataPtr(), w.lo, w.hi, DEL.dataPtr(), d.lo, d.hi,
                        MASK.dataPtr(), m.lo, m.hi, &nm, &ns);
}
inline void k2_post_clip (const amrex::Box& bx, const amrex::FArrayBox& RHOP, amrex::FArrayBox& ZZP, amrex::FArrayBox& TMP,
                          const amrex::IArrayBox& MASK, s4k1::MassConstants const& C) {
```

### 3.4 A second kernel: the clip terms, top half (K1 then K2)
K1:
```cpp
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        for (int d = 0; d < 7; ++d) { T(i,j,k,d) = 0.0; }
        CF(i,j,k) = 0;
        const Real QMIN = QMIN_IN;
        const Real QMAX = SPECIES ? RHOP(i,j,k) : QMAX_IN;                   // RHO_ZZ_MAX = RHOP(I,J,K) (:887)
        const Real QC = Q(i,j,k,NQ);
        if (SPECIES) {                                                       // :885, :888 (SOLID first)
            if (MASK(i,j,k,0) != 0) return;
            if (QC >= QMIN && QC <= QMAX) return;
        } else {                                                             // :809, :811 (range first)
            if (QC >= QMIN && QC <= QMAX) return;
            if (MASK(i,j,k,0) != 0) return;
        }
        Real Q_CUT, SIGN_FACTOR;
        int flag = CF_CLIPPED;
        if (QC < QMIN) { Q_CUT = QMIN; SIGN_FACTOR =  1.0; flag |= CF_LOW; }
        else           { Q_CUT = QMAX; SIGN_FACTOR = -1.0; }
        const Real VC1 = DY*DZ;                                              // VC1(-3:3), uniform grid (:801-807)
        const Real VC  = DX*VC1;                                             // VC(-3:3) (:822-828)
        const Real MASS_C = amrex::Math::abs(Q_CUT - QC)*VC;                 // :830
        Real MASS_N[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};               // MASS_N(-3:3) -> [d+3]
```
K2:
```fortran
integer :: i, j, k, d, flag, i0, i1, j0, j1, k0, k1
real(eb) :: QMAX, QC, Q_CUT, SIGN_FACTOR, VC1, VC, MASS_C, SUM_MASS_N, CONST, MASS_N(-3:3)
i0 = lo(1); i1 = hi(1); j0 = lo(2); j1 = hi(2); k0 = lo(3); k1 = hi(3)
#if defined(__NVCOMPILER)
! nvfortran 26.9 rejects has_device_addr (syntax error); ADR-001 fallback is_device_ptr on the
! explicit-shape dummies (OpenMP 5.1: non-c_ptr is_device_ptr items are treated as has_device_addr)
!$omp target teams loop collapse(3) is_device_ptr(Q, RHOP, MASK, T, CF) &
#else
!$omp target teams loop collapse(3) has_device_addr(Q, RHOP, MASK, T, CF) &
#endif
!$omp private(d, flag, QMAX, QC, Q_CUT, SIGN_FACTOR, VC1, VC, MASS_C, SUM_MASS_N, CONST, MASS_N)
do k = k0, k1
   do j = j0, j1
      do i = i0, i1
         do d = 0, 6
            T(i,j,k,d) = 0._eb
         enddo
         CF(i,j,k) = 0
         ! S4b: no private copy of a loop-invariant bound (was QMIN = QMIN_IN; QMAX = QMAX_IN at the loop top): with
         ! private(QMIN) nvfortran 26.9 hoists the invariant store to the teams level (-Minfo=mp "implicit private(qmin)"
         ! at the target line) and the device threads use a wrong QMIN: every cell clipped on the GPU (isolated with
         ! prototypes/s4_cuda_mass/k2_repro). QMIN_IN is read directly; QMAX depends on the cell, so it cannot be
         ! hoisted. Same values, same comparisons as before.
         QMAX = merge(RHOP(i,j,k), QMAX_IN, SPECIES /= 0)                      ! RHO_ZZ_MAX (:887)
         QC = Q(i,j,k,NQI)
         if (SPECIES /= 0) then                                                ! :885, :888
            if (MASK(i,j,k,0) /= 0) cycle
            if (QC >= QMIN_IN .and. QC <= QMAX) cycle
         else                                                                  ! :809, :811
            if (QC >= QMIN_IN .and. QC <= QMAX) cycle
            if (MASK(i,j,k,0) /= 0) cycle
         endif
         flag = CF_CLIPPED
         if (QC < QMIN_IN) then
            Q_CUT = QMIN_IN; SIGN_FACTOR = 1._eb; flag = ior(flag, CF_LOW)
         else
            Q_CUT = QMAX; SIGN_FACTOR = -1._eb
         endif
         VC1 = DY*DZ
```

## 4. Questions (answer each in one or two lines)
Please answer independently, then send the answers to the Chief Architect. "Reads" means: could you follow what the code does and find a mistake in it without help.

1. **Reading.** For §3.2 and §3.3, which one can you read fastest, and which one do you trust more after reading it? (K1 / K2 / equal.)
2. **Matching the FDS source.** Which is easier to compare line by line with §3.1 to check it is the same physics?
3. **Editing.** You must add one term to `RHS`. In which style would you be more confident of doing it correctly without a C++ or GPU expert? Name what you would worry about in each.
4. **Hidden rules.** K2 has two rules that the compiler does not enforce and the CPU tests cannot catch (§2, failures 1 and 2). Is that acceptable for a Fortran team that will add kernels for years? K1 has none, but uses C++ syntax. Which risk is smaller for the project?
5. **Plumbing.** In K2 every array comes with a bounds vector (§3.3 argument list, §3.3 wrapper). Does this hurt readability enough to matter, or do you stop seeing it after one kernel?
6. **Directives.** Is the single `!$omp target teams loop` line, with `private` and `is_device_ptr`, something you could maintain, given the FDS code already uses `!$OMP DO`?
7. **Toolchain.** K2 depends on one vendor compiler (`nvfortran`) that had two defects in one kernel and is not officially tested with AMReX. K1 uses the compiler AMReX supports. How much weight do you give this against readability? (none / some / decisive)
8. **Decision.** Which style should the physics kernels use: K1, K2, or K2 with CUDA Fortran allowed for profiled hot loops (each keeping its OpenMP version)? One line of reasons.

## 5. What happens next
The Chief Architect records the answers and the reasons in ADR-001 ("P1 readability review record"), the owner signs off the outcome, and ADR-001 is set to Accepted with the rejected style listed under "Rejected alternatives". If the answers split, the owner decides. GPU timing (stream ordering, sync cost, a profile) is an implementation item and not part of this review; it can only reopen the kernel style if K2 is picked and the measured sync cost is unacceptable, with K1 as the fallback.
