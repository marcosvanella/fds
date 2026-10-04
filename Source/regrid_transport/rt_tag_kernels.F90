! rt_tag_kernels.F90: tagging kernels of the regrid_transport library (R4), K2 style (ADR-001): Fortran, flat explicit-shape arrays in the AMReX
! index space, no module variables, no I/O, no AMReX types. Directive pattern of the project's K2 kernels (s4_omp.inc / s5gen_k2.F90, docs/tools/README.md): only the
! macro tail is a macro, the `!$omp` sentinel stays literal.
!   host build    : `parallel do collapse(2)` (a comment when OpenMP is off).
!   offload build : RT_OFFLOAD or S4_OFFLOAD defined, or nvfortran -mp=gpu. A region without a callee is `target teams loop collapse(3)` (S4_LOOP); a region that calls a
!                   declare-target routine is `target teams distribute parallel do collapse(3)` (S5_LOOP_CALLEE, rule 7 of ADR-001 v0.7). The array dummies are
!                   device addresses (is_device_ptr on nvfortran, has_device_addr elsewhere); no map of field data (the AMReX arena owns it); scalars are firstprivate.
!   Switches (same as the generated kernels): S4K2_COLL (collapse depth, default 3), S4K2_THREADS / S4K2_TEAMS (thread_limit / num_teams), S5_FORCE_DPD.
! The CMake option RT_TAG_OFFLOAD builds the library with the offload variant (same entry points and arguments). Checked by tests/tag_kernel_check.F90 (host default vs offload source, bitwise;
! the 38 cases ran identically on an RTX 4070 with nvfortran -mp=gpu -gpu=mem:managed).
!
! Upstream provenance (K2 header note): none. FDS has no refinement tagging; these are new kernels of the AMR code (FR-011, D-058), so there is no upstream file and line to name.
! Passive scalars: the criterion reads one cell field; a species mass fraction is passed as rho*Z with rho (use_den /= 0) and is divided inside, any passive scalar is just another field.
! Cylindrical terms: none; the difference criterion is undivided (no cell sizes) and the kernels see only index space, so Cartesian and cylindrical meshes give the same tags.
!
! A tag is one byte (AMReX TagBox: CLEAR = 0, BUF = 1, SET = 2; only SET cells seed the tag buffer, so the value written here must be SET; TagOps.cpp static_asserts it). The kernels only SET tags (logical OR of the criteria).
!   rt_tag_cells    : threshold or undivided-difference criterion on one cell field, with TAG_KEEP hysteresis from a "covered by the finer level" mask
!   rt_tag_box      : sets all tags inside an index box (user boxes, static finer &MESH footprints)
!   rt_tag_copy_box : dst = src inside an index box (used to clip tags to the refinable region)
!   rt_tag_count    : number of set tags in a box (host only: a count is a sum, and the approved clause list has no sum reduction; a device build counts in C++)
module rt_tag_kernels
use, intrinsic :: iso_c_binding, only: c_int, c_double, c_signed_char, c_long
implicit none
private
integer, parameter :: eb = c_double
integer(c_signed_char), parameter :: RT_TAG_SET = 2_c_signed_char   ! amrex::TagBox::SET
public :: rt_tag_cells, rt_tag_box, rt_tag_copy_box, rt_tag_count

#if !defined(S4K2_COLL)
#define S4K2_COLL 3
#endif
#if defined(S4K2_THREADS) && defined(S4K2_TEAMS)
#define S4K2_CL thread_limit(S4K2_THREADS) num_teams(S4K2_TEAMS)
#elif defined(S4K2_THREADS)
#define S4K2_CL thread_limit(S4K2_THREADS)
#elif defined(S4K2_TEAMS)
#define S4K2_CL num_teams(S4K2_TEAMS)
#else
#define S4K2_CL
#endif
#if defined(RT_OFFLOAD) && !defined(S4_OFFLOAD)
#  define S4_OFFLOAD
#endif
#if defined(S4_OFFLOAD) || defined(__NVCOMPILER_OPENMP_GPU)
#  if defined(S5_FORCE_DPD)
#    define S4_LOOP target teams distribute parallel do collapse(S4K2_COLL) S4K2_CL
#  else
#    define S4_LOOP target teams loop collapse(S4K2_COLL) S4K2_CL
#  endif
#  define S5_LOOP_CALLEE target teams distribute parallel do collapse(S4K2_COLL) S4K2_CL
#  if defined(__NVCOMPILER)
#    define S4_DEV(l) is_device_ptr l
#  else
#    define S4_DEV(l) has_device_addr l
#  endif
#else
#  define S4_LOOP parallel do collapse(2)
#  define S5_LOOP_CALLEE parallel do collapse(2)
#  define S4_DEV(l)
#endif

contains

! Cell value used by the criterion: q, or q/den when use_den /= 0 (species mass fraction from rho*Z and rho).
pure real(eb) function cell_value(q, d, use_den) result(v)
!$omp declare target
real(eb), intent(in) :: q, d
integer(c_int), intent(in) :: use_den
if (use_den /= 0) then
   v = q/d
else
   v = q
end if
end function cell_value

! |v0 - vn|, divided by max(|v0|,|vn|) when rel /= 0 (zero if both are zero).
pure real(eb) function ndiff(v0, vn, rel) result(df)
!$omp declare target
real(eb), intent(in) :: v0, vn
integer(c_int), intent(in) :: rel
real(eb) :: dm
df = abs(v0 - vn)
if (rel /= 0) then
   dm = max(abs(v0), abs(vn))
   if (dm > 0._eb) then
      df = df/dm
   else
      df = 0._eb
   end if
end if
end function ndiff

! mode 0 (ABOVE): tag where  v - base > thr_eff.
! mode 1 (DIFF) : tag where  max over the active directions d of |v - v_nb| (nb = the two neighbours in d) > thr_eff; undivided (no division by dx).
!                 rel /= 0: each difference is divided by max(|v|, |v_nb|) (zero if both are zero).
! thr_eff = thr*keepfac where cov(i,j,k) /= 0 (the cell is under the next finer level; keepfac <= 1 is TAG_KEEP), else thr. use_cov = 0: no mask.
! v = q, or q/den when use_den /= 0 (den has the bounds of q; when unused the caller passes q again).
! q (and den) need one valid ghost layer in every active direction; dirs(d) = 0 for a single-cell direction (not read).
subroutine rt_tag_cells(lo, hi, q, qlo, qhi, den, use_den, cov, clo, chi, use_cov, tag, tlo, thi, &
                        mode, thr, base, keepfac, rel, dirs) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), qlo(3), qhi(3), clo(3), chi(3), tlo(3), thi(3)
integer(c_int), intent(in) :: use_den, use_cov, mode, rel, dirs(3)
real(eb), intent(in) :: q(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3))
real(eb), intent(in) :: den(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3))
integer(c_int), intent(in) :: cov(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
integer(c_signed_char), intent(inout) :: tag(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3))
real(eb), intent(in), value :: thr, base, keepfac
integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, dx, dy, dz
real(eb) :: te, v0, best
logical :: hit
ilo = lo(1); ihi = hi(1); jlo = lo(2); jhi = hi(2); klo = lo(3); khi = hi(3)   ! bounds and direction flags as scalars: array dummies in a region would be copied
dx = dirs(1); dy = dirs(2); dz = dirs(3)
!$omp S5_LOOP_CALLEE private(te, v0, best, hit) S4_DEV((q, den, cov, tag))
do k = klo, khi
do j = jlo, jhi
do i = ilo, ihi
   te = thr
   if (use_cov /= 0) then
      if (cov(i,j,k) /= 0) te = thr*keepfac
   end if
   v0 = cell_value(q(i,j,k), den(i,j,k), use_den)
   if (mode == 0) then
      hit = (v0 - base) > te
   else
      best = 0._eb
      if (dx /= 0) then
         best = max(best, ndiff(v0, cell_value(q(i-1,j,k), den(i-1,j,k), use_den), rel))
         best = max(best, ndiff(v0, cell_value(q(i+1,j,k), den(i+1,j,k), use_den), rel))
      end if
      if (dy /= 0) then
         best = max(best, ndiff(v0, cell_value(q(i,j-1,k), den(i,j-1,k), use_den), rel))
         best = max(best, ndiff(v0, cell_value(q(i,j+1,k), den(i,j+1,k), use_den), rel))
      end if
      if (dz /= 0) then
         best = max(best, ndiff(v0, cell_value(q(i,j,k-1), den(i,j,k-1), use_den), rel))
         best = max(best, ndiff(v0, cell_value(q(i,j,k+1), den(i,j,k+1), use_den), rel))
      end if
      hit = best > te
   end if
   if (hit) tag(i,j,k) = RT_TAG_SET
end do
end do
end do
end subroutine rt_tag_cells

subroutine rt_tag_box(lo, hi, tag, tlo, thi, blo, bhi) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), tlo(3), thi(3), blo(3), bhi(3)
integer(c_signed_char), intent(inout) :: tag(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3))
integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
ilo = max(lo(1), blo(1)); ihi = min(hi(1), bhi(1))
jlo = max(lo(2), blo(2)); jhi = min(hi(2), bhi(2))
klo = max(lo(3), blo(3)); khi = min(hi(3), bhi(3))
!$omp S4_LOOP S4_DEV((tag))
do k = klo, khi
do j = jlo, jhi
do i = ilo, ihi
   tag(i,j,k) = RT_TAG_SET
end do
end do
end do
end subroutine rt_tag_box

subroutine rt_tag_copy_box(lo, hi, src, slo, shi, dst, dlo, dhi, blo, bhi) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), slo(3), shi(3), dlo(3), dhi(3), blo(3), bhi(3)
integer(c_signed_char), intent(in) :: src(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
integer(c_signed_char), intent(inout) :: dst(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3))
integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
ilo = max(lo(1), blo(1)); ihi = min(hi(1), bhi(1))
jlo = max(lo(2), blo(2)); jhi = min(hi(2), bhi(2))
klo = max(lo(3), blo(3)); khi = min(hi(3), bhi(3))
!$omp S4_LOOP S4_DEV((src, dst))
do k = klo, khi
do j = jlo, jhi
do i = ilo, ihi
   dst(i,j,k) = src(i,j,k)
end do
end do
end do
end subroutine rt_tag_copy_box

! Host only (no target region): tags are in host memory here. A device build sums with an AMReX reduction in C++.
subroutine rt_tag_count(lo, hi, tag, tlo, thi, n) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), tlo(3), thi(3)
integer(c_signed_char), intent(in) :: tag(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3))
integer(c_long), intent(out) :: n
integer :: i, j, k
integer(c_long) :: s
s = 0_c_long
do k = lo(3), hi(3)
do j = lo(2), hi(2)
do i = lo(1), hi(1)
   if (tag(i,j,k) /= 0_c_signed_char) s = s + 1_c_long
end do
end do
end do
n = s
end subroutine rt_tag_count

end module rt_tag_kernels
