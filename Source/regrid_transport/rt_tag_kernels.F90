! rt_tag_kernels.F90: tagging kernels of the regrid_transport library (R4), K2 style (ADR-001): Fortran, flat explicit-shape arrays in the AMReX
! index space, no module variables, no I/O, no AMReX types. Host build: `parallel do` (or plain loops without OpenMP). Offload build
! (RT_OFFLOAD defined, or nvfortran -mp=gpu): `target teams distribute parallel do` with the device data given by is_device_ptr / has_device_addr, no map
! of the big arrays (the AMReX arena owns them). The offload build is NOT tested here (no GPU compiler run); the host path is.
!
! A tag is one byte (AMReX TagBox: 0 = clear, 1 = set). The kernels only SET tags (logical OR of the criteria); they never clear one.
!   rt_tag_cells    : threshold or undivided-difference criterion on one cell field, with TAG_KEEP hysteresis from a "covered by the finer level" mask
!   rt_tag_box      : sets all tags inside an index box (user boxes, static finer &MESH footprints)
!   rt_tag_copy_box : dst = src inside an index box (used to clip tags to the refinable region)
!   rt_tag_count    : number of set tags in a box
module rt_tag_kernels
use, intrinsic :: iso_c_binding, only: c_int, c_double, c_signed_char, c_long
implicit none
private
integer, parameter :: eb = c_double
public :: rt_tag_cells, rt_tag_box, rt_tag_copy_box, rt_tag_count

#if defined(RT_OFFLOAD) || defined(__NVCOMPILER_OPENMP_GPU)
#  define RT_LOOP target teams distribute parallel do collapse(3)
#  if defined(__NVCOMPILER)
#    define RT_DEV(l) is_device_ptr l
#  else
#    define RT_DEV(l) has_device_addr l
#  endif
#else
#  define RT_LOOP parallel do collapse(2)
#  define RT_DEV(l)
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
integer :: i, j, k
real(eb) :: te, v0, best
logical :: hit
!$omp RT_LOOP private(i, te, v0, best, hit) RT_DEV((q, den, cov, tag))
do k = lo(3), hi(3)
do j = lo(2), hi(2)
do i = lo(1), hi(1)
   te = thr
   if (use_cov /= 0) then
      if (cov(i,j,k) /= 0) te = thr*keepfac
   end if
   v0 = cell_value(q(i,j,k), den(i,j,k), use_den)
   if (mode == 0) then
      hit = (v0 - base) > te
   else
      best = 0._eb
      if (dirs(1) /= 0) then
         best = max(best, ndiff(v0, cell_value(q(i-1,j,k), den(i-1,j,k), use_den), rel))
         best = max(best, ndiff(v0, cell_value(q(i+1,j,k), den(i+1,j,k), use_den), rel))
      end if
      if (dirs(2) /= 0) then
         best = max(best, ndiff(v0, cell_value(q(i,j-1,k), den(i,j-1,k), use_den), rel))
         best = max(best, ndiff(v0, cell_value(q(i,j+1,k), den(i,j+1,k), use_den), rel))
      end if
      if (dirs(3) /= 0) then
         best = max(best, ndiff(v0, cell_value(q(i,j,k-1), den(i,j,k-1), use_den), rel))
         best = max(best, ndiff(v0, cell_value(q(i,j,k+1), den(i,j,k+1), use_den), rel))
      end if
      hit = best > te
   end if
   if (hit) tag(i,j,k) = 1_c_signed_char
end do
end do
end do
end subroutine rt_tag_cells

subroutine rt_tag_box(lo, hi, tag, tlo, thi, blo, bhi) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), tlo(3), thi(3), blo(3), bhi(3)
integer(c_signed_char), intent(inout) :: tag(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3))
integer :: i, j, k
!$omp RT_LOOP private(i) RT_DEV((tag))
do k = max(lo(3), blo(3)), min(hi(3), bhi(3))
do j = max(lo(2), blo(2)), min(hi(2), bhi(2))
do i = max(lo(1), blo(1)), min(hi(1), bhi(1))
   tag(i,j,k) = 1_c_signed_char
end do
end do
end do
end subroutine rt_tag_box

subroutine rt_tag_copy_box(lo, hi, src, slo, shi, dst, dlo, dhi, blo, bhi) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), slo(3), shi(3), dlo(3), dhi(3), blo(3), bhi(3)
integer(c_signed_char), intent(in) :: src(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3))
integer(c_signed_char), intent(inout) :: dst(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3))
integer :: i, j, k
!$omp RT_LOOP private(i) RT_DEV((src, dst))
do k = max(lo(3), blo(3)), min(hi(3), bhi(3))
do j = max(lo(2), blo(2)), min(hi(2), bhi(2))
do i = max(lo(1), blo(1)), min(hi(1), bhi(1))
   dst(i,j,k) = src(i,j,k)
end do
end do
end do
end subroutine rt_tag_copy_box

subroutine rt_tag_count(lo, hi, tag, tlo, thi, n) bind(c)
integer(c_int), intent(in) :: lo(3), hi(3), tlo(3), thi(3)
integer(c_signed_char), intent(in) :: tag(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3))
integer(c_long), intent(out) :: n
integer :: i, j, k
integer(c_long) :: s
s = 0_c_long
!$omp parallel do collapse(2) private(i) reduction(+:s)
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
