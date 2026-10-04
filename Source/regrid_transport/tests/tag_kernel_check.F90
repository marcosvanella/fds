! tag_kernel_check.F90: standalone check of the tagging kernels (rt_tag_kernels.F90), built twice by CMake: the host build (default) and the
! offload-source build (RT_OFFLOAD: `target teams distribute parallel do`; where no accelerator is present the target regions run on the host).
! Each case is compared with an independent cell-by-cell reference written here; a result line (count and checksum of the tag bytes) is printed
! for every case so that the two builds can be compared bitwise by tests/compare_tag_kernel_outputs.cmake. Arrays are allocatable so that a build with
! managed memory (nvfortran -gpu=mem:managed) hands the kernels pointers valid on the device.
program tag_kernel_check
use, intrinsic :: iso_c_binding, only: c_int, c_double, c_signed_char, c_long
use rt_tag_kernels
implicit none
integer, parameter :: eb = c_double
integer :: nfail, ncase
integer(c_int) :: lo(3), hi(3), qlo(3), qhi(3)
real(eb), allocatable :: q(:,:,:), den(:,:,:)
integer(c_int), allocatable :: cov(:,:,:)
integer(c_signed_char), allocatable :: tag(:,:,:), ref(:,:,:), tag2(:,:,:)
nfail = 0
ncase = 0

call run_geometry([1, 1, 1], [17, 13, 11], [1, 1, 1], 11_c_int)    ! 3-D
call run_geometry([-3, 0, 2], [20, 0, 17], [1, 0, 1], 7_c_int)     ! 2-D: one cell in y, no y ghost layer, negative low index
call run_geometry([0, 0, 0], [30, 0, 0], [1, 0, 0], 5_c_int)       ! 1-D
call run_boxes()

if (nfail == 0) then
   write (*, '(a,i0,a)') 'TAGKERNEL PASS ', ncase, ' cases'
else
   write (*, '(a,i0,a,i0,a)') 'TAGKERNEL FAIL ', nfail, ' of ', ncase, ' cases'
   stop 1
end if

contains

subroutine lcg_fill(a, seed)
real(eb), intent(out) :: a(:,:,:)
integer(c_int), intent(in) :: seed
integer(c_long) :: s
integer :: i, j, k
s = int(seed, c_long)*2654435761_c_long + 12345_c_long
do k = 1, size(a,3)
do j = 1, size(a,2)
do i = 1, size(a,1)
   s = iand(s*6364136223846793005_c_long + 1442695040888963407_c_long, huge(1_c_long))
   a(i,j,k) = real(iand(ishft(s, -20), 1048575_c_long), eb)/1048576._eb
end do
end do
end do
end subroutine lcg_fill

function checksum(t) result(h)
integer(c_signed_char), intent(in) :: t(:,:,:)
integer(c_long) :: h
integer :: i, j, k
h = 1469598103_c_long
do k = 1, size(t,3)
do j = 1, size(t,2)
do i = 1, size(t,1)
   h = iand(ieor(h, int(t(i,j,k), c_long)+1_c_long)*1099511_c_long + 7_c_long, 2147483647_c_long)
end do
end do
end do
end function checksum

! independent reference: the criterion written out per direction, with the unit step (no helper functions shared with the kernel)
subroutine reference(mode, thr, base, keepfac, use_den, use_cov, rel, dirs, r)
integer(c_int), intent(in) :: mode, use_den, use_cov, rel, dirs(3)
real(eb), intent(in) :: thr, base, keepfac
integer(c_signed_char), intent(inout) :: r(:,:,:)
integer :: i, j, k, d, s
integer :: o(3)
real(eb) :: te, v0, vn, df, best
do k = lo(3), hi(3)
do j = lo(2), hi(2)
do i = lo(1), hi(1)
   te = thr
   if (use_cov /= 0 .and. cov(i,j,k) /= 0) te = thr*keepfac
   v0 = q(i,j,k); if (use_den /= 0) v0 = q(i,j,k)/den(i,j,k)
   best = 0._eb
   if (mode == 0) then
      if (v0 - base > te) r(i-lo(1)+1, j-lo(2)+1, k-lo(3)+1) = 2_c_signed_char
   else
      do d = 1, 3
         if (dirs(d) == 0) cycle
         do s = -1, 1, 2
            o = [i, j, k]; o(d) = o(d) + s
            vn = q(o(1),o(2),o(3)); if (use_den /= 0) vn = q(o(1),o(2),o(3))/den(o(1),o(2),o(3))
            df = abs(v0 - vn)
            if (rel /= 0) then
               if (max(abs(v0), abs(vn)) > 0._eb) then
                  df = df/max(abs(v0), abs(vn))
               else
                  df = 0._eb
               end if
            end if
            best = max(best, df)
         end do
      end do
      if (best > te) r(i-lo(1)+1, j-lo(2)+1, k-lo(3)+1) = 2_c_signed_char
   end if
end do
end do
end do
end subroutine reference

subroutine run_geometry(blo, bhi, dirs3, seed)
integer, intent(in) :: blo(3), bhi(3), dirs3(3)
integer(c_int), intent(in) :: seed
integer(c_int) :: dirs(3), mode, use_den, use_cov, rel
integer :: g(3), n(3), nset, ncell
real(eb) :: thr
integer(c_long) :: cnt, cs
lo = int(blo, c_int); hi = int(bhi, c_int)
g = merge(1, 0, dirs3 /= 0)
qlo = lo - int(g, c_int); qhi = hi + int(g, c_int)
dirs = int(dirs3, c_int)
if (allocated(q)) deallocate (q, den, cov, tag, ref, tag2)
allocate (q(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3)), den(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3)))
allocate (cov(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3)))
allocate (tag(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3)), tag2(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3)))
allocate (ref(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
n = int(hi - lo + 1)
ncell = product(n)
call lcg_fill(q, seed)
call lcg_fill(den, seed + 100_c_int)
den = den + 0.5_eb                 ! positive
q(qlo(1):qhi(1), qlo(2):qhi(2), qlo(3):qhi(3)) = q*den
cov = 0_c_int
cov(lo(1):hi(1):2, :, :) = 1_c_int
do mode = 0, 1
do use_den = 0, 1
do use_cov = 0, 1
do rel = 0, 1
   if (mode == 0 .and. rel == 1) cycle
   thr = merge(0.6_eb, 0.3_eb, mode == 0)
   if (use_den == 0 .and. mode == 0) thr = 0.5_eb
   if (rel == 1) thr = 0.8_eb
   tag = 0_c_signed_char
   tag(lo(1), lo(2), lo(3)) = 2_c_signed_char      ! a tag set before the call must survive
   ref = 0_c_signed_char
   ref(lo(1), lo(2), lo(3)) = 2_c_signed_char
   ncase = ncase + 1
   call rt_tag_cells(lo, hi, q, qlo, qhi, den, use_den, cov, qlo, qhi, use_cov, tag, qlo, qhi, mode, thr, 0.1_eb, 0.5_eb, rel, dirs)
   call reference(mode, thr, 0.1_eb, 0.5_eb, use_den, use_cov, rel, dirs, ref)
   call rt_tag_count(lo, hi, tag, qlo, qhi, cnt)
   cs = checksum(tag(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
   nset = count(ref /= 0_c_signed_char)
   write (*, '(a,i0,a,i0,i0,i0,i0,a,i0,a,i0)') 'TK geom=', seed, ' mode/den/cov/rel=', mode, use_den, use_cov, rel, ' count=', cnt, ' sum=', cs
   if (any(tag(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) /= ref)) then
      nfail = nfail + 1
      write (*, '(a)') 'TAGKERNEL mismatch with the independent reference'
   end if
   if (cnt /= nset .or. cnt == 0 .or. cnt == ncell) then
      nfail = nfail + 1
      write (*, '(a,i0,a,i0,a,i0)') 'TAGKERNEL count ', cnt, ' reference ', nset, ' cells ', ncell      ! all-set or none-set cases prove nothing
   end if
   ! negative control: a reference with a different threshold must differ from the kernel result (the comparison can fail)
   ref = 0_c_signed_char
   ref(lo(1), lo(2), lo(3)) = 2_c_signed_char
   call reference(mode, thr*1.5_eb, 0.1_eb, 0.5_eb, use_den, use_cov, rel, dirs, ref)
   if (all(tag(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) == ref)) then
      nfail = nfail + 1
      write (*, '(a)') 'TAGKERNEL negative control did not detect a changed threshold'
   end if
end do
end do
end do
end do
end subroutine run_geometry

subroutine run_boxes()
integer(c_int) :: blo(3), bhi(3), clo(3), chi(3)
integer(c_long) :: cnt, cs
integer :: nexp
lo = [0, -2, 1]; hi = [9, 7, 12]
if (allocated(ref)) deallocate (ref)
if (allocated(tag)) deallocate (tag, tag2)
allocate (tag(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), tag2(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
tag = 0_c_signed_char
blo = [2, -1, 3]; bhi = [20, 4, 5]               ! box sticks out of the tile in x: clipped to the tile
ncase = ncase + 1
call rt_tag_box(lo, hi, tag, lo, hi, blo, bhi)
call rt_tag_count(lo, hi, tag, lo, hi, cnt)
nexp = (9 - 2 + 1)*(4 - (-1) + 1)*(5 - 3 + 1)
cs = checksum(tag)
write (*, '(a,i0,a,i0)') 'TK box count=', cnt, ' sum=', cs
if (cnt /= nexp) then
   nfail = nfail + 1
   write (*, '(a,i0,a,i0)') 'TAGKERNEL rt_tag_box count ', cnt, ' expected ', nexp
end if
tag2 = 0_c_signed_char
clo = [4, 0, 2]; chi = [6, 3, 4]
ncase = ncase + 1
call rt_tag_copy_box(lo, hi, tag, lo, hi, tag2, lo, hi, clo, chi)     ! intersection of the two boxes only
call rt_tag_count(lo, hi, tag2, lo, hi, cnt)
cs = checksum(tag2)
nexp = (6 - 4 + 1)*(3 - 0 + 1)*(4 - 3 + 1)                   ! tag box has z in 3..5, copy box z in 2..4: z 3..4
write (*, '(a,i0,a,i0)') 'TK copy count=', cnt, ' sum=', cs
if (cnt /= nexp) then
   nfail = nfail + 1
   write (*, '(a,i0,a,i0)') 'TAGKERNEL rt_tag_copy_box count ', cnt, ' expected ', nexp
end if
end subroutine run_boxes

end program tag_kernel_check
