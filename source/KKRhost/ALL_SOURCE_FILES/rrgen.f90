module mod_rrgen

contains

! 02.08.95 *************************************************************
subroutine rrgen(bv1, lsurf, rr, nrd)
  ! **********************************************************************
  ! *                                                                    *
  ! * generates a number of real space vectors to construct the          *
  ! * clusters representing the local surrounding of the atoms in        *
  ! * routine CLSGEN99                                                   *
  ! *                                                                    *
  ! **********************************************************************
  use mod_datatypes, only: dp
  use mod_vmul
  use mod_vadd
  use mod_veq
  use mod_scalpr
  use mod_dsort
  implicit none
  ! ..
  ! .. Scalar arguments ..
  logical :: lsurf
  integer :: nrd
  ! ..
  ! .. Array arguments ..
  real (kind=dp) :: bv1(3, 3), rr(3, 0:nrd)
  ! ..
  ! .. Local scalars ..
  real (kind=dp) :: epsshl, r, r1, r2, r3, rmax, rr2, rs
  integer :: i, j, k, n1, n2, n3, pos, iprint
  integer :: nr
  ! ..
  ! .. Local arrays
  real (kind=dp) :: rabs(nrd), rr1(3, nrd), v(3), vx(3), vy(3), vz(3), vx0(3), &
    vy0(3), vz0(3)
  integer :: ind(nrd)
  ! ..
  ! .. Data Statements ..
  data epsshl/1.0e-5_dp/
  ! ..................................................................
  write (1337, '(5X,A,/)') '< RRGEN > : generation of real space mesh RR(NR)'

  iprint = 0

  call scalpr(bv1(1,1), bv1(1,1), r1)
  call scalpr(bv1(1,2), bv1(1,2), r2)
  call scalpr(bv1(1,3), bv1(1,3), r3)
  rmax = 5.e0_dp

  r1 = sqrt(r1)
  r2 = sqrt(r2)
  r3 = sqrt(r3)
  r = 1.5e0_dp*rmax + sqrt(r1*r1+r2*r2+r3*r3) + epsshl
  rs = r*r
  n1 = nint(r/r1)
  n2 = nint(r/r2)
  if (.not. lsurf) n3 = nint(r/r3)

  n1 = min(12, n1)
  n2 = min(12, n2)
  if (.not. lsurf) n3 = min(12, n3)

  n1 = max(2, n1)
  n2 = max(2, n2)
  if (.not. lsurf) n3 = max(2, n3)

  if (lsurf) n3 = 0

  write (1337, 100) r
  write (1337, 110) rs
  if (lsurf) then
    write (1337, 120) n1, n2
  else
    write (1337, 130) n1, n2, n3
  end if

  nr = 0
  rr(1, 0) = 0.0e0_dp
  rr(2, 0) = 0.0e0_dp
  rr(3, 0) = 0.0e0_dp

  call vmul(bv1(1,1), real(-n1-1, kind=dp), vx0(1))
  call vmul(bv1(1,2), real(-n2-1, kind=dp), vy0(1))
  call vmul(bv1(1,3), real(-n3-1, kind=dp), vz0(1))
  call veq(vx0, vx)
  ! **********************************************************************
  do i = -n1, n1
    call vadd(vx, bv1(1,1), vx)
    call veq(vy0, vy)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do j = -n2, n2
      call vadd(vy, bv1(1,2), vy)
      call veq(vz0, vz)
      ! ----------------------------------------------------------------------
      do k = -n3, n3
        call vadd(vz, bv1(1,3), vz)
        call vadd(vx, vy, v)
        call vadd(v, vz, v)
        call scalpr(v, v, rr2)

        if (((rr2<=rs) .or. (abs(i)+abs(j)+abs(k)<=6)) .and. (rr2>epsshl)) &
          then
          nr = nr + 1

          if (nr>nrd) then
            write (6, *) 'Dimension ERROR. Please, change the ', &
              'parameter NRD in inc.p to ', nr, nrd
            stop
          end if

          rr1(1, nr) = v(1)
          rr1(2, nr) = v(2)
          rr1(3, nr) = v(3)
          rabs(nr) = sqrt(rr2)
        end if
      end do
      ! ----------------------------------------------------------------------
    end do
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do

  ! store changed nr dimension
  nrd = nr
  ! **********************************************************************

  write (1337, 140) nr + 1

  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  if (iprint>0) then
    write (1337, 150)
    write (1337, 170) 0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
  end if
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  call dsort(rabs, ind, nr, pos)
  do i = 1, nr
    pos = ind(i)
    rr(1, i) = rr1(1, pos)
    rr(2, i) = rr1(2, pos)
    rr(3, i) = rr1(3, pos)
    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    if (iprint>0) write (1337, 170) i, rr(1, i), rr(2, i), rr(3, i), rabs(pos)
    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  end do

  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  if (iprint>0) write (1337, 160)
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

100 format (10x, 'Radius R        : ', f15.6, ' (ALAT    units)')
110 format (10x, '       R**2     : ', f15.6, ' (ALAT**2 units)')
120 format (10x, 'mesh divisions  : ', 5x, 2i5)
130 format (10x, 'mesh divisions  : ', 3i5)
140 format (10x, 'vectors created : ', i15)
150 format (/, 10x, 60('+'), /, 18x, &
    'generated real-space mesh-points (ALAT units)', /, 10x, 60('+'), /, 13x, &
    'index      x           y           z          distance  ', /, 10x, &
    60('-'))
160 format (10x, 60('+'))
170 format (10x, i6, 3f12.3, f15.4)
end subroutine rrgen

end module mod_rrgen
