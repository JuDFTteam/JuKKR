subroutine potcut(imt1, irc1, ins, lmpot, r, vm2z, vspsme, vins, z1, irmd, &
  irmind)
! **********************************************************************
! * set potential equal zero between muffin-tin and outer sphere       *
! **********************************************************************
  implicit none
!     ..
!     .. Scalar Arguments ..
  double precision :: z1
  integer :: irmd, irmind
  integer :: imt1, ins, irc1, lmpot
!     ..
!     .. Array Arguments ..
  double precision :: r(*), vins(irmind:irmd, *), vm2z(*), vspsme(*)
!     ..
!     .. Local Scalars ..
  integer :: ir, ist, lm
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: max
!     ..
  write (1337, *) 'potcut: potential equal 2*Z/R between MT ', &
    'and outer sphere'
  do ir = imt1 + 1, irc1
    vm2z(ir) = 2.0d0*z1/r(ir)
    vspsme(ir) = 2.0d0*z1/r(ir)
  end do

  if (ins>=1) then
    ist = max(irmind, imt1+1)
    do ir = ist, irc1
      do lm = 2, lmpot
        vins(ir, lm) = 0.0d0
      end do
    end do
  end if
end subroutine ! SUBROUTINE POTCUT
