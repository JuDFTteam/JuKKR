subroutine potcut(imt1, irc1, ins, lmpot, r, vm2z, vspsme, vins, z1, irmd, &
  irmind)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! * set potential equal zero between muffin-tin and outer sphere       *
  ! **********************************************************************
  implicit none
  ! ..
  ! .. Scalar Arguments ..
  real (kind=dp) :: z1
  integer :: irmd, irmind
  integer :: imt1, ins, irc1, lmpot
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: r(*), vins(irmind:irmd, *), vm2z(*), vspsme(*)
  ! ..
  ! .. Local Scalars ..
  integer :: ir, ist, lm
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: max
  ! ..
  write (1337, *) 'potcut: potential equal 2*Z/R between MT ', &
    'and outer sphere'
  do ir = imt1 + 1, irc1
    vm2z(ir) = 2.0e0_dp*z1/r(ir)
    vspsme(ir) = 2.0e0_dp*z1/r(ir)
  end do

  if (ins>=1) then
    ist = max(irmind, imt1+1)
    do ir = ist, irc1
      do lm = 2, lmpot
        vins(ir, lm) = 0.0e0_dp
      end do
    end do
  end if
end subroutine potcut              ! SUBROUTINE POTCUT
