!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Set potential equal zero between muffin-tin and outer sphere
!> Author:
!> Set potential equal zero between muffin-tin and outer sphere
!------------------------------------------------------------------------------------
module mod_potcut
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Set potential equal zero between muffin-tin and outer sphere
  !> Author: 
  !> Category: potential, KKRhost
  !> Deprecated: False 
  !> Set potential equal zero between muffin-tin and outer sphere
  !-------------------------------------------------------------------------------
  subroutine potcut(imt1, irc1, ins, lmpot, r, vm2z, vspsme, vins, z1, irmd, irmind)

    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in) :: ins      !! (LPOT+1)**2
    integer, intent(in) :: imt1     !! R point at MT radius
    integer, intent(in) :: irc1     !! R point for potential cutting
    integer, intent(in) :: irmd     !! Maximum number of radial points
    integer, intent(in) :: lmpot    !! (LPOT+1)**2
    integer, intent(in) :: irmind   !! irmd - irnsd
    real (kind=dp), intent(in) :: z1 !! Nuclear charge
    real (kind=dp), dimension(*), intent(in) :: r  !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension(*), intent(inout) :: vm2z
    real (kind=dp), dimension(*), intent(inout) :: vspsme
    real (kind=dp), dimension(irmind:irmd, *), intent(inout) :: vins   !! Non-spherical part of the potential
    ! ..
    ! .. Local Scalars ..
    integer :: ir, ist, lm
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: max
    ! ..
    write (1337, *) 'potcut: potential equal 2*Z/R between MT ', 'and outer sphere'
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
  end subroutine potcut            ! SUBROUTINE POTCUT

end module mod_potcut
