!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_getbr3
  
  private
  public :: getbr3

contains

  !-------------------------------------------------------------------------------
  !> Summary: Construct auxiliary third Bravais vector for slab calculation
  !> Author: Phivos Mavropoulos
  !> Date: Oct. 2016
  !> Category: KKRhost, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> In case of slab geometry, define a 3rd Bravais vector along
  !> the z direction so that a periodic slab is constructed.
  !> Also set the z-component of the two first BRAVAIS to zero.
  !> Use the "outside region" defined by tleft, tright to make sure
  !> that the bravais vector exceeds the physically relevant thickness.
  !> Purpose: to use the 3D-Madelung routine in slab systems.
  !> Assumes cartesian coordinates in z direction.
  !-------------------------------------------------------------------------------
  subroutine getbr3(nembd1, nlbasis, alat, tleft, nrbasis, tright, bravais, recbv, volume0)
    use :: mod_datatypes, only: dp
    implicit none

    integer :: nembd1, nlbasis, nrbasis
    real (kind=dp) :: alat, tleft(3, nembd1), tright(3, nembd1)
    real (kind=dp) :: bravais(3, 3), recbv(3, 3), volume0
    integer :: i1, i2
    real (kind=dp) :: zmax, z1, z2, diff
    

    ! Large initial value
    zmax = 1.e100_dp
    do i1 = 1, nlbasis
      z1 = tleft(3, i1)
      do i2 = 1, nrbasis
        z2 = tright(3, i2)
        diff = abs(z1-z2)
        ! Smallest distance between left and right sites:
        if (diff<zmax) zmax = diff
      end do
    end do

    bravais(3, 1) = 0.e0_dp
    bravais(3, 2) = 0.e0_dp
    bravais(1:2, 3) = 0.e0_dp
    bravais(3, 3) = zmax     ! Units of a

    volume0 = volume0*zmax*alat

    recbv(3, 1) = 0.e0_dp
    recbv(3, 2) = 0.e0_dp
    recbv(1:2, 3) = 0.e0_dp
    recbv(3, 3) = 1.e0_dp/zmax     ! Units of 2pi/a

    return
  end subroutine getbr3

end module mod_getbr3
