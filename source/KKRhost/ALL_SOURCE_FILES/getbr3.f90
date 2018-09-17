module mod_getbr3
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine getbr3(nembd1, nlbasis, alat, tleft, nrbasis, tright, bravais, recbv, volume0)
    ! In case of slab geometry, define a 3rd Bravais vector along
    ! the z direction so that a periodic slab is constructed.
    ! Also set the z-component of the two first BRAVAIS to zero.
    ! Use the "outside region" defined by tleft, tright to make sure
    ! that the bravais vector exceeds the physically relevant thickness.
    ! Purpose: to use the 3D-Madelung routine in slab systems.
    ! Assumes cartesian coordinates in z direction.
    ! Phivos Mavropoulos, Oct. 2016


    implicit none

    integer :: nembd1, nlbasis, nrbasis
    real (kind=dp) :: alat, tleft(3, nembd1), tright(3, nembd1)
    ! Large initial value
    real (kind=dp) :: bravais(3, 3), recbv(3, 3), volume0
    ! Smallest distance between left and right sites:
    integer :: i1, i2
    real (kind=dp) :: zmax, z1, z2, diff

    ! Units of a
    zmax = 1.e100_dp
    do i1 = 1, nlbasis
      z1 = tleft(3, i1)
      do i2 = 1, nrbasis
        z2 = tright(3, i2)
        diff = abs(z1-z2)

        if (diff<zmax) zmax = diff
      end do
    end do
    ! Units of 2pi/a
    bravais(3, 1) = 0.e0_dp
    bravais(3, 2) = 0.e0_dp
    bravais(1:2, 3) = 0.e0_dp
    bravais(3, 3) = zmax
    ! In case of slab geometry, define a 3rd Bravais vector along
    volume0 = volume0*zmax*alat
    ! the z direction so that a periodic slab is constructed.
    recbv(3, 1) = 0.e0_dp
    recbv(3, 2) = 0.e0_dp
    recbv(1:2, 3) = 0.e0_dp
    recbv(3, 3) = 1.e0_dp/zmax     ! Also set the z-component of the two first
    ! BRAVAIS to zero.
    ! Use the "outside region" defined by tleft, tright to make sure
    return
  end subroutine getbr3

end module mod_getbr3
