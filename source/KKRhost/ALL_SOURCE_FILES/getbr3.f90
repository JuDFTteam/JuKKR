subroutine getbr3(nembd1, nlbasis, alat, tleft, nrbasis, tright, bravais, &
  recbv, volume0)
!     In case of slab geometry, define a 3rd Bravais vector along
!     the z direction so that a periodic slab is constructed.
!     Also set the z-component of the two first BRAVAIS to zero.
!     Use the "outside region" defined by tleft, tright to make sure
!     that the bravais vector exceeds the physically relevant thickness.
!     Purpose: to use the 3D-Madelung routine in slab systems.
!     Assumes cartesian coordinates in z direction.
!     Phivos Mavropoulos, Oct. 2016


  implicit none
  include 'inc.p'

  integer :: nembd1, nlbasis, nrbasis
  real *8 :: alat, tleft(3, nembd1), tright(3, nembd1)
! Large initial value
  real *8 :: bravais(3, 3), recbv(3, 3), volume0
! Smallest distance between left and right sites:
  integer :: i1, i2
  real *8 :: zmax, z1, z2, diff, twopi

  twopi = 8.d0*datan(1.d0)
! Units of a
  zmax = 1.d100 
  do i1 = 1, nlbasis
    z1 = tleft(3, i1)
    do i2 = 1, nrbasis
      z2 = tright(3, i2)
      diff = abs(z1-z2)

      if (diff<zmax) zmax = diff
    end do
  end do
! Units of 2pi/a
  bravais(3, 1) = 0.d0
  bravais(3, 2) = 0.d0
  bravais(1:2, 3) = 0.d0
  bravais(3, 3) = zmax 
!     In case of slab geometry, define a 3rd Bravais vector along
  volume0 = volume0*zmax*alat
!     the z direction so that a periodic slab is constructed.
  recbv(3, 1) = 0.d0
  recbv(3, 2) = 0.d0
  recbv(1:2, 3) = 0.d0
  recbv(3, 3) = 1.d0/zmax !     Also set the z-component of the two first BRAVAIS to zero.
!     Use the "outside region" defined by tleft, tright to make sure
  return
end subroutine
