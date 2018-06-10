SUBROUTINE getbr3(  &
        nembd1,nlbasis,alat,tleft,nrbasis,tright,  &
        bravais,recbv,volume0)
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
! Input:
      INTEGER NEMBD1,NLBASIS,NRBASIS
      REAL*8 ALAT,TLEFT(3,NEMBD1),TRIGHT(3,NEMBD1)
! Input/Output
      REAL*8 BRAVAIS(3,3),RECBV(3,3),VOLUME0
! Inside:
      INTEGER I1,I2
      REAL*8 ZMAX,Z1,Z2,DIFF,TWOPI

twopi = 8.d0*DATAN(1.d0)

zmax = 1.d100  ! Large initial value
DO i1 = 1,nlbasis
  z1 = tleft(3,i1)
  DO i2 = 1,nrbasis
    z2 = tright(3,i2)
    diff = ABS(z1-z2)
! Smallest distance between left and right sites:
    IF (diff < zmax) zmax = diff
  END DO
END DO

bravais(3,1) = 0.d0
bravais(3,2) = 0.d0
bravais(1:2,3) = 0.d0
bravais(3,3) = zmax     ! Units of a

volume0 = volume0 * zmax * alat

recbv(3,1) = 0.d0
recbv(3,2) = 0.d0
recbv(1:2,3) = 0.d0
recbv(3,3) = 1.d0 / zmax ! Units of 2pi/a

RETURN
END SUBROUTINE getbr3
