      SUBROUTINE GETBR3(
     >                  NEMBD1,NLBASIS,ALAT,TLEFT,NRBASIS,TRIGHT,
     X                  BRAVAIS,RECBV,VOLUME0)
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

      TWOPI = 8.D0*DATAN(1.D0)

      ZMAX = 1.D100  ! Large initial value
      DO I1 = 1,NLBASIS
         Z1 = TLEFT(3,I1)
         DO I2 = 1,NRBASIS
            Z2 = TRIGHT(3,I2)
            DIFF = ABS(Z1-Z2)
            ! Smallest distance between left and right sites:
            IF (DIFF.LT.ZMAX) ZMAX = DIFF 
         ENDDO
      ENDDO

      BRAVAIS(3,1) = 0.D0
      BRAVAIS(3,2) = 0.D0
      BRAVAIS(1:2,3) = 0.D0
      BRAVAIS(3,3) = ZMAX     ! Units of a

      VOLUME0 = VOLUME0 * ZMAX * ALAT

      RECBV(3,1) = 0.D0
      RECBV(3,2) = 0.D0
      RECBV(1:2,3) = 0.D0
      RECBV(3,3) = 1.D0 / ZMAX ! Units of 2pi/a

      RETURN
      END
