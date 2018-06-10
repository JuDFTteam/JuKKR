SUBROUTINE potcut(imt1,irc1,ins,lmpot,r,vm2z,vspsme,vins,z1,  &
        irmd,irmind)
! **********************************************************************
! * set potential equal zero between muffin-tin and outer sphere       *
! **********************************************************************
IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
DOUBLE PRECISION :: z1
INTEGER :: irmd,irmind
INTEGER :: imt1,ins,irc1,lmpot
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION :: r(*),vins(irmind:irmd,*),vm2z(*),vspsme(*)
!     ..
!     .. Local Scalars ..
INTEGER :: ir,ist,lm
!     ..
!     .. Intrinsic Functions ..
INTRINSIC MAX
!     ..
WRITE(1337,*) 'potcut: potential equal 2*Z/R between MT ', 'and outer sphere'
DO ir = imt1 + 1,irc1
  vm2z(ir) = 2.0D0*z1/r(ir)
  vspsme(ir) = 2.0D0*z1/r(ir)
END DO

IF (ins >= 1) THEN
  ist = MAX(irmind,imt1+1)
  DO ir = ist,irc1
    DO lm = 2,lmpot
      vins(ir,lm) = 0.0D0
    END DO
  END DO
END IF
END                           ! SUBROUTINE POTCUT
