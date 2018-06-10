! **********************************************************************
SUBROUTINE bofm(pl1,pl2,BLOCK,nsize,gin,almd)
! **********************************************************************

IMPLICIT NONE
!.. Scalar Arguments ..
INTEGER ALMD,NSIZE,PL1,PL2
!..
!.. Array Arguments ..
DOUBLE COMPLEX BLOCK(NSIZE,NSIZE),GIN(ALMD,ALMD)
!..
!.. Local Scalars ..
INTEGER I1,I1S,I2,I2S
!..
i1s = (pl1-1)*nsize
i2s = (pl2-1)*nsize
DO i1 = 1,nsize
  DO i2 = 1,nsize
    BLOCK(i1,i2) = gin(i1s+i1,i2s+i2)
  END DO
END DO

RETURN

END SUBROUTINE bofm
