! **********************************************************************
SUBROUTINE btom(pl1,pl2,BLOCK,nsize,gin,almd,lsub)
!     This subroutine copies or subtracts a block to a matrix
! **********************************************************************
      IMPLICIT NONE
!.. Scalar Arguments ..
      INTEGER ALMD,NSIZE,PL1,PL2
      LOGICAL LSUB
!..
!.. Array Arguments ..
      DOUBLE COMPLEX BLOCK(NSIZE,NSIZE),GIN(ALMD,ALMD)
!..
!.. Local Scalars ..
      INTEGER I1,I1S,I2,I2S
!     ..
i1s = (pl1-1)*nsize
i2s = (pl2-1)*nsize
IF (lsub) THEN
  DO i1 = 1,nsize
    DO i2 = 1,nsize
      gin(i1s+i1,i2s+i2) = gin(i1s+i1,i2s+i2) - BLOCK(i1,i2)
    END DO
  END DO
ELSE
  DO i1 = 1,nsize
    DO i2 = 1,nsize
      gin(i1s+i1,i2s+i2) = BLOCK(i1,i2)
    END DO
  END DO
END IF

RETURN

END SUBROUTINE btom
