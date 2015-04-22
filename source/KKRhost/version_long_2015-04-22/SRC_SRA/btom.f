C **********************************************************************
      SUBROUTINE BTOM(PL1,PL2,BLOCK,NSIZE,GIN,ALMD,LSUB)
C     This subroutine copies or subtracts a block to a matrix
C **********************************************************************
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER ALMD,NSIZE,PL1,PL2
      LOGICAL LSUB
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX BLOCK(NSIZE,NSIZE),GIN(ALMD,ALMD)
C     ..
C     .. Local Scalars ..
      INTEGER I1,I1S,I2,I2S
C     ..
      I1S = (PL1-1)*NSIZE
      I2S = (PL2-1)*NSIZE
      IF (LSUB) THEN
        DO I1 = 1,NSIZE
          DO I2 = 1,NSIZE
            GIN(I1S+I1,I2S+I2) = GIN(I1S+I1,I2S+I2) - BLOCK(I1,I2)
          END DO
        END DO
      ELSE
        DO I1 = 1,NSIZE
          DO I2 = 1,NSIZE
            GIN(I1S+I1,I2S+I2) = BLOCK(I1,I2)
          END DO
        END DO
      END IF

      RETURN

      END
