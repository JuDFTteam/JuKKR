C **********************************************************************
      SUBROUTINE BOFM(PL1,PL2,BLOCK,NSIZE,GIN,ALMD)
C **********************************************************************
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER ALMD,NSIZE,PL1,PL2
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX BLOCK(NSIZE,NSIZE),GIN(ALMD,ALMD)
C     ..
C     .. Local Scalars ..
      INTEGER I1,I1S,I2,I2S
C     ..
      I1S = (PL1-1)*NSIZE
      I2S = (PL2-1)*NSIZE
      DO I1 = 1,NSIZE
        DO I2 = 1,NSIZE
          BLOCK(I1,I2) = GIN(I1S+I1,I2S+I2)
        END DO
      END DO

      RETURN

      END
