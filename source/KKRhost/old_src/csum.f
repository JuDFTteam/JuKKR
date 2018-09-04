c 19.10.95 *************************************************************
      COMPLEX*16 FUNCTION CSUM(N,V,IV)
c **********************************************************************
c        sum up the first N elements of the double complex
c        array V(*) with a stepwidth of IV
c ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = (0D0,0D0)
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      CSUM = VSUM
      RETURN
      END
