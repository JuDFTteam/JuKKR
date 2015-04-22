      DOUBLE PRECISION FUNCTION SSUM(N,V,IV)
c **********************************************************************
c        sum up the first N elements of the double precision
c        array V(*) with a stepwidth of IV
c ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = 0.0D0
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      SSUM = VSUM
      END
