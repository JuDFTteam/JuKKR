C ************************************************************************
      SUBROUTINE CINIT(N,A)
C ************************************************************************
c-----------------------------------------------------------------------
c     initialize the first n values of a complex array a with zero
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX A(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
      DOUBLE COMPLEX CZERO
      parameter(CZERO=(0.0d0,0.0d0))
C     ..
      DO 10 I = 1,N
        A(I) = CZERO
   10 CONTINUE
      END
