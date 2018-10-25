      SUBROUTINE RINIT(N,A)
C **********************************************************************
C * Setting the first N values of a double precision array A to zero   *
C **********************************************************************
C     ..
C     .. Arguments ..
      INTEGER N
      DOUBLE PRECISION A(*)
C     ..
C     .. Locals ..
      INTEGER I,M,MP1
      DOUBLE PRECISION DZERO
C     ..
      DATA DZERO / 0.0D0 /
C     ..
C     ..
      M = MOD(N,5)
      IF ( M.NE.0 ) THEN
         DO I = 1,M
            A(I) = DZERO
         END DO
         IF ( N.LT.5 ) RETURN
      END IF
      MP1 = M + 1
      DO I = MP1,N,5
        A(I  ) = DZERO
        A(I+1) = DZERO
        A(I+2) = DZERO
        A(I+3) = DZERO
        A(I+4) = DZERO
      END DO
C
      END
