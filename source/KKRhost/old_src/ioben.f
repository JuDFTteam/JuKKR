      INTEGER FUNCTION IOBEN(R)
c-----------------------------------------------------------------------
c
c                             --   --
c     Calculates the function |  r  |  (next upper or equal integer)
c                             |     |
c
c     Descrition of input parameters:
c
c       r : real number to look for
c
c                                           Rudolf Berrendorf, July 1992
c                                           last update: February 1994
c-----------------------------------------------------------------------


C     .. Scalar Arguments ..

      DOUBLE PRECISION R
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,NINT
C     ..

C     .. Parameters ..
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1D-6)
C     ..
      IF ((NINT(R)-R).LT.EPS) THEN
        IF (ABS(NINT(R)-R).LT.EPS) THEN
          IOBEN = NINT(R)
        ELSE
          IOBEN = NINT(R+1.0)
        END IF
      ELSE
        IF (ABS(INT(R)-R).LT.EPS) THEN
          IOBEN = INT(R)
        ELSE
          IOBEN = INT(R+1.0)
        END IF
      END IF

      RETURN
      END
