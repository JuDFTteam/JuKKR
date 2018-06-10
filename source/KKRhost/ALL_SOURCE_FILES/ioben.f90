INTEGER FUNCTION ioben(r)
!-----------------------------------------------------------------------

!                             --   --
!     Calculates the function |  r  |  (next upper or equal integer)
!                             |     |

!     Descrition of input parameters:

!       r : real number to look for

!                                           Rudolf Berrendorf, July 1992
!                                           last update: February 1994
!-----------------------------------------------------------------------

implicit none

!.. Scalar Arguments ..

      DOUBLE PRECISION R
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,INT,NINT
!..

!.. Parameters ..
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1D-6)
!..

IF ((nint(r)-r) < eps) THEN
  IF (ABS(nint(r)-r) < eps) THEN
    ioben = nint(r)
  ELSE
    ioben = nint(r+1.0)
  END IF
ELSE
  IF (ABS(INT(r)-r) < eps) THEN
    ioben = INT(r)
  ELSE
    ioben = INT(r+1.0)
  END IF
END IF

RETURN
END FUNCTION ioben
