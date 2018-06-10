integer function ioben(r)
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

  double precision :: r
!..
!.. Intrinsic Functions ..
  intrinsic :: abs, int, nint
!..

!.. Parameters ..
  double precision :: eps
  parameter (eps=1d-6)
!..

  if ((nint(r)-r)<eps) then
    if (abs(nint(r)-r)<eps) then
      ioben = nint(r)
    else
      ioben = nint(r+1.0)
    end if
  else
    if (abs(int(r)-r)<eps) then
      ioben = int(r)
    else
      ioben = int(r+1.0)
    end if
  end if

  return
end function
