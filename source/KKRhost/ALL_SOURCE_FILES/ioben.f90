    Integer Function ioben(r)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------

!                             --   --
!     Calculates the function |  r  |  (next upper or equal integer)
!                             |     |

!     Descrition of input parameters:

!       r : real number to look for

!                                           Rudolf Berrendorf, July 1992
!                                           last update: February 1994
!-----------------------------------------------------------------------

      Implicit None

!.. Scalar Arguments ..

      Real (Kind=dp) :: r
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, int, nint
!..

!.. Parameters ..
      Real (Kind=dp) :: eps
      Parameter (eps=1E-6_dp)
!..

      If ((nint(r)-r)<eps) Then
        If (abs(nint(r)-r)<eps) Then
          ioben = nint(r)
        Else
          ioben = nint(r+1.0_dp)
        End If
      Else
        If (abs(int(r)-r)<eps) Then
          ioben = int(r)
        Else
          ioben = int(r+1.0_dp)
        End If
      End If

      Return
    End Function
