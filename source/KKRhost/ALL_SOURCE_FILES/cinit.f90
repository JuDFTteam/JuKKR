    Subroutine cinit(n, a)
      Use mod_datatypes, Only: dp
! **********************************************************************
! * Setting the first N values of a double complex array A to zero     *
! **********************************************************************
!     ..
!     .. Arguments ..
      Integer :: n
      Complex (Kind=dp) :: a(*)
!     ..
!     .. Locals
      Integer :: i, m, mp1
      Complex (Kind=dp) :: czero
!     ..
      Data czero/(0.0E0_dp, 0.0E0_dp)/
!     ..
      m = mod(n, 5)
      If (m/=0) Then
        Do i = 1, m
          a(i) = czero
        End Do
        If (n<5) Return
      End If
      mp1 = m + 1
      Do i = mp1, n, 5
        a(i) = czero
        a(i+1) = czero
        a(i+2) = czero
        a(i+3) = czero
        a(i+4) = czero
      End Do

    End Subroutine
