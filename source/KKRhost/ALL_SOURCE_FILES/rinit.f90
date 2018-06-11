    Subroutine rinit(n, a)
      Use mod_datatypes, Only: dp
! **********************************************************************
! * Setting the first N values of a double precision array A to zero   *
! **********************************************************************
!..
!.. Arguments ..
      Integer :: n
      Real (Kind=dp) :: a(*)
!..
!.. Locals ..
      Integer :: i, m, mp1
      Real (Kind=dp) :: dzero
!..
      Data dzero/0.0E0_dp/
!..
!..
      m = mod(n, 5)
      If (m/=0) Then
        Do i = 1, m
          a(i) = dzero
        End Do
        If (n<5) Return
      End If
      mp1 = m + 1
      Do i = mp1, n, 5
        a(i) = dzero
        a(i+1) = dzero
        a(i+2) = dzero
        a(i+3) = dzero
        a(i+4) = dzero
      End Do

    End Subroutine
