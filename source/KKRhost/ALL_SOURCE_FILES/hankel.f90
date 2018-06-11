    Subroutine hankel(h, l, arg)
      Use mod_datatypes, Only: dp
!  this subroutine uses the explicit formulas for the hankel
!  functions. for higher l-values these formulas may lead to
!  loss of significant figures. This subroutine should be used
!  only for core states.
      Implicit None
!.. Scalar Arguments ..
      Complex (Kind=dp) :: arg
      Integer :: l
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: h(*)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: a1, a2, a3, a4
!..
!.. Intrinsic Functions ..
      Intrinsic :: exp
!..
!.. Parameters ..
      Complex (Kind=dp) :: ci
      Parameter (ci=(0.0E0_dp,1.0E0_dp))
!     ..
      h(1) = -exp(arg*ci)/arg
      If (l/=1) Then
        a1 = (1.E0_dp, 0.E0_dp) - arg*ci
        h(2) = h(1)*a1/arg
        If (l/=2) Then
          a1 = 3.E0_dp*a1
          a2 = arg*arg
          h(3) = h(1)*(a1-a2)/a2
          If (l/=3) Then
            a1 = 5.E0_dp*a1
            a3 = a2*arg*ci
            a4 = a2*arg
            a2 = 6.E0_dp*a2
            h(4) = h(1)*(a1-a2+a3)/a4
            If (l/=4) Then
              a1 = 7.E0_dp*a1
              a2 = 7.5E0_dp*a2
              a3 = 10.E0_dp*a3
              a4 = a4*arg
              h(5) = h(1)*(a1-a2+a3+a4)/a4
              If (l/=5) Then
                h(6) = (9.0E0_dp, 0.0E0_dp)*h(5)/arg - h(4)
                If (l/=6) Then
                  Write (6, Fmt=100) l
                  Stop 'HANKEL'

                End If

              End If

            End If

          End If

        End If

      End If

      Return


100   Format (2X, ' hankel :  l=', I2, ' is too large')
    End Subroutine
