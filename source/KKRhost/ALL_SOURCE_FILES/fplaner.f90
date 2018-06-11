    Subroutine fplaner(alpha, g, r)
      Use mod_datatypes, Only: dp
! ************************************************
! This sub calculates the derivatives of the real
! space contribution to the ewald sum .

!              l
!             d     erfc(lamda*sqrt(d*d+z*z))
!      lim    --   ------------------------
!      z->0     l        sqrt(d*d+z*z)
!             dz

!  Up to l = 4 (l=1,3,5,7 etc vanish)



! ************************************************

      Implicit None
      Real (Kind=dp) :: alpha, g(0:4), r
      Integer :: l
      Real (Kind=dp) :: derfc, lamda, er, ex, pi, pref, sqpi

      Do l = 0, 4
        g(l) = 0.E0_dp
      End Do
      pi = 4.E0_dp*atan(1.E0_dp)
      sqpi = sqrt(pi)
      er = derfc(alpha)
      ex = exp(-alpha*alpha)
      lamda = alpha/r

      g(0) = er/r

      pref = sqrt(5.E0_dp/pi)/4.E0_dp
      g(2) = -pref*(er/r/r/r+ex*2.E0_dp*lamda/r/r/sqpi)

      pref = 3.E0_dp*sqrt(9.E0_dp/pi)/16.E0_dp/9.E0_dp
      g(4) = pref*(9.E0_dp*er+ex*(12.E0_dp*alpha**3+18.E0_dp*alpha)/sqpi)/r/r/ &
        r/r/r

    End Subroutine
