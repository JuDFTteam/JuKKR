    Subroutine fplaneg(lamda, g, pref, lmax, ga, vol)
      Use mod_datatypes, Only: dp
! **************************************************************************
! This sub calculates the derivatives of the inverce
! space contribution to the ewald sum


!     l     (                                                      )
!    d   pi*( exp(gz)*erfc(g/lam/lam + 2*z)*lam/2                  )  /
!           (             +   exp(-gz)*erfc(g/lam/lam - 2*z)*lam/2 ) / (g Vol)
!    ---    (                                                      )
!      l
!    dz

!   And the limit z -> 0 is taken (lam is the lamda parameter

! *********************************************************************
      Implicit None
      Integer :: lmax
      Real (Kind=dp) :: alpha, g(0:4), pref(0:lmax)
      Integer :: l
      Real (Kind=dp) :: derfc, lamda, er, ex, pi, ga, vol
      Real (Kind=dp) :: sqpi

      Do l = 0, 4
        g(l) = 0.E0_dp
      End Do
      Do l = 0, lmax
        pref(l) = 0.E0_dp
      End Do
      pi = 4.E0_dp*atan(1.E0_dp)
      sqpi = sqrt(pi)
      alpha = ga/2.E0_dp/lamda
      er = derfc(alpha)
      ex = exp(-alpha*alpha)

      If (abs(ga)>1.E-6_dp) Then
        g(0) = 2.E0_dp*pi/vol*er/ga
      Else
        g(0) = 0.E0_dp
      End If

      pref(2) = sqrt(5.E0_dp/pi)/2.E0_dp/2.E0_dp
      g(2) = pref(2)/vol*(2.E0_dp*pi*ga*er-ex*4.E0_dp*sqpi*lamda)

      pref(4) = 3.E0_dp*sqrt(9.E0_dp/pi)/16.E0_dp/9.E0_dp
      g(4) = pref(4)/vol*(2.E0_dp*pi*ga*ga*ga*er-ex*4.E0_dp*sqpi*(lamda*ga*ga- &
        2.E0_dp*lamda**3))

!      pref(6) = sqrt(13.d0/pi)/32.d0/13.d0
!      g(6) = pref(6)/vol*(2.d0*pi*ga**5*er - ex*4.d0*sqpi*
!     &                               (lamda*ga**4
!     &                               -2.d0*ga*ga*lamda*lamda
!     &                               +12.d0*lamda**5))
    End Subroutine
