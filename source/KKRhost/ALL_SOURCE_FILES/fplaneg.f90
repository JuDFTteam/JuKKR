subroutine fplaneg(lamda, g, pref, lmax, ga, vol)
  use :: mod_datatypes, only: dp
  ! **************************************************************************
  ! This sub calculates the derivatives of the inverce
  ! space contribution to the ewald sum


  ! l     (                                                      )
  ! d   pi*( exp(gz)*erfc(g/lam/lam + 2*z)*lam/2                  )  /
  ! (             +   exp(-gz)*erfc(g/lam/lam - 2*z)*lam/2 ) / (g Vol)
  ! ---    (                                                      )
  ! l
  ! dz

  ! And the limit z -> 0 is taken (lam is the lamda parameter

  ! *********************************************************************
  implicit none
  integer :: lmax
  real (kind=dp) :: alpha, g(0:4), pref(0:lmax)
  integer :: l
  real (kind=dp) :: derfc, lamda, er, ex, pi, ga, vol
  real (kind=dp) :: sqpi

  do l = 0, 4
    g(l) = 0.e0_dp
  end do
  do l = 0, lmax
    pref(l) = 0.e0_dp
  end do
  pi = 4.e0_dp*atan(1.e0_dp)
  sqpi = sqrt(pi)
  alpha = ga/2.e0_dp/lamda
  er = derfc(alpha)
  ex = exp(-alpha*alpha)

  if (abs(ga)>1.e-6_dp) then
    g(0) = 2.e0_dp*pi/vol*er/ga
  else
    g(0) = 0.e0_dp
  end if

  pref(2) = sqrt(5.e0_dp/pi)/2.e0_dp/2.e0_dp
  g(2) = pref(2)/vol*(2.e0_dp*pi*ga*er-ex*4.e0_dp*sqpi*lamda)

  pref(4) = 3.e0_dp*sqrt(9.e0_dp/pi)/16.e0_dp/9.e0_dp
  g(4) = pref(4)/vol*(2.e0_dp*pi*ga*ga*ga*er-ex*4.e0_dp*sqpi*(lamda*ga*ga- &
    2.e0_dp*lamda**3))

  ! pref(6) = sqrt(13.d0/pi)/32.d0/13.d0
  ! g(6) = pref(6)/vol*(2.d0*pi*ga**5*er - ex*4.d0*sqpi*
  ! &                               (lamda*ga**4
  ! &                               -2.d0*ga*ga*lamda*lamda
  ! &                               +12.d0*lamda**5))
end subroutine fplaneg
