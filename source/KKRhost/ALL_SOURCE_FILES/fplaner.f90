subroutine fplaner(alpha, g, r)
  use :: mod_datatypes, only: dp
  ! ************************************************
  ! This sub calculates the derivatives of the real
  ! space contribution to the ewald sum .

  ! l
  ! d     erfc(lamda*sqrt(d*d+z*z))
  ! lim    --   ------------------------
  ! z->0     l        sqrt(d*d+z*z)
  ! dz

  ! Up to l = 4 (l=1,3,5,7 etc vanish)



  ! ************************************************

  implicit none
  real (kind=dp) :: alpha, g(0:4), r
  integer :: l
  real (kind=dp) :: derfc, lamda, er, ex, pi, pref, sqpi

  do l = 0, 4
    g(l) = 0.e0_dp
  end do
  pi = 4.e0_dp*atan(1.e0_dp)
  sqpi = sqrt(pi)
  er = derfc(alpha)
  ex = exp(-alpha*alpha)
  lamda = alpha/r

  g(0) = er/r

  pref = sqrt(5.e0_dp/pi)/4.e0_dp
  g(2) = -pref*(er/r/r/r+ex*2.e0_dp*lamda/r/r/sqpi)

  pref = 3.e0_dp*sqrt(9.e0_dp/pi)/16.e0_dp/9.e0_dp
  g(4) = pref*(9.e0_dp*er+ex*(12.e0_dp*alpha**3+18.e0_dp*alpha)/sqpi)/r/r/r/r/ &
    r

end subroutine fplaner
