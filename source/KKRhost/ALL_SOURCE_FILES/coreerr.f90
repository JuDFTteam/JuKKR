module mod_coreerr

contains

subroutine coreerr(err, var, s, nsol, pow, qow, piw, qiw)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! *   CALCULATE THE MISMATCH OF THE RADIAL WAVE FUNCTIONS AT THE     *
  ! *   POINT  NMATCH  FOR OUT- AND INWARD INTEGRATION                 *
  ! *                                                                  *
  ! ********************************************************************
  implicit none
  ! Dummy arguments
  integer :: nsol, s
  real (kind=dp) :: err(4), piw(2, 2), pow(2, 2), qiw(2, 2), qow(2, 2), var(4)

  ! Local variables
  integer :: t

  err(1) = pow(s, s) - piw(s, s)*var(2)
  err(2) = qow(s, s) - qiw(s, s)*var(2)

  if (nsol==1) return

  t = 3 - s

  err(1) = err(1) + pow(s, t)*var(3) - piw(s, t)*var(2)*var(4)
  err(2) = err(2) + qow(s, t)*var(3) - qiw(s, t)*var(2)*var(4)
  err(3) = pow(t, s) - piw(t, s)*var(2) + pow(t, t)*var(3) - &
    piw(t, t)*var(2)*var(4)
  err(4) = qow(t, s) - qiw(t, s)*var(2) + qow(t, t)*var(3) - &
    qiw(t, t)*var(2)*var(4)

end subroutine coreerr

end module mod_coreerr
