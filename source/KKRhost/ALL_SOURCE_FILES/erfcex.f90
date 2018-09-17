module mod_erfcex
  use :: mod_datatypes, only: dp
  private :: dp

contains

  function erfcex(z)
    ! -----------------------------------------------------------------------

    ! calculates complementary errorfunction times sqrt(pi)
    ! times exp(z*z)  by continued fractions

    ! -----------------------------------------------------------------------
    implicit none
    real (kind=dp) :: erfcex

    ! .. scalar arguments ..
    real (kind=dp) :: z

    ! .. local scalars ..
    real (kind=dp) :: bound, erf1, exzz, f, fa, q, ratio, sqrtpi, term, u, ua, v, x, xa, y, z2, zz

    ! .. intrinsic functions ..
    intrinsic :: abs, atan, exp, sqrt


    bound = 3.e-11_dp
    sqrtpi = sqrt(4.0e0_dp*atan(1.0e0_dp))
    zz = z*z

    ! ---> choose algorithm

    if (z<1.5e0_dp) then

      ! this exponential was outside the if statement
      ! but for large arguments the exponent blow up
      ! changes made 21/10/99

      exzz = exp(zz)

      z2 = 2.0e0_dp*zz
      erf1 = z
      ratio = 1.0e0_dp
      term = z
100   continue
      ratio = ratio + 2.0e0_dp
      term = term*z2/ratio
      erf1 = erf1 + term
      if (term>bound) go to 100
      erfcex = sqrtpi*exzz - 2.0e0_dp*erf1

    else

      ! ---> continued fraction expansion : abramowitz p. 298, eq. (7.1.14)

      u = 1.0e0_dp
      v = 0.0e0_dp
      x = z
      y = 1.0e0_dp
      q = 0.5e0_dp
      f = (u+v*q)/(x+y*q)
110   continue
      ua = u
      u = u*z + v*q
      v = ua
      xa = x
      x = x*z + y*q
      y = xa
      q = q + 0.5e0_dp
      fa = f
      f = (u+v*q)/(x+y*q)
      if (abs(fa-f)>bound*f) go to 110
      erfcex = f
    end if

  end function erfcex

end module mod_erfcex
