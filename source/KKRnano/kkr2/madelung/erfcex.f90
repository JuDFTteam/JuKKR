  double precision function erfcex(z)
  implicit none
    double precision, intent(in) :: z
    !-----------------------------------------------------------------------
    !     calculates complementary errorfunction times sqrt(pi)
    !      times exp(z*z)  by continued fractions
    !-----------------------------------------------------------------------
    double precision, parameter :: bound = 3.d-11
    double precision :: erf1, exzz, f, fa, q, ratio, sqrtpi, term, u, ua, v, x, xa, y, tz2, z2

    sqrtpi = sqrt(4.d0*atan(1.d0))
    z2 = z*z

    !---> choose algorithm
    if (z < 1.5d0) then
      ! this exponential was outside the if statement but for large arguments the exponent blow up

      exzz = exp(z2)

      tz2 = 2.d0*z2
      erf1 = z
      ratio = 1.d0
      term = z
      
   10 continue
      ratio = ratio + 2.d0
      term = term*tz2/ratio
      erf1 = erf1 + term
      if (term > bound) go to 10
      
      erfcex = sqrtpi*exzz - 2.d0*erf1

    else
      ! continued fraction expansion : abramowitz p. 298, eq. (7.1.14)
      u = 1.d0
      v = 0.d0
      x = z
      y = 1.d0
      q = 0.5d0
      f = (u + v*q)/(x + y*q)
      
   20 continue
      ua = u
      u = u*z + v*q
      v = ua
      xa = x
      x = x*z + y*q
      y = xa
      q = q + 0.5d0
      fa = f
      f = (u + v*q)/(x + y*q)
      if (abs(fa - f) > bound*f) go to 20
      
      erfcex = f
    endif

  endfunction erfcex
