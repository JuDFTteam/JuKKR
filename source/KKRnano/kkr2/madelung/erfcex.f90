  double precision function erfcex(z)
  implicit none
    double precision, intent(in) :: z
    !-----------------------------------------------------------------------
    !     calculates complementary errorfunction times sqrt(pi)
    !      times exp(z*z)  by continued fractions
    !-----------------------------------------------------------------------
    double precision, parameter :: bound = 3.d-11
    double precision :: erf1, exzz, f, fa, q, denom, sqrtpi, term, u, ua, v, x, xa, y, tz2, z2
    logical :: run
    
    sqrtpi = sqrt(4.d0*atan(1.d0))
    z2 = z*z
    
    !---> choose algorithm
    if (z < 1.5d0) then
      ! this exponential was outside the if statement but for large arguments the exponent blow up

      exzz = exp(z2)

      tz2 = 2.d0*z2
      erf1 = z
      denom = 1.d0
      term = z
      
      run = .true. ! always do at least one iteration
      do while (run)
        denom = denom + 2.d0
        term = term*tz2/denom
        erf1 = erf1 + term
        run = (term > bound)
      enddo ! while
      
      erfcex = sqrtpi*exzz - 2.d0*erf1

    else
      ! continued fraction expansion : abramowitz p. 298, eq. (7.1.14)
      u = 1.d0
      v = 0.d0
      x = z
      y = 1.d0
      q = 0.5d0
      f = (u + v*q)/(x + y*q)
      
      run = .true. ! always do at least one iteration
      do while (run)
        ua = u
        u = u*z + v*q
        v = ua
        xa = x
        x = x*z + y*q
        y = xa
        q = q + 0.5d0
        fa = f
        f = (u + v*q)/(x + y*q)
        run = (abs(fa - f) > bound*f)
      enddo ! while
      
      erfcex = f
    endif

  endfunction erfcex
