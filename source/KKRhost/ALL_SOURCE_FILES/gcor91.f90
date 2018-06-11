    Subroutine gcor91(a, a1, b1, b2, b3, b4, p, rs, gg, ggrs)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------
!called by corlsd
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!.. Scalar Arguments ..
      Real (Kind=dp) :: a, a1, b1, b2, b3, b4, gg, ggrs, p, rs
!..
!.. Local Scalars ..
      Real (Kind=dp) :: p1, q0, q1, q2, q3, rs12, rs32, rsp
!..
!.. Intrinsic Functions ..
      Intrinsic :: log, sqrt
!..
      p1 = p + 1.E0_dp
      q0 = -2.E0_dp*a*(1.E0_dp+a1*rs)
      rs12 = sqrt(rs)
      rs32 = rs12**3
      rsp = rs**p
      q1 = 2.E0_dp*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
      q2 = log(1.E0_dp+1.E0_dp/q1)
      gg = q0*q2
      q3 = a*(b1/rs12+2.E0_dp*b2+3.E0_dp*b3*rs12+2.E0_dp*b4*p1*rsp)
      ggrs = -2.E0_dp*a*a1*q2 - q0*q3/(q1**2+q1)
      Return
    End Subroutine
