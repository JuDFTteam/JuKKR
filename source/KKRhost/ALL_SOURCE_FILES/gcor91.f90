SUBROUTINE gcor91(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)
!-----------------------------------------------------------------
!called by corlsd
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!.. Scalar Arguments ..
      DOUBLE PRECISION A,A1,B1,B2,B3,B4,GG,GGRS,P,RS
!..
!.. Local Scalars ..
      DOUBLE PRECISION P1,Q0,Q1,Q2,Q3,RS12,RS32,RSP
!..
!.. Intrinsic Functions ..
      INTRINSIC LOG,SQRT
!..
p1 = p + 1.d0
q0 = -2.d0*a* (1.d0+a1*rs)
rs12 = SQRT(rs)
rs32 = rs12**3
rsp = rs**p
q1 = 2.d0*a* (b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
q2 = LOG(1.d0+1.d0/q1)
gg = q0*q2
q3 = a* (b1/rs12+2.d0*b2+3.d0*b3*rs12+2.d0*b4*p1*rsp)
ggrs = -2.d0*a*a1*q2 - q0*q3/ (q1**2+q1)
RETURN
END SUBROUTINE gcor91
