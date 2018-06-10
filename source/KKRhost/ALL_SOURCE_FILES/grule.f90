! **********************************************************************
SUBROUTINE grule(n,x,w)

!***********************************************************************

!     determines the (n+1)/2 nonnegative points x(i) and
!     the corresponding weights w(i) of the n-point
!     gauss-legendre integration rule, normalized to the
!     interval [-1,1]. the x(i) appear in descending order.

!     this routine is from 'methods of numerical integration',
!     p.j. davis and p. rabinowitz, page 369.

!***********************************************************************

!.. Scalar Arguments ..
INTEGER N
!..
!.. Array Arguments ..
DOUBLE PRECISION W(*),X(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION D1,D2PN,D3PN,D4PN,DEN,DP,DPN,E1,FX,H,P,PI,PK, &
                 PKM1,PKP1,T,T1,U,V,X0
INTEGER I,IT,K,M
!..
!.. Intrinsic Functions ..
INTRINSIC COS,ATAN

pi = 4.d0*ATAN(1.d0)
m  = (n+1)/2
e1 = n* (n+1)
DO  i = 1,m
  t = (4*i-1)*pi/ (4*n+2)
  x0 = (1.0D0- (1.0D0-1.0D0/n)/ (8.0D0*n*n))*COS(t)
  
!--->    iterate on the value  (m.w. jan. 1982)
  
  DO  it = 1,2
    pkm1 = 1.
    pk = x0
    DO  k = 2,n
      t1 = x0*pk
      pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1
      pkm1 = pk
      pk = pkp1
    END DO
    den = 1. - x0*x0
    d1 = n* (pkm1-x0*pk)
    dpn = d1/den
    d2pn = (2.*x0*dpn-e1*pk)/den
    d3pn = (4.*x0*d2pn+ (2.-e1)*dpn)/den
    d4pn = (6.*x0*d3pn+ (6.-e1)*d2pn)/den
    u = pk/dpn
    v = d2pn/dpn
    h = -u* (1.+.5*u* (v+u* (v*v-u*d3pn/ (3.*dpn))))
    p = pk + h* (dpn+.5*h* (d2pn+h/3.* (d3pn+.25*h*d4pn)))
    dp = dpn + h* (d2pn+.5*h* (d3pn+h*d4pn/3.))
    h = h - p/dp
    x0 = x0 + h
  END DO
  x(i) = x0
  fx = d1 - h*e1* (pk+.5*h* (dpn+h/3.* (d2pn+.25*h* (d3pn+ .2*h*d4pn))))
  w(i) = 2.* (1.-x(i)*x(i))/ (fx*fx)
END DO
IF (m+m > n) x(m) = 0.
END SUBROUTINE grule
