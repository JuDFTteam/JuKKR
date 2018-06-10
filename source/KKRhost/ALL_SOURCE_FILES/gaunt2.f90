!***********************************************************************
SUBROUTINE gaunt2(w,yr,n)
! ************************************************************************
!     sets up values needed for gaunt
!        m. weinert  january 1982

!     changed for calculating with real spherical harmonics
!                                           b.drittler  july 1987

!     W(N)        integration weights on 4*LMAXD points in the intervall
!                 (-1,0) (from routine GRULE)

!     YR(N,L,M)   spherical harmonics on 4*LMAXD points to angular
!                 momentum indices (l,m) scaled with a factor
!                 of RF=(4*pi)**(1/3)

!-----------------------------------------------------------------------
      IMPLICIT NONE
!.. Arguments
      INTEGER N
      DOUBLE PRECISION W(*),YR(N,0:N,0:N)
!..
!.. Local Scalars ..
      DOUBLE PRECISION A,CD,CTH,FAC,FPI,RF,STH,T
      INTEGER K,L,LOMAX,M
!..
!.. Local Arrays ..
      DOUBLE PRECISION P(0:N+1,0:N),X(N)
!..
!.. External Subroutines ..
      EXTERNAL GRULE
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!     ..
fpi = 16D0*ATAN(1D0)
rf = fpi**(1D0/3D0)
lomax = n

!--->    obtain gauss-legendre points and weights

CALL grule(2*n,x,w)

!--->    generate associated legendre functions for m.ge.0

DO k = 1,n
  cth = x(k)
  sth = SQRT(1.d0-cth*cth)
  fac = 1.d0
  
!--->    loop over m values
  
  DO m = 0,lomax
    fac = - DBLE(2*m-1)*fac
    p(m,m) = fac
    p(m+1,m) = DBLE(2*m+1)*cth*fac
    
!--->    recurse upward in l
    
    DO l = m + 2,lomax
      p(l,m) = ( DBLE(2*l-1)*cth*p(l-1,m)  &
          - DBLE(l+m-1)    *p(l-2,m) ) / DBLE(l-m)
    END DO
    
    fac = fac*sth
  END DO
  
!--->    multiply in the normalization factors
  
  DO l = 0,lomax
    a = rf*SQRT((2*l+1)/fpi)
    cd = 1.d0
    yr(k,l,0) = a*p(l,0)
    
    DO m = 1,l
      t = DBLE( (l+1-m)* (l+m))
      cd = cd/t
      yr(k,l,m) = a*SQRT(2.d0*cd)*p(l,m)
    END DO
  END DO
END DO
END SUBROUTINE gaunt2
