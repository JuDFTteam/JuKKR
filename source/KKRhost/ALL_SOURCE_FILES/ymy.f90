! **********************************************************************
SUBROUTINE ymy(v1,v2,v3,r,ylm,lmax)
! **********************************************************************
!    this subroutine calculates real spherical harmonics with the
!     normalization : <y|y> =1
!    returns also r = length of vector v

!     generate the complex spherical harmonics for the vector v
!     using a stable upward recursion in l.  (see notes
!     by m. weinert.)
!                                  m.weinert  1982

!     converted to real spherical harmonics .
!                                  b.drittler 1987
!-----------------------------------------------------------------------

      IMPLICIT NONE
!.. Parameters ..
      DOUBLE PRECISION SZERO
      PARAMETER (SZERO=1.0D-20)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION, intent(out) :: R
      DOUBLE PRECISION, intent(in) :: V1,V2,V3
      INTEGER, intent(in) :: LMAX
!..
!.. Array Arguments ..
      DOUBLE PRECISION, intent(out) :: YLM((2*LMAX+1)**2)
!..
!.. Local Scalars ..
      DOUBLE PRECISION A,CD,CPH,CTH,FAC,FPI,PI,RTWO,SGM,SPH,STH,T,XY,XYZ
      INTEGER I,L,M
!..
!.. Local Arrays ..
      DOUBLE PRECISION C(0:LMAX),P(0:LMAX,0:LMAX),S(0:LMAX)
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!..
!.. External Subroutines ..
      EXTERNAL RCSTOP
!     ..
pi = 4.d0*ATAN(1.d0)
fpi = 4.d0*pi
rtwo = SQRT(2.0D0)

!--->    calculate sin and cos of theta and phi

xy = v1**2 + v2**2
xyz = xy + v3**2

r = SQRT(xyz)
IF (xyz <= 0.0D0) THEN
  WRITE(*,*) xyz
  CALL rcstop('ylm=0   ')
  
ELSE
  
  IF (xy > szero*xyz) THEN
    xy = SQRT(xy)
    xyz = SQRT(xyz)
    cth = v3/xyz
    sth = xy/xyz
    cph = v1/xy
    sph = v2/xy
    
  ELSE
    
    sth = 0.0D0
    cth = 1.0D0
    IF (v3 < 0) cth = -1.0D0
    cph = 1.0D0
    sph = 0.0D0
  END IF
  
!--->    generate associated legendre functions for m.ge.0
!        loop over m values
  
  fac = 1.0D0
  DO  m = 0,lmax - 1
    fac = - (2*m-1)*fac
    p(m,m) = fac
    p(m+1,m) = (2*m+1)*cth*fac
    
!--->    recurse upward in l
    
    DO  l = m + 2,lmax
      p(l,m) = ((2*l-1)*cth*p(l-1,m)- (l+m-1)*p(l-2,m))/ (l-m)
    END DO
    fac = fac*sth
  END DO
  p(lmax,lmax) = - (2*lmax-1)*fac
  
!--->    determine sin and cos of phi
  
  s(0) = 0.0D0
  s(1) = sph
  c(0) = 1.0D0
  c(1) = cph
  DO  m = 2,lmax
    s(m) = 2*cph*s(m-1) - s(m-2)
    c(m) = 2*cph*c(m-1) - c(m-2)
  END DO
  
!--->    multiply in the normalization factors
  
  i = 0
  DO  l = 0,lmax
    i = i + l + 1
    a = SQRT((2*l+1)/fpi)
    cd = 1
    ylm(i) = a*p(l,0)
    sgm = -rtwo
    DO  m = 1,l
      t = (l+1-m)* (l+m)
      cd = cd/t
      t = a*SQRT(cd)
      ylm(i+m) = sgm*t*p(l,m)*c(m)
      ylm(i-m) = sgm*t*p(l,m)*s(m)
      sgm = -sgm
    END DO
    i = i + l
  END DO
  
END IF

RETURN

END SUBROUTINE ymy
