DOUBLE PRECISION FUNCTION erfcex(z)
!-----------------------------------------------------------------------

!     calculates complementary errorfunction times sqrt(pi)
!      times exp(z*z)  by continued fractions

!-----------------------------------------------------------------------
implicit none

!.. scalar arguments ..
double precision z

!.. local scalars ..
double precision bound,erf1,exzz,f,fa,q,ratio,sqrtpi,term,u,ua, &
                 v,x,xa,y,z2,zz

!.. intrinsic functions ..
intrinsic abs,atan,exp,sqrt


bound=3.d-11
sqrtpi = SQRT(4.0D0*ATAN(1.0D0))
zz = z*z

!---> choose algorithm

IF (z < 1.5D0) THEN
  
!     this exponential was outside the if statement
!     but for large arguments the exponent blow up
!     changes made 21/10/99
  
  exzz = EXP(zz)
  
  z2 = 2.0D0*zz
  erf1 = z
  ratio = 1.0D0
  term = z
  10    CONTINUE
  ratio = ratio + 2.0D0
  term = term*z2/ratio
  erf1 = erf1 + term
  IF (term > bound) GO TO 10
  erfcex = sqrtpi*exzz - 2.0D0*erf1
  
ELSE
  
!---> continued fraction expansion : abramowitz p. 298, eq. (7.1.14)
  
  u = 1.0D0
  v = 0.0D0
  x = z
  y = 1.0D0
  q = 0.5D0
  f = (u+v*q)/ (x+y*q)
  20    CONTINUE
  ua = u
  u = u*z + v*q
  v = ua
  xa = x
  x = x*z + y*q
  y = xa
  q = q + 0.5D0
  fa = f
  f = (u+v*q)/ (x+y*q)
  IF (ABS(fa-f) > bound*f) GO TO 20
  erfcex = f
END IF

END FUNCTION erfcex
