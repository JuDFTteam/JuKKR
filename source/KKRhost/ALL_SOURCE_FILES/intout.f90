SUBROUTINE intout(g,f,v,e,l,nne,k2,dg,a,b,z,nsra)
implicit none
!.. Scalar Arguments ..
      DOUBLE PRECISION A,B,DG,E,Z
      INTEGER K2,L,NNE,NSRA
!..
!.. Array Arguments ..
      DOUBLE PRECISION F(*),G(*),V(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION AA,ALFA,B1,B2,BB,BETA,CVLIGHT,DET,DF1,DF2,DF3, &
                       DG1,DG2,DG3,DR,EA,FLLP1,H83,P12,P21,PHI,PP,QQ,R, &
                       R1,R2,R3,R83SQ,RPB,S,SG,SGM1,U,X,Y,ZZ
      INTEGER I,I1,K,KM1,N
!..
!.. Local Arrays ..
      DOUBLE PRECISION D(2,3),PX(20),QX(20)
!..
!.. Intrinsic Functions ..
      INTRINSIC DSIGN,EXP,SQRT
!..
!     ..
zz = z + z
cvlight = 274.0720442D0
IF (nsra == 1) cvlight = 1.0D0
ea = EXP(a)
fllp1 = l* (l+1.d0)
r83sq = 64.d0/9.d0
r1 = 1.d0/9.d0
r2 = -5.d0*r1
r3 = 19.d0*r1
h83 = 8.d0/3.d0
aa = -zz/cvlight
bb = fllp1 - aa*aa
p21 = (v(1)-e)/cvlight
p12 = cvlight - p21
px(1) = 0.d0
qx(1) = 0.d0
IF (z <= 20.d0 .OR. nsra == 1) THEN
  s = 1.d0*l
  px(2) = 0.d0
  px(3) = 1.d0
  DO  k = 2,9
    px(k+2) = ((v(1)-e)*px(k)-zz*px(k+1))/ (k+l+l)/ (k-1.d0)
  END DO
  DO  k = 2,10
    qx(k) = px(k+1)* (l+k-2.d0)/cvlight
  END DO
  
ELSE
  s = SQRT(fllp1+1.d0-aa*aa)
  px(2) = 1.d0
  qx(2) = (1.d0-s)/aa
  DO  i = 3,10
    n = i - 2
    alfa = p12*qx(i-1)
    beta = p21*aa*px(i-1)
    IF (l /= 0) beta = beta - p12*aa*px(i-1) - p12*p21*px(i-2) +  &
        (n+s)*p12*qx(i-1)
    det = n* (n+s+s)*aa
    px(i) = (alfa* (n+s+1.d0)*aa-aa*beta)/det
    qx(i) = (beta* (n+s-1.d0)-bb*alfa)/det
  END DO
END IF

g(1) = 0.d0
f(1) = 0.d0
rpb = b
DO  k = 2,4
  rpb = rpb*ea
  r = rpb - b
  dr = a*rpb
  phi = (e+zz/r-v(k))*dr/cvlight
  u = dr*cvlight + phi
  IF (nsra == 1) u = dr
  x = -dr/r
  y = -fllp1*x*x/u + phi
  pp = px(10)
  qq = qx(10)
  DO  i1 = 3,10
    i = 12 - i1
    pp = pp*r + px(i)
    qq = qq*r + qx(i)
  END DO
  g(k) = (r**s)*pp
  f(k) = (r**s)*qq
  sg = DSIGN(1.0D0,g(k))
  sgm1 = DSIGN(1.0D0,g(k-1))
  IF (sg*sgm1 < 0.d0) nne = nne + 1
  d(1,k-1) = u*f(k) - x*g(k)
  d(2,k-1) = x*f(k) - y*g(k)
END DO
dg1 = d(1,1)
dg2 = d(1,2)
dg3 = d(1,3)
df1 = d(2,1)
df2 = d(2,2)
df3 = d(2,3)
DO  k = 5,k2
  km1 = k - 1
  rpb = rpb*ea
  r = rpb - b
  dr = a*rpb
  phi = (e+zz/r-v(k))*dr/cvlight
  u = dr*cvlight + phi
  IF (nsra == 1) u = dr
  x = -dr/r
  y = -fllp1*x*x/u + phi
  det = r83sq - x*x + u*y
  b1 = g(km1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
  b2 = f(km1)*h83 + r1*df1 + r2*df2 + r3*df3
  g(k) = (b1* (h83-x)+b2*u)/det
  f(k) = (b2* (h83+x)-b1*y)/det
  sg = DSIGN(1.0D0,g(k))
  sgm1 = DSIGN(1.0D0,g(km1))
  IF (sg*sgm1 < 0.d0) nne = nne + 1
  dg1 = dg2
  dg2 = dg3
  dg3 = u*f(k) - x*g(k)
  df1 = df2
  df2 = df3
  df3 = x*f(k) - y*g(k)
END DO
dg = dg3
END SUBROUTINE intout
