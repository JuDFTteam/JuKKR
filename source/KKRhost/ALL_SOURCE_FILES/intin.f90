SUBROUTINE intin(g,f,v,e,l,nne,valu,slop,k1,k2,kc,dg,a,b,z,nsra)
!.. Scalar Arguments ..
      DOUBLE PRECISION A,B,DG,E,SLOP,VALU,Z
      INTEGER K1,K2,KC,L,NNE,NSRA
!..
!.. Array Arguments ..
      DOUBLE PRECISION F(*),G(*),V(*)
!..
!.. Local Scalars ..
     DOUBLE PRECISION AF1,AF2,AF3,AG1,AG2,AG3,B1,B2,CVLIGHT,DET,DF1, &
                      DF2,DF3,DG1,DG2,DG3,DR,EA,FF,FLLP1,GG,H83,PHI,Q, &
                      R,R1,R2,R3,R83SQ,RPB,SDG3,SG,SGP1,U,VB,X,Y,ZZ
      INTEGER I,K,KP1
!..
!.. Local Arrays ..
      DOUBLE PRECISION D(2,3)
!..
!.. Intrinsic Functions ..
      INTRINSIC DSIGN,EXP,SQRT
!     ..
zz = z + z
cvlight = 274.0720442D0
IF (nsra == 1) cvlight = 1.0D0
fllp1 = l* (l+1.d0)
r83sq = 64.d0/9.d0
r1 = 1.d0/9.d0
r2 = -5.d0*r1
r3 = 19.d0*r1
h83 = -8.d0/3.d0
ea = EXP(a)
rpb = b*EXP(a*k1-a)
r = rpb - b
dr = a*rpb
phi = (e+zz/r-v(k1))*dr/cvlight
u = dr*cvlight + phi
IF (nsra == 1) u = dr
x = -dr/r
y = -fllp1*x*x/u + phi
g(k1) = valu
f(k1) = (slop*dr+x*valu)/u
q = 1.d0/SQRT(ea)
ag1 = slop*dr
af1 = x*f(k1) - y*g(k1)
k = k1
dg3 = ag1
IF (k2 /= k1) THEN
  DO  i = 1,3
    kp1 = k
    k = k - 1
    rpb = rpb*q
    dr = rpb*a
    r = rpb - b
    gg = g(kp1) - .5D0*ag1
    ff = f(kp1) - .5D0*af1
    vb = (3.d0*v(kp1)+6.d0*v(k)-v(k-1))*.125D0
    phi = (e+zz/r-vb)*dr/cvlight
    u = dr*cvlight + phi
    IF (nsra == 1) u = dr
    x = -dr/r
    y = -fllp1*x*x/u + phi
    ag2 = u*ff - x*gg
    af2 = x*ff - y*gg
    gg = g(kp1) - .5D0*ag2
    ff = f(kp1) - .5D0*af2
    ag3 = u*ff - x*gg
    af3 = x*ff - y*gg
    rpb = rpb*q
    dr = a*rpb
    r = rpb - b
    phi = (e+zz/r-v(k))*dr/cvlight
    u = dr*cvlight + phi
    IF (nsra == 1) u = dr
    x = -dr/r
    y = -fllp1*x*x/u + phi
    gg = g(kp1) - ag3
    ff = f(kp1) - af3
    g(k) = g(kp1) - (ag1+2.d0* (ag2+ag3)+u*ff-x*gg)/6.d0
    f(k) = f(kp1) - (af1+2.d0* (af2+af3)+x*ff-y*gg)/6.d0
    sg = DSIGN(1.0D0,g(k))
    sgp1 = DSIGN(1.0D0,g(kp1))
    IF (sg*sgp1 < 0.d0) nne = nne + 1
    ag1 = u*f(k) - x*g(k)
    af1 = x*f(k) - y*g(k)
    IF (k == k2) THEN
      GO TO 30
      
    ELSE
      d(1,i) = ag1
      d(2,i) = af1
    END IF
    
  END DO
  q = 1.d0/ea
  dg1 = d(1,1)
  dg2 = d(1,2)
  dg3 = d(1,3)
  df1 = d(2,1)
  df2 = d(2,2)
  df3 = d(2,3)
  20   CONTINUE
  kp1 = k
  k = k - 1
  rpb = rpb*q
  dr = a*rpb
  r = rpb - b
  phi = (e+zz/r-v(k))*dr/cvlight
  u = dr*cvlight + phi
  IF (nsra == 1) u = dr
  x = -dr/r
  y = -fllp1*x*x/u + phi
  det = r83sq - x*x + u*y
  b1 = g(kp1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
  b2 = f(kp1)*h83 + r1*df1 + r2*df2 + r3*df3
  g(k) = (b1* (h83-x)+b2*u)/det
  f(k) = (b2* (h83+x)-b1*y)/det
  sg = DSIGN(1.0D0,g(k))
  sgp1 = DSIGN(1.0D0,g(kp1))
  IF (sg*sgp1 < 0.d0) nne = nne + 1
  dg1 = dg2
  df1 = df2
  dg2 = dg3
  df2 = df3
  dg3 = u*f(k) - x*g(k)
  df3 = x*f(k) - y*g(k)
  sdg3 = DSIGN(1.0D0,dg3)
  IF (k > k2 .AND. sg*sdg3 < 0.d0) GO TO 20
END IF
30 kc = k
dg = dg3
END SUBROUTINE intin
