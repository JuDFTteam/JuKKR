subroutine intin(g, f, v, e, l, nne, valu, slop, k1, k2, kc, dg, a, b, z, &
  nsra)
!.. Scalar Arguments ..
  double precision :: a, b, dg, e, slop, valu, z
  integer :: k1, k2, kc, l, nne, nsra
!..
!.. Array Arguments ..
  double precision :: f(*), g(*), v(*)
!..
!.. Local Scalars ..
  double precision :: af1, af2, af3, ag1, ag2, ag3, b1, b2, cvlight, det, df1, &
    df2, df3, dg1, dg2, dg3, dr, ea, ff, fllp1, gg, h83, phi, q, r, r1, r2, &
    r3, r83sq, rpb, sdg3, sg, sgp1, u, vb, x, y, zz
  integer :: i, k, kp1
!..
!.. Local Arrays ..
  double precision :: d(2, 3)
!..
!.. Intrinsic Functions ..
  intrinsic :: dsign, exp, sqrt
!     ..
  zz = z + z
  cvlight = 274.0720442d0
  if (nsra==1) cvlight = 1.0d0
  fllp1 = l*(l+1.d0)
  r83sq = 64.d0/9.d0
  r1 = 1.d0/9.d0
  r2 = -5.d0*r1
  r3 = 19.d0*r1
  h83 = -8.d0/3.d0
  ea = exp(a)
  rpb = b*exp(a*k1-a)
  r = rpb - b
  dr = a*rpb
  phi = (e+zz/r-v(k1))*dr/cvlight
  u = dr*cvlight + phi
  if (nsra==1) u = dr
  x = -dr/r
  y = -fllp1*x*x/u + phi
  g(k1) = valu
  f(k1) = (slop*dr+x*valu)/u
  q = 1.d0/sqrt(ea)
  ag1 = slop*dr
  af1 = x*f(k1) - y*g(k1)
  k = k1
  dg3 = ag1
  if (k2/=k1) then
    do i = 1, 3
      kp1 = k
      k = k - 1
      rpb = rpb*q
      dr = rpb*a
      r = rpb - b
      gg = g(kp1) - .5d0*ag1
      ff = f(kp1) - .5d0*af1
      vb = (3.d0*v(kp1)+6.d0*v(k)-v(k-1))*.125d0
      phi = (e+zz/r-vb)*dr/cvlight
      u = dr*cvlight + phi
      if (nsra==1) u = dr
      x = -dr/r
      y = -fllp1*x*x/u + phi
      ag2 = u*ff - x*gg
      af2 = x*ff - y*gg
      gg = g(kp1) - .5d0*ag2
      ff = f(kp1) - .5d0*af2
      ag3 = u*ff - x*gg
      af3 = x*ff - y*gg
      rpb = rpb*q
      dr = a*rpb
      r = rpb - b
      phi = (e+zz/r-v(k))*dr/cvlight
      u = dr*cvlight + phi
      if (nsra==1) u = dr
      x = -dr/r
      y = -fllp1*x*x/u + phi
      gg = g(kp1) - ag3
      ff = f(kp1) - af3
      g(k) = g(kp1) - (ag1+2.d0*(ag2+ag3)+u*ff-x*gg)/6.d0
      f(k) = f(kp1) - (af1+2.d0*(af2+af3)+x*ff-y*gg)/6.d0
      sg = dsign(1.0d0, g(k))
      sgp1 = dsign(1.0d0, g(kp1))
      if (sg*sgp1<0.d0) nne = nne + 1
      ag1 = u*f(k) - x*g(k)
      af1 = x*f(k) - y*g(k)
      if (k==k2) then
        go to 110

      else
        d(1, i) = ag1
        d(2, i) = af1
      end if

    end do
    q = 1.d0/ea
    dg1 = d(1, 1)
    dg2 = d(1, 2)
    dg3 = d(1, 3)
    df1 = d(2, 1)
    df2 = d(2, 2)
    df3 = d(2, 3)
100 continue
    kp1 = k
    k = k - 1
    rpb = rpb*q
    dr = a*rpb
    r = rpb - b
    phi = (e+zz/r-v(k))*dr/cvlight
    u = dr*cvlight + phi
    if (nsra==1) u = dr
    x = -dr/r
    y = -fllp1*x*x/u + phi
    det = r83sq - x*x + u*y
    b1 = g(kp1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
    b2 = f(kp1)*h83 + r1*df1 + r2*df2 + r3*df3
    g(k) = (b1*(h83-x)+b2*u)/det
    f(k) = (b2*(h83+x)-b1*y)/det
    sg = dsign(1.0d0, g(k))
    sgp1 = dsign(1.0d0, g(kp1))
    if (sg*sgp1<0.d0) nne = nne + 1
    dg1 = dg2
    df1 = df2
    dg2 = dg3
    df2 = df3
    dg3 = u*f(k) - x*g(k)
    df3 = x*f(k) - y*g(k)
    sdg3 = dsign(1.0d0, dg3)
    if (k>k2 .and. sg*sdg3<0.d0) go to 100
  end if
110 kc = k
  dg = dg3
end subroutine
