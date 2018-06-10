subroutine intout(g, f, v, e, l, nne, k2, dg, a, b, z, nsra)
  implicit none
!.. Scalar Arguments ..
  double precision :: a, b, dg, e, z
  integer :: k2, l, nne, nsra
!..
!.. Array Arguments ..
  double precision :: f(*), g(*), v(*)
!..
!.. Local Scalars ..
  double precision :: aa, alfa, b1, b2, bb, beta, cvlight, det, df1, df2, df3, &
    dg1, dg2, dg3, dr, ea, fllp1, h83, p12, p21, phi, pp, qq, r, r1, r2, r3, &
    r83sq, rpb, s, sg, sgm1, u, x, y, zz
  integer :: i, i1, k, km1, n
!..
!.. Local Arrays ..
  double precision :: d(2, 3), px(20), qx(20)
!..
!.. Intrinsic Functions ..
  intrinsic :: dsign, exp, sqrt
!..
!     ..
  zz = z + z
  cvlight = 274.0720442d0
  if (nsra==1) cvlight = 1.0d0
  ea = exp(a)
  fllp1 = l*(l+1.d0)
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
  if (z<=20.d0 .or. nsra==1) then
    s = 1.d0*l
    px(2) = 0.d0
    px(3) = 1.d0
    do k = 2, 9
      px(k+2) = ((v(1)-e)*px(k)-zz*px(k+1))/(k+l+l)/(k-1.d0)
    end do
    do k = 2, 10
      qx(k) = px(k+1)*(l+k-2.d0)/cvlight
    end do

  else
    s = sqrt(fllp1+1.d0-aa*aa)
    px(2) = 1.d0
    qx(2) = (1.d0-s)/aa
    do i = 3, 10
      n = i - 2
      alfa = p12*qx(i-1)
      beta = p21*aa*px(i-1)
      if (l/=0) beta = beta - p12*aa*px(i-1) - p12*p21*px(i-2) + &
        (n+s)*p12*qx(i-1)
      det = n*(n+s+s)*aa
      px(i) = (alfa*(n+s+1.d0)*aa-aa*beta)/det
      qx(i) = (beta*(n+s-1.d0)-bb*alfa)/det
    end do
  end if

  g(1) = 0.d0
  f(1) = 0.d0
  rpb = b
  do k = 2, 4
    rpb = rpb*ea
    r = rpb - b
    dr = a*rpb
    phi = (e+zz/r-v(k))*dr/cvlight
    u = dr*cvlight + phi
    if (nsra==1) u = dr
    x = -dr/r
    y = -fllp1*x*x/u + phi
    pp = px(10)
    qq = qx(10)
    do i1 = 3, 10
      i = 12 - i1
      pp = pp*r + px(i)
      qq = qq*r + qx(i)
    end do
    g(k) = (r**s)*pp
    f(k) = (r**s)*qq
    sg = dsign(1.0d0, g(k))
    sgm1 = dsign(1.0d0, g(k-1))
    if (sg*sgm1<0.d0) nne = nne + 1
    d(1, k-1) = u*f(k) - x*g(k)
    d(2, k-1) = x*f(k) - y*g(k)
  end do
  dg1 = d(1, 1)
  dg2 = d(1, 2)
  dg3 = d(1, 3)
  df1 = d(2, 1)
  df2 = d(2, 2)
  df3 = d(2, 3)
  do k = 5, k2
    km1 = k - 1
    rpb = rpb*ea
    r = rpb - b
    dr = a*rpb
    phi = (e+zz/r-v(k))*dr/cvlight
    u = dr*cvlight + phi
    if (nsra==1) u = dr
    x = -dr/r
    y = -fllp1*x*x/u + phi
    det = r83sq - x*x + u*y
    b1 = g(km1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
    b2 = f(km1)*h83 + r1*df1 + r2*df2 + r3*df3
    g(k) = (b1*(h83-x)+b2*u)/det
    f(k) = (b2*(h83+x)-b1*y)/det
    sg = dsign(1.0d0, g(k))
    sgm1 = dsign(1.0d0, g(km1))
    if (sg*sgm1<0.d0) nne = nne + 1
    dg1 = dg2
    dg2 = dg3
    dg3 = u*f(k) - x*g(k)
    df1 = df2
    df2 = df3
    df3 = x*f(k) - y*g(k)
  end do
  dg = dg3
end subroutine
