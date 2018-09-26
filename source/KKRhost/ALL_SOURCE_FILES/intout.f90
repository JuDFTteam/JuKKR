module mod_intout
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  subroutine intout(g, f, v, e, l, nne, k2, dg, a, b, z, nsra)
    implicit none
    ! .. Scalar Arguments ..
    real (kind=dp) :: a, b, dg, e, z
    integer :: k2, l, nne, nsra
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: f(*), g(*), v(*)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: aa, alfa, b1, b2, bb, beta, cvlight, det, df1, df2, df3, dg1, dg2, dg3, dr, ea, fllp1, h83, p12, p21, phi, pp, qq, r, r1, r2, r3, r83sq, rpb, s, sg, sgm1, u, &
      x, y, zz
    integer :: i, i1, k, km1, n
    ! ..
    ! .. Local Arrays ..
    real (kind=dp) :: d(2, 3), px(20), qx(20)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: sign, exp, sqrt
    ! ..
    ! ..
    zz = z + z
    cvlight = 274.0720442e0_dp
    if (nsra==1) cvlight = 1.0e0_dp
    ea = exp(a)
    fllp1 = l*(l+1.e0_dp)
    r83sq = 64.e0_dp/9.e0_dp
    r1 = 1.e0_dp/9.e0_dp
    r2 = -5.e0_dp*r1
    r3 = 19.e0_dp*r1
    h83 = 8.e0_dp/3.e0_dp
    aa = -zz/cvlight
    bb = fllp1 - aa*aa
    p21 = (v(1)-e)/cvlight
    p12 = cvlight - p21
    px(1) = 0.e0_dp
    qx(1) = 0.e0_dp
    if (z<=20.e0_dp .or. nsra==1) then
      s = 1.e0_dp*l
      px(2) = 0.e0_dp
      px(3) = 1.e0_dp
      do k = 2, 9
        px(k+2) = ((v(1)-e)*px(k)-zz*px(k+1))/(k+l+l)/(k-1.e0_dp)
      end do
      do k = 2, 10
        qx(k) = px(k+1)*(l+k-2.e0_dp)/cvlight
      end do

    else
      s = sqrt(fllp1+1.e0_dp-aa*aa)
      px(2) = 1.e0_dp
      qx(2) = (1.e0_dp-s)/aa
      do i = 3, 10
        n = i - 2
        alfa = p12*qx(i-1)
        beta = p21*aa*px(i-1)
        if (l/=0) beta = beta - p12*aa*px(i-1) - p12*p21*px(i-2) + (n+s)*p12*qx(i-1)
        det = n*(n+s+s)*aa
        px(i) = (alfa*(n+s+1.e0_dp)*aa-aa*beta)/det
        qx(i) = (beta*(n+s-1.e0_dp)-bb*alfa)/det
      end do
    end if

    g(1) = 0.e0_dp
    f(1) = 0.e0_dp
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
      sg = sign(1.0e0_dp, g(k))
      sgm1 = sign(1.0e0_dp, g(k-1))
      if (sg*sgm1<0.e0_dp) nne = nne + 1
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
      sg = sign(1.0e0_dp, g(k))
      sgm1 = sign(1.0e0_dp, g(km1))
      if (sg*sgm1<0.e0_dp) nne = nne + 1
      dg1 = dg2
      dg2 = dg3
      dg3 = u*f(k) - x*g(k)
      df1 = df2
      df2 = df3
      df3 = x*f(k) - y*g(k)
    end do
    dg = dg3
  end subroutine intout

end module mod_intout
