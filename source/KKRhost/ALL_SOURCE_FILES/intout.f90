    Subroutine intout(g, f, v, e, l, nne, k2, dg, a, b, z, nsra)
      Use mod_datatypes, Only: dp
      Implicit None
!.. Scalar Arguments ..
      Real (Kind=dp) :: a, b, dg, e, z
      Integer :: k2, l, nne, nsra
!..
!.. Array Arguments ..
      Real (Kind=dp) :: f(*), g(*), v(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: aa, alfa, b1, b2, bb, beta, cvlight, det, df1, df2, &
        df3, dg1, dg2, dg3, dr, ea, fllp1, h83, p12, p21, phi, pp, qq, r, r1, &
        r2, r3, r83sq, rpb, s, sg, sgm1, u, x, y, zz
      Integer :: i, i1, k, km1, n
!..
!.. Local Arrays ..
      Real (Kind=dp) :: d(2, 3), px(20), qx(20)
!..
!.. Intrinsic Functions ..
      Intrinsic :: sign, exp, sqrt
!..
!     ..
      zz = z + z
      cvlight = 274.0720442E0_dp
      If (nsra==1) cvlight = 1.0E0_dp
      ea = exp(a)
      fllp1 = l*(l+1.E0_dp)
      r83sq = 64.E0_dp/9.E0_dp
      r1 = 1.E0_dp/9.E0_dp
      r2 = -5.E0_dp*r1
      r3 = 19.E0_dp*r1
      h83 = 8.E0_dp/3.E0_dp
      aa = -zz/cvlight
      bb = fllp1 - aa*aa
      p21 = (v(1)-e)/cvlight
      p12 = cvlight - p21
      px(1) = 0.E0_dp
      qx(1) = 0.E0_dp
      If (z<=20.E0_dp .Or. nsra==1) Then
        s = 1.E0_dp*l
        px(2) = 0.E0_dp
        px(3) = 1.E0_dp
        Do k = 2, 9
          px(k+2) = ((v(1)-e)*px(k)-zz*px(k+1))/(k+l+l)/(k-1.E0_dp)
        End Do
        Do k = 2, 10
          qx(k) = px(k+1)*(l+k-2.E0_dp)/cvlight
        End Do

      Else
        s = sqrt(fllp1+1.E0_dp-aa*aa)
        px(2) = 1.E0_dp
        qx(2) = (1.E0_dp-s)/aa
        Do i = 3, 10
          n = i - 2
          alfa = p12*qx(i-1)
          beta = p21*aa*px(i-1)
          If (l/=0) beta = beta - p12*aa*px(i-1) - p12*p21*px(i-2) + &
            (n+s)*p12*qx(i-1)
          det = n*(n+s+s)*aa
          px(i) = (alfa*(n+s+1.E0_dp)*aa-aa*beta)/det
          qx(i) = (beta*(n+s-1.E0_dp)-bb*alfa)/det
        End Do
      End If

      g(1) = 0.E0_dp
      f(1) = 0.E0_dp
      rpb = b
      Do k = 2, 4
        rpb = rpb*ea
        r = rpb - b
        dr = a*rpb
        phi = (e+zz/r-v(k))*dr/cvlight
        u = dr*cvlight + phi
        If (nsra==1) u = dr
        x = -dr/r
        y = -fllp1*x*x/u + phi
        pp = px(10)
        qq = qx(10)
        Do i1 = 3, 10
          i = 12 - i1
          pp = pp*r + px(i)
          qq = qq*r + qx(i)
        End Do
        g(k) = (r**s)*pp
        f(k) = (r**s)*qq
        sg = sign(1.0E0_dp, g(k))
        sgm1 = sign(1.0E0_dp, g(k-1))
        If (sg*sgm1<0.E0_dp) nne = nne + 1
        d(1, k-1) = u*f(k) - x*g(k)
        d(2, k-1) = x*f(k) - y*g(k)
      End Do
      dg1 = d(1, 1)
      dg2 = d(1, 2)
      dg3 = d(1, 3)
      df1 = d(2, 1)
      df2 = d(2, 2)
      df3 = d(2, 3)
      Do k = 5, k2
        km1 = k - 1
        rpb = rpb*ea
        r = rpb - b
        dr = a*rpb
        phi = (e+zz/r-v(k))*dr/cvlight
        u = dr*cvlight + phi
        If (nsra==1) u = dr
        x = -dr/r
        y = -fllp1*x*x/u + phi
        det = r83sq - x*x + u*y
        b1 = g(km1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b2 = f(km1)*h83 + r1*df1 + r2*df2 + r3*df3
        g(k) = (b1*(h83-x)+b2*u)/det
        f(k) = (b2*(h83+x)-b1*y)/det
        sg = sign(1.0E0_dp, g(k))
        sgm1 = sign(1.0E0_dp, g(km1))
        If (sg*sgm1<0.E0_dp) nne = nne + 1
        dg1 = dg2
        dg2 = dg3
        dg3 = u*f(k) - x*g(k)
        df1 = df2
        df2 = df3
        df3 = x*f(k) - y*g(k)
      End Do
      dg = dg3
    End Subroutine
