    Subroutine intin(g, f, v, e, l, nne, valu, slop, k1, k2, kc, dg, a, b, z, &
      nsra)
      Use mod_datatypes, Only: dp
!.. Scalar Arguments ..
      Real (Kind=dp) :: a, b, dg, e, slop, valu, z
      Integer :: k1, k2, kc, l, nne, nsra
!..
!.. Array Arguments ..
      Real (Kind=dp) :: f(*), g(*), v(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: af1, af2, af3, ag1, ag2, ag3, b1, b2, cvlight, det, &
        df1, df2, df3, dg1, dg2, dg3, dr, ea, ff, fllp1, gg, h83, phi, q, r, &
        r1, r2, r3, r83sq, rpb, sdg3, sg, sgp1, u, vb, x, y, zz
      Integer :: i, k, kp1
!..
!.. Local Arrays ..
      Real (Kind=dp) :: d(2, 3)
!..
!.. Intrinsic Functions ..
      Intrinsic :: sign, exp, sqrt
!     ..
      zz = z + z
      cvlight = 274.0720442E0_dp
      If (nsra==1) cvlight = 1.0E0_dp
      fllp1 = l*(l+1.E0_dp)
      r83sq = 64.E0_dp/9.E0_dp
      r1 = 1.E0_dp/9.E0_dp
      r2 = -5.E0_dp*r1
      r3 = 19.E0_dp*r1
      h83 = -8.E0_dp/3.E0_dp
      ea = exp(a)
      rpb = b*exp(a*k1-a)
      r = rpb - b
      dr = a*rpb
      phi = (e+zz/r-v(k1))*dr/cvlight
      u = dr*cvlight + phi
      If (nsra==1) u = dr
      x = -dr/r
      y = -fllp1*x*x/u + phi
      g(k1) = valu
      f(k1) = (slop*dr+x*valu)/u
      q = 1.E0_dp/sqrt(ea)
      ag1 = slop*dr
      af1 = x*f(k1) - y*g(k1)
      k = k1
      dg3 = ag1
      If (k2/=k1) Then
        Do i = 1, 3
          kp1 = k
          k = k - 1
          rpb = rpb*q
          dr = rpb*a
          r = rpb - b
          gg = g(kp1) - .5E0_dp*ag1
          ff = f(kp1) - .5E0_dp*af1
          vb = (3.E0_dp*v(kp1)+6.E0_dp*v(k)-v(k-1))*.125E0_dp
          phi = (e+zz/r-vb)*dr/cvlight
          u = dr*cvlight + phi
          If (nsra==1) u = dr
          x = -dr/r
          y = -fllp1*x*x/u + phi
          ag2 = u*ff - x*gg
          af2 = x*ff - y*gg
          gg = g(kp1) - .5E0_dp*ag2
          ff = f(kp1) - .5E0_dp*af2
          ag3 = u*ff - x*gg
          af3 = x*ff - y*gg
          rpb = rpb*q
          dr = a*rpb
          r = rpb - b
          phi = (e+zz/r-v(k))*dr/cvlight
          u = dr*cvlight + phi
          If (nsra==1) u = dr
          x = -dr/r
          y = -fllp1*x*x/u + phi
          gg = g(kp1) - ag3
          ff = f(kp1) - af3
          g(k) = g(kp1) - (ag1+2.E0_dp*(ag2+ag3)+u*ff-x*gg)/6.E0_dp
          f(k) = f(kp1) - (af1+2.E0_dp*(af2+af3)+x*ff-y*gg)/6.E0_dp
          sg = sign(1.0E0_dp, g(k))
          sgp1 = sign(1.0E0_dp, g(kp1))
          If (sg*sgp1<0.E0_dp) nne = nne + 1
          ag1 = u*f(k) - x*g(k)
          af1 = x*f(k) - y*g(k)
          If (k==k2) Then
            Go To 110

          Else
            d(1, i) = ag1
            d(2, i) = af1
          End If

        End Do
        q = 1.E0_dp/ea
        dg1 = d(1, 1)
        dg2 = d(1, 2)
        dg3 = d(1, 3)
        df1 = d(2, 1)
        df2 = d(2, 2)
        df3 = d(2, 3)
100     Continue
        kp1 = k
        k = k - 1
        rpb = rpb*q
        dr = a*rpb
        r = rpb - b
        phi = (e+zz/r-v(k))*dr/cvlight
        u = dr*cvlight + phi
        If (nsra==1) u = dr
        x = -dr/r
        y = -fllp1*x*x/u + phi
        det = r83sq - x*x + u*y
        b1 = g(kp1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b2 = f(kp1)*h83 + r1*df1 + r2*df2 + r3*df3
        g(k) = (b1*(h83-x)+b2*u)/det
        f(k) = (b2*(h83+x)-b1*y)/det
        sg = sign(1.0E0_dp, g(k))
        sgp1 = sign(1.0E0_dp, g(kp1))
        If (sg*sgp1<0.E0_dp) nne = nne + 1
        dg1 = dg2
        df1 = df2
        dg2 = dg3
        df2 = df3
        dg3 = u*f(k) - x*g(k)
        df3 = x*f(k) - y*g(k)
        sdg3 = sign(1.0E0_dp, dg3)
        If (k>k2 .And. sg*sdg3<0.E0_dp) Go To 100
      End If
110   kc = k
      dg = dg3
    End Subroutine
