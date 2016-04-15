      subroutine intin(g, f, v, e, l, nne, valu, slop, k1, k2, kc, dg, a, b, z, nsra)
      implicit none
      integer, intent(in) :: k1, k2, l, nsra
      integer, intent(inout) :: nne
      integer, intent(out) :: kc
      double precision, intent(in) :: a, b, e, slop, valu, z
      double precision, intent(in) :: v(*)
      double precision, intent(out) :: f(*),  g(*), dg

!     .. locals ..
      double precision :: af1, af2, af3, ag1, ag2, ag3, b1, b2, cvlight, det, df1
      double precision :: df2, df3, dg1, dg2, dg3, dr, ea, ff, fllp1, gg, phi, q
      double precision :: r, rpb, sdg3, sg, sgp1, u, vb, x, y, zz
      double precision, parameter :: r83sq=64.d0/9.d0, r1=1.d0/9.d0, r2=-5.d0*r1, r3=19.d0*r1, h83=-8.d0/3.d0
      integer :: i, ir, jr
      double precision :: d(2,3)

      zz = z + z
      cvlight = 274.0720442d0; if (nsra == 1) cvlight = 1.d0
      fllp1 = l*(l + 1.d0)
      
      ea = exp(a)
      rpb = b*exp(a*k1 - a)
      r = rpb - b
      dr = a*rpb
      phi = (e + zz/r - v(k1))*dr/cvlight
      u = dr*cvlight + phi
      if (nsra == 1) u = dr
      x = -dr/r
      y = -fllp1*x*x/u + phi
      g(k1) = valu
      f(k1) = (slop*dr + x*valu)/u
      q = 1.d0/sqrt(ea)
      ag1 = slop*dr
      af1 = x*f(k1) - y*g(k1)
      ir = k1
      dg3 = ag1
      if (k2 /= k1) then
        do i = 1, 3
          jr = ir
          ir = ir - 1
          rpb = rpb*q
          dr = rpb*a
          r = rpb - b
          gg = g(jr) - .5d0*ag1
          ff = f(jr) - .5d0*af1
          vb = (3.d0*v(jr) + 6.d0*v(ir) - v(ir-1))*.125d0
          phi = (e + zz/r - vb)*dr/cvlight
          u = dr*cvlight + phi
          if (nsra == 1) u = dr
          x = -dr/r
          y = -fllp1*x*x/u + phi
          ag2 = u*ff - x*gg
          af2 = x*ff - y*gg
          gg = g(jr) - .5d0*ag2
          ff = f(jr) - .5d0*af2
          ag3 = u*ff - x*gg
          af3 = x*ff - y*gg
          rpb = rpb*q
          dr = a*rpb
          r = rpb - b
          phi = (e + zz/r - v(ir))*dr/cvlight
          u = dr*cvlight + phi
          if (nsra == 1) u = dr
          x = -dr/r
          y = -fllp1*x*x/u + phi
          gg = g(jr) - ag3
          ff = f(jr) - af3
          g(ir) = g(jr) - (ag1 + 2.d0*(ag2 + ag3) + u*ff - x*gg)/6.d0
          f(ir) = f(jr) - (af1 + 2.d0*(af2 + af3) + x*ff - y*gg)/6.d0
          sg   = dsign(1.d0,g(ir))
          sgp1 = dsign(1.d0,g(jr))
          if (sg*sgp1 < 0.d0) nne = nne + 1
          ag1 = u*f(ir) - x*g(ir)
          af1 = x*f(ir) - y*g(ir)
          if (ir == k2) then
            kc = ir
            dg = dg3
            return
          else
            d(1,i) = ag1
            d(2,i) = af1
          endif
        enddo
        q = 1.d0/ea
        dg1 = d(1,1)
        dg2 = d(1,2)
        dg3 = d(1,3)
        df1 = d(2,1)
        df2 = d(2,2)
        df3 = d(2,3)
   20   continue
        jr = ir
        ir = jr - 1
        rpb = rpb*q
        dr = a*rpb
        r = rpb - b
        phi = (e + zz/r - v(ir))*dr/cvlight
        u = dr*cvlight + phi
        if (nsra == 1) u = dr
        x = -dr/r
        y = -fllp1*x*x/u + phi
        det = r83sq - x*x + u*y
        b1 = g(jr)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b2 = f(jr)*h83 + r1*df1 + r2*df2 + r3*df3
        g(ir) = (b1*(h83 - x) + b2*u)/det
        f(ir) = (b2*(h83 + x) - b1*y)/det
        sg   = dsign(1.d0,g(ir))
        sgp1 = dsign(1.d0,g(jr))
        if (sg*sgp1 < 0.d0) nne = nne + 1
        dg1 = dg2
        df1 = df2
        dg2 = dg3
        df2 = df3
        dg3 = u*f(ir) - x*g(ir)
        df3 = x*f(ir) - y*g(ir)
        sdg3 = dsign(1.d0,dg3)
        if (ir > k2 .and. sg*sdg3 < 0.d0) goto 20
      endif
      kc = ir
      dg = dg3
      endsubroutine intin
      
