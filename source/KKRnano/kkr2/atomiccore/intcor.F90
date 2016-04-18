      subroutine intcor(f1, f2, rho, g, f, v, value, slope, l, nn, e, sm, nre, vlnc, a, b, z, rn, nr, tol, irm, ipr, nitmax, nsra)
      implicit none
      external :: hankel, intin, intout

      double precision, intent(in) :: a, b, f1, f2, rn, slope, tol, value, z
      double precision, intent(inout) :: e
      double precision, intent(out) :: sm
      integer, intent(in) :: ipr, irm, l, nitmax, nn, nr, nsra
      integer, intent(out) :: nre
      logical, intent(in) :: vlnc
      double precision, intent(in) :: v(*)
      double precision, intent(out) :: f(*), g(*), rho(*)

!     .. locals ..
      double precision :: cvlight, de, dg1, dg2, dpsi1, dpsi2, drdikc, e1, e2, ea
      double precision :: gkc2, pi, pkc1, pkc2, psi1, psi2, q, qkc1, qkc2, ratio
      double precision :: re, rkc, rpb, slop, tsrme, valu, vme, xxx
      integer :: k, k2, kc, iiter, nne, run
      double complex :: arg, cappai, dofe, hl(6)

      character(len=*), parameter :: F9020="(2i3, 2i4, 1p, 3d16.8, 1x, 2d9.2)", &
      F9000="(' l=', i3, '  nn=', i2, '  nr=', i4, '  f1/e/f2=', 3f10.3, /, ' tol=', 1p, d12.3, '  value/slope=', 2d12.3)", &
      F9030="(/, ' **** int: 0-pressure bcs not real')", F9050="(' *** int: stop after', i4, ' iterations')", &
      F9010="(13x, '  no boundary condition had to be used')", &
      F9040="(' state', i2, ', ', i1, ':', i4, 'x, ', i5, '/', i3, ',  bc=', 1p, 2d12.3, /, 14x, 'e=', d14.6, '   de=', d11.2, '   sm=', d12.4)"

      pi = 4.d0*atan(1.d0)
      cvlight = 274.0720442d0; if (nsra == 1) cvlight = 1.d0
      ea = exp(a)
      iiter = 0
      e1 = f1
      e2 = f2
      rho(1:irm) = 0.d0

      if (ipr == 2) write(6, fmt=F9000) l, nn, nr, f1, e, f2, tol, value, slope

      run = 1
      do while (run > 0)

        iiter = iiter + 1
        if (iiter > nitmax) then
          write(6, fmt=F9050) nitmax
          stop 'intcor'
        endif

        g(1:irm) = 0.d0
        f(1:irm) = 0.d0

        if (e <= e1 .or. e >= e2) e = 0.5d0*(e1 + e2)
        nre = nr
  !     write(6, *) e, e1, e2
        if (e > -1.d-8) return

        tsrme = 2.d0*sqrt(-e)
        re = (log(-tsrme*e/1.d-8)/tsrme - (z+z)/e)*2.d0

        nre = log(re/b + 1.d0)/a + 1.d0
        nre = (nre/2)*2 + 1
        nre = min0(nre, nr)
        nre = max0(nre, 35)

        xxx = 1.d0
        valu = 1.d-1
        slop = -1.d-1
        if (nre < nr .and. iiter == 1 .and. ipr /= 0) write(6, fmt=F9010)
        if (nre >= nr) then
          valu = value
          slop = slope
          if (.not. vlnc) then
  !--->     single site boundary condition
            vme = -e
            if (nsra == 1) then
              cappai = dcmplx(0.d0, dsqrt(vme))
            else
              cappai = dcmplx(0.d0, dsqrt((1.d0 - vme/cvlight/cvlight)*vme))
            endif
            arg = cappai*rn
            call hankel(hl, l+2, arg)
            dofe = real(l+1)/rn - cappai*hl(l+2)/hl(l+1)
            valu = 1.d-10
            slop = valu*dofe
            ! todo: valu and slop are real, maybe this part can be simplified ...
            ! since the Hankel functions are evaluated for purely imaginary arguments
          endif

        endif
        k2 = 30; if (nn == 0) k2 = nre/3
        nne = 0 ! init the number of nodes

        call intin(g, f, v, e, l, nne, valu, slop, nre, k2, kc, dg2, a, b, z, nsra)

        rkc = b*exp(a*kc - a) - b
        drdikc = a*(rkc + b)
        gkc2 = g(kc)
        psi2 = g(kc)
        dpsi2 = dg2/drdikc
        qkc2 = psi2*psi2 + dpsi2*dpsi2*rkc*rkc
        pkc2 = 0.5d0 - atan(rkc*dpsi2/psi2)/pi

        call intout(g, f, v, e, l, nne, kc, dg1, a, b, z, nsra)

        psi1 = g(kc)
        dpsi1 = dg1/drdikc
        qkc1 = psi1*psi1 + dpsi1*dpsi1*rkc*rkc
        pkc1 = 0.5d0 - atan(rkc*dpsi1/psi1)/pi
        if (nne == 9) nne = 0 ! why?
        if (nne == nn) then
  !         ratio1 = gkc2/g(kc)
  !         ratio = sqrt(qkc2/qkc1)
  !         if (ratio1 < 0.d0) ratio = -ratio
          ratio = sign(sqrt(qkc2/qkc1), gkc2*g(kc))
          g(1:kc) = g(1:kc)*ratio
          if (nsra == 1) then
            f(1:nre) = 0.d0
          else
            f(1:kc) = f(1:kc)*ratio
          endif
          sm = 0.d0
          rpb = b/ea
          q = ea*ea
          do k = 2, nre - 1, 2
            rpb = rpb*q
            sm = sm + rpb*(g(k)*g(k) + f(k)*f(k))
          enddo ! k
          rpb = b
          sm = sm + sm
          do k = 3, nre - 2, 2
            rpb = rpb*q
            sm = sm + rpb*(g(k)*g(k) + f(k)*f(k))
          enddo ! k
          sm = 2.d0*sm + rpb*q*(g(nre)*g(nre) + f(nre)*f(nre))
          sm = a*sm/3.d0
          de = pi*qkc2*(pkc2 - pkc1)/sm/rkc
          if (iiter >= nitmax-10 .or. ipr == 2) write(6, fmt=F9020) iiter, nne, nre, kc, e1, e, e2, de
          if (de > 0.d0) e1 = e
          if (de < 0.d0) e2 = e
          e = e + de
          if (abs(de) <= tol .or. iiter >= nitmax) run = 0 ! converged

        else
          if (iiter >= nitmax-10 .or. ipr == 2) write(6, fmt=F9020) iiter, nne, nre, kc, e1, e, e2
          if (nne > nn) e2 = e
          if (nne < nn) e1 = e
          e = .5d0* (e1+e2)

        endif
      enddo ! while(run > 0)

      e = e - de
      do k = 1, nre
        rho(k) = g(k)*g(k) + f(k)*f(k)
      enddo ! k
      if (xxx <= 0.d0) write(6, fmt=F9030)
      if (iiter >= nitmax-10 .or. ipr >= 1 .or. xxx <= 0.d0) &
        write(6, fmt=F9040) l, nn, iiter, kc, nre, valu, slop, e, de, sm
      endsubroutine intcor
