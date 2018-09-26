module mod_intcor

contains

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  subroutine intcor(f1, f2, rho, g, f, v, value, slope, l, nn, e, sum, nre, vlnc, a, b, z, rn, nr, tol, irm, ipr, nitmax, nsra)
    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    use :: mod_intin, only: intin
    use :: mod_intout, only: intout
    use :: mod_hankel, only: hankel
    implicit none
    ! .. Scalar Arguments ..
    real (kind=dp) :: a, b, e, f1, f2, rn, slope, sum, tol, value, z
    integer :: ipr, irm, l, nitmax, nn, nr, nre, nsra
    logical :: vlnc
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: f(*), g(*), rho(*), v(*)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: arg, cappai, dofe
    real (kind=dp) :: cvlight, de, dg1, dg2, dpsi1, dpsi2, drdikc, e1, e2, ea, gkc2, pi, pkc1, pkc2, psi1, psi2, q, qkc1, qkc2, ratio, ratio1, re, rkc, rpb, slop, tsrme, valu, vme, &
      xxx, zz
    integer :: ir, k, k2, kc, niter, nne, nrem1, nrem2
    ! ..
    ! .. Local Arrays ..
    complex (kind=dp) :: hl(6)
    ! ..
    pi = 4.d0*atan(1.d0)
    zz = z + z
    cvlight = 274.0720442d0
    if (nsra==1) cvlight = 1.0d0
    ea = exp(a)
    niter = 0
    e1 = f1
    e2 = f2
    if (ipr==2 .and. (t_inc%i_write>0)) write (1337, fmt=120) l, nn, nr, f1, e, f2, tol, value, slope
100 continue
    niter = niter + 1
    do ir = 1, irm
      g(ir) = 0.0d0
      f(ir) = 0.0d0
    end do
    if (niter>nitmax) then
      go to 110

    else
      if (e<=e1 .or. e>=e2) e = .5d0*(e1+e2)
      nre = nr
      if (e<=-1.d-8) then
        tsrme = 2.d0*sqrt(-e)
        re = (log(-tsrme*e/1.d-8)/tsrme-zz/e)*2.d0
        nre = nint(log(re/b+1.d0)/a+1.d0)
        nre = (nre/2)*2 + 1
        nre = min0(nre, nr)
        nre = max0(nre, 35)
      end if
      xxx = 1.d0
      valu = 1.d-1
      slop = -1.d-1
      if (nre<nr .and. niter==1 .and. ipr/=0 .and. (t_inc%i_write>0)) write (1337, fmt=130)
      if (nre>=nr) then
        valu = value
        slop = slope
        if (.not. vlnc) then
          ! --->   single site  boundary condition
          vme = -e
          if (nsra==1) then
            cappai = cmplx(0.d0, sqrt(vme), kind=dp)
          else
            cappai = cmplx(0.d0, sqrt((1.d0-vme/cvlight/cvlight)*vme), kind=dp)
          end if
          arg = cappai*rn
          call hankel(hl, l+2, arg)
          dofe = real(l+1, kind=dp)/rn - cappai*hl(l+2)/hl(l+1)
          valu = 1.d-10
          slop = real(valu*dofe, kind=dp)
        end if

      end if
      k2 = 30
      if (nn==0) k2 = nre/3
      nne = 0

      call intin(g, f, v, e, l, nne, valu, slop, nre, k2, kc, dg2, a, b, z, nsra)

      rkc = b*exp(a*kc-a) - b
      drdikc = a*(rkc+b)
      gkc2 = g(kc)
      psi2 = g(kc)
      dpsi2 = dg2/drdikc
      qkc2 = psi2*psi2 + dpsi2*dpsi2*rkc*rkc
      pkc2 = .5d0 - atan(rkc*dpsi2/psi2)/pi

      call intout(g, f, v, e, l, nne, kc, dg1, a, b, z, nsra)

      psi1 = g(kc)
      dpsi1 = dg1/drdikc
      qkc1 = psi1*psi1 + dpsi1*dpsi1*rkc*rkc
      pkc1 = .5d0 - atan(rkc*dpsi1/psi1)/pi
      if (nne==9) nne = 0
      if (nne==nn) then
        ratio1 = gkc2/g(kc)
        ratio = sqrt(qkc2/qkc1)
        if (ratio1<0.d0) ratio = -ratio
        do k = 1, kc
          g(k) = g(k)*ratio
          f(k) = f(k)*ratio
        end do
        sum = 0.d0
        if (nsra==1) then
          do k = 1, nre
            f(k) = 0.0d0
          end do
        end if
        rpb = b/ea
        q = ea*ea
        nrem1 = nre - 1
        do k = 2, nrem1, 2
          rpb = rpb*q
          sum = sum + rpb*(g(k)*g(k)+f(k)*f(k))
        end do
        rpb = b
        sum = sum + sum
        nrem2 = nre - 2
        do k = 3, nrem2, 2
          rpb = rpb*q
          sum = sum + rpb*(g(k)*g(k)+f(k)*f(k))
        end do
        sum = sum + sum + rpb*q*(g(nre)*g(nre)+f(nre)*f(nre))
        sum = a*sum/3.d0
        de = pi*qkc2*(pkc2-pkc1)/sum/rkc
        if (niter>=nitmax-10 .or. ipr==2 .and. (t_inc%i_write>0)) write (1337, fmt=140) niter, nne, nre, kc, e1, e, e2, de
        if (de>0.d0) e1 = e
        if (de<0.d0) e2 = e
        e = e + de
        if (abs(de)>tol .and. niter<nitmax) go to 100

      else
        if (niter>=nitmax-10 .or. ipr==2 .and. (t_inc%i_write>0)) write (1337, fmt=140) niter, nne, nre, kc, e1, e, e2
        if (nne>nn) e2 = e
        if (nne<nn) e1 = e
        e = .5d0*(e1+e2)
        go to 100

      end if

    end if

    e = e - de
    do k = 1, nre
      rho(k) = g(k)*g(k) + f(k)*f(k)
    end do
    if (xxx<=0.d0 .and. (t_inc%i_write>0)) write (1337, fmt=150)
    if (niter>=nitmax-10 .or. ipr>=1 .or. xxx<=0.d0 .and. (t_inc%i_write>0)) write (1337, fmt=160) l, nn, niter, kc, nre, valu, slop, e, de, sum
    return

110 write (6, fmt=170) nitmax
    stop 'INTCOR'


120 format (' l=', i3, '  nn=', i2, '  nr=', i4, '  f1/e/f2=', 3f10.3, /, ' tol=', 1p, d12.3, '  value/slope=', 2d12.3)
130 format (13x, '  no boundary condition had to be used')
140 format (2i3, 2i4, 1p, 3d16.8, 1x, 2d9.2)
150 format (/, ' **** int: 0-pressure bcs not real')
160 format (' state', i2, ',', i1, ':', i4, 'x,', i5, '/', i3, ',  bc=', 1p, 2d12.3, /, 14x, 'e=', d14.6, '   de=', d11.2, '   sum=', d12.4)
170 format (' *** int: stop after', i4, ' iterations')
  end subroutine intcor

end module mod_intcor
