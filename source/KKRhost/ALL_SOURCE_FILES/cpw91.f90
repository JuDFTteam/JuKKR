module mod_cpw91
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Correlation of PW91
  !> Author: 
  !> Category: KKRhost, xc-potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> gga91 correlation
  !> 
  !> Input
  !>   rs: seitz radius
  !>   zta: relative spin polarization
  !>   t: abs(grad d)/(d*2.*ks*gz)
  !>   uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*gz)**3)
  !>   vv: (laplacian d)/(d * (2*ks*gz)**2)
  !>   ww: (grad d)*(gradzta)/(d * (2*ks*gz)**2
  !> Output
  !>   h: nonlocal part of correlation energy per electron
  !>   dvcup,-dn: nonlocal parts of correlation potentials.
  !>
  !> with ks=sqrt(4*kf/pai), gz=[(1+zta)**(2/3)+(1-zta)**(2/3)]/2, &
  !> kf=cbrt(3*pai**2*d).
  !-------------------------------------------------------------------------------
  subroutine cpw91(fk, sk, gz, ec, ecrs, eczta, rs, zta, t, uu, vv, ww, h, dvcup, dvcdn)

    implicit none
    ! .. Scalar Arguments ..
    real (kind=dp) :: dvcdn, dvcup, ec, ecrs, eczta, fk, gz, h, rs, sk, t, uu, vv, ww, zta
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: a4, alf, argmx, b, b2, b2fac, bec, bet, bg, c1, c2, c3, c4, c5, c6, cc, cc0, ccrs, coeff, comm, cx, delt, fact0, fact1, fact2, fact3, fact4, fact5, gm, gz3, &
      gz4, h0, h0b, h0bt, h0rs, h0rst, h0t, h0tt, h0z, h0zt, h1, h1rs, h1rst, h1t, h1tt, h1z, h1zt, hrs, hrst, ht, htt, hz, hzt, pon, pref, q4, q5, q6, q7, q8, q9, r0, r1, r2, r3, &
      r4, rs2, rs3, rsthrd, t2, t4, t6, thrd2, thrdm, xnu
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: exp, log
    ! ..
    ! .. Save statement ..
    save :: xnu, cc0, cx, alf, c1, c2, c3, c4, c5, c6, a4, thrdm, thrd2
    ! ..
    ! .. Data statements ..
    ! -----------------------------------------------------------------
    data xnu, cc0, cx, alf/15.75592e0_dp, 0.004235e0_dp, -0.001667212e0_dp, 0.09e0_dp/
    data c1, c2, c3, c4/0.002568e0_dp, 0.023266e0_dp, 7.389e-6_dp, 8.723e0_dp/
    data c5, c6, a4/0.472e0_dp, 7.389e-2_dp, 100.e0_dp/
    data thrdm, thrd2/ -0.333333333333e0_dp, 0.666666666667e0_dp/
    ! ..
    ! -----------------------------------------------------------------
    argmx = 174.0e0_dp
    bet = xnu*cc0
    delt = 2.e0_dp*alf/bet
    gz3 = gz**3
    gz4 = gz3*gz
    pon = -delt*ec/(gz3*bet)
    if (pon>argmx) then
      b = 0.e0_dp
    else
      b = delt/(exp(pon)-1.e0_dp)
    end if
    b2 = b*b
    t2 = t*t
    t4 = t2*t2
    t6 = t4*t2
    rs2 = rs*rs
    rs3 = rs2*rs
    q4 = 1.e0_dp + b*t2
    q5 = 1.e0_dp + b*t2 + b2*t4
    q6 = c1 + c2*rs + c3*rs2
    q7 = 1.e0_dp + c4*rs + c5*rs2 + c6*rs3
    cc = -cx + q6/q7
    r0 = (sk/fk)**2
    r1 = a4*r0*gz4
    coeff = cc - cc0 - 3.e0_dp*cx/7.e0_dp
    r2 = xnu*coeff*gz3
    r3 = exp(-r1*t2)
    h0 = gz3*(bet/delt)*log(1.e0_dp+delt*q4*t2/q5)
    h1 = r3*r2*t2
    h = h0 + h1
    ! local correlation option:
    ! h = 0.0d0

    ! energy done. now the potential:

    ccrs = (c2+2.e0_dp*c3*rs)/q7 - q6*(c4+2.e0_dp*c5*rs+3.e0_dp*c6*rs2)/q7**2
    rsthrd = rs/3.e0_dp
    r4 = rsthrd*ccrs/coeff
    gm = ((1.e0_dp+zta)**thrdm-(1.e0_dp-zta)**thrdm)/3.e0_dp
    if (pon>argmx) then
      b2fac = 0.e0_dp
    else
      b2fac = b2*(delt/b+1.e0_dp)
    end if
    bg = -3.e0_dp*ec*b2fac/(bet*gz4)
    bec = b2fac/(bet*gz3)
    q8 = q5*q5 + delt*q4*q5*t2
    q9 = 1.e0_dp + 2.e0_dp*b*t2
    h0b = -bet*gz3*b*t6*(2.e0_dp+b*t2)/q8
    h0rs = -rsthrd*h0b*bec*ecrs
    fact0 = 2.e0_dp*delt - 6.e0_dp*b
    fact1 = q5*q9 + q4*q9*q9
    h0bt = 2.e0_dp*bet*gz3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
    h0rst = rsthrd*t2*h0bt*bec*ecrs
    h0z = 3.e0_dp*gm*h0/gz + h0b*(bg*gm+bec*eczta)
    h0t = 2.e0_dp*bet*gz3*q9/q8
    h0zt = 3.e0_dp*gm*h0t/gz + h0bt*(bg*gm+bec*eczta)
    fact2 = q4*q5 + b*t2*(q4*q9+q5)
    fact3 = 2.e0_dp*b*q5*q9 + delt*fact2
    h0tt = 4.e0_dp*bet*gz3*t*(2.e0_dp*b/q8-(q9*(fact3/q8))/q8)
    h1rs = r3*r2*t2*(-r4+r1*t2/3.e0_dp)
    fact4 = 2.e0_dp - r1*t2
    h1rst = r3*r2*t2*(2.e0_dp*r4*(1.e0_dp-r1*t2)-thrd2*r1*t2*fact4)
    h1z = gm*r3*r2*t2*(3.e0_dp-4.e0_dp*r1*t2)/gz
    h1t = 2.e0_dp*r3*r2*(1.e0_dp-r1*t2)
    h1zt = 2.e0_dp*gm*r3*r2*(3.e0_dp-11.e0_dp*r1*t2+4.e0_dp*r1*r1*t4)/gz
    h1tt = 4.e0_dp*r3*r2*r1*t*(-2.e0_dp+r1*t2)
    hrs = h0rs + h1rs
    hrst = h0rst + h1rst
    ht = h0t + h1t
    htt = h0tt + h1tt
    hz = h0z + h1z
    hzt = h0zt + h1zt
    comm = h + hrs + hrst + t2*ht/6.e0_dp + 7.e0_dp*t2*t*htt/6.e0_dp
    pref = hz - gm*t2*ht/gz
    fact5 = gm*(2.e0_dp*ht+t*htt)/gz
    comm = comm - pref*zta - uu*htt - vv*ht - ww*(hzt-fact5)
    dvcup = comm + pref
    dvcdn = comm - pref

    ! local correlation option:
    ! dvcup = 0.0d0
    ! dvcdn = 0.0d0

    return
  end subroutine cpw91

end module mod_cpw91
