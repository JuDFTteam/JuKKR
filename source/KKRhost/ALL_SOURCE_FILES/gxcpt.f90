subroutine gxcpt(idspr, ro, zta, agr, agru, agrd, g2r, g2ru, g2rd, gggr, &
  gggru, gggrd, grgru, grgrd, gzgr, xcptu, xcptd, xced, vxlu, vxld, vclu, &
  vcld, xedl, cedl, vxgu, vxgd, vcgu, vcgd, xedg, cedg)
!.....-----------------------------------------------------------------
!.....gxcp: exchange-correlation potential in ry. also total-energy.
!.....-----------------------------------------------------------------
  implicit none
!.. Scalar Arguments ..
  double precision :: agr, agrd, agru, cedg, cedl, g2r, g2rd, g2ru, gggr, &
    gggrd, gggru, grgrd, grgru, gzgr, ro, vcgd, vcgu, vcld, vclu, vxgd, vxgu, &
    vxld, vxlu, xced, xcptd, xcptu, xedg, xedl, zta
  integer :: idspr
!..
!.. Local Scalars ..
  double precision :: a, a1, a2, a3, af, alc, alf, alfc, ap, b, b1, b1f, b1p, &
    b2, b2f, b2p, b3, bcr, beta, bf, bp, brs, bx, bxd, bxu, bz41, c, c1, c113, &
    c115, c13, c1415, c2, c23, c2915, c2q23, c3, c32, c43, c53, c56, c76, c83, &
    ca, ccf, ccp, ce, cef, cep, cf, cgz, cp, crdc, crf, cro, crp, crr1, crr2, &
    d, dacdr, dbdr, dbrod, dbrou, dcdr, dd, decdrf, decdrp, df, dfdz, dlta, &
    dp, dsdfd, dsdfu, dspr, dsprs, dvdr1, dvdr2, dvdrd, dvdru, ec, ecf, ecp, &
    ecrs, eczta, ef3vi, expfai, f1d, f1u, f2d, f2u, f3d, f3u, fai, fai2, fd, &
    fdd0, fk, fu, fz, g, gf, gp, gr2, gr2d, gr2u, gz, gz2, gz3, hugef, huges, &
    pi, q, q1, q2, q3, r, rnc, ro113, ro13, ro2, ro43, ro76, ro83, rod, rod13, &
    rod23, rod3, rod43, rod53, rou, rou13, rou23, rou3, rou43, rou53, rs, rs2, &
    rs3, sd, sd2, sd3, sd4, sd6, sidfd, sidfu, sk, sml, ssfc, su, su2, su3, &
    su4, su6, tc, td, tksg, tu, uc, ud, uu, vc, vc13, vc45d, vc45u, vc6, vccf, &
    vcf, vcl1, vcl2, vcp, vxp, vz, wc, x, x0, x01, x02, x03, xedgd, xedgu, &
    xedld, xedlu, xf, xl, xl0, xl01, xl02, xl03, xl1, xl2, xl3, xld, xld1, &
    xld2, xld3, xlf, xp, xs, zt13m, zt13p, zta3, zta4
  integer :: ibh, ica, icg, iex, igd, igh, igl, imj, ip9, ipg, ivg, ivn, ixlf
!..
!.. External Subroutines ..
  external :: corlsd, cpw91, exch91
!..
!.. Intrinsic Functions ..
  intrinsic :: acos, atan, exp, log, sqrt
!..
!.. Statement Functions ..
  double precision :: fbet, fdedr, fdfdz, ffz, fncecl, fncecs, fncf, fncvcl, &
    fncvcs, fvnec, fvq
!..
!.. Save statement ..
  save :: gp, gf, b1p, b1f, b2p, b2f, cp, cf, dp, df, ap, bp, af, bf, a1, x01, &
    b1, c1, a2, x02, b2, c2, a3, x03, b3, c3, fdd0, huges, hugef, dspr, igl, &
    igh, imj, ibh, ica, icg, ivn, ipg, ivg, ip9, igd, ixlf, iex, xlf
!..
!.. Statement Function definitions ..
  fncf(x) = (1.d0+x*x*x)*log(1.d0+1.d0/x) + x/2.d0 - x*x - 0.333333333d0
  fncecl(r, g, b1, b2) = g/(1.d0+b1*sqrt(r)+b2*r)
  fncvcl(ce, r, b1, b2) = ce*(1.d0+1.16666667d0*b1*sqrt(r)+1.33333333d0*b2*r)/ &
    (1.d0+b1*sqrt(r)+b2*r)
  fncecs(r, a, b, c, d) = a*log(r) + b + c*r*log(r) + d*r
  fncvcs(r, a, b, c, d) = a*log(r) + (b-a/3.d0) + 0.666666667d0*c*r*log(r) + &
    (2.d0*d-c)*r/3.d0
  ffz(zta) = 1.923661051d0*((1.d0+zta)**1.3333333333d0+(1.d0-zta)** &
    1.3333333333d0-2.d0)
  fdfdz(zta) = 2.564881401d0*((1.d0+zta)**.333333333333d0-(1.d0-zta)** &
    .333333333333d0)
  fvq(b, c) = sqrt(4.d0*c-b**2)
  fvnec(a, x, xl, x0, xl0, b, q) = a*(log(x*x/xl)+2.d0*b/q*atan(q/(2.d0*x+ &
    b))-b*x0/xl0*(log((x-x0)**2/xl)+2.d0*(b+2.d0*x0)/q*atan(q/(2.d0*x+b))))
  fbet(fdd0, ecf, ecp, alc) = fdd0*(ecf-ecp)/alc - 1.d0
  fdedr(ro, x, a, x0, xl, xl0, xld, b, q) = -x/(6.d0*ro)*a* &
    ((2.d0*xl-x*xld)/(x*xl)-b*(4.d0/(xld**2+q**2)+x0/xl0*((2.d0*xl-(x-x0)*xld) &
    /((x-x0)*xl)-4.d0*(b+2.d0*x0)/(xld**2+q**2))))
!..
!.. Data statements ..
  data gp, gf, b1p, b1f, b2p, b2f, cp, cf, dp, df/ -.2846d0, -.1686d0, &
    1.0529d0, 1.3981d0, 0.3334d0, 0.2611d0, 0.0040d0, 0.0014d0, -.0232d0, &
    -.0096d0/
  data ap, bp, af, bf/0.0622d0, -.096d0, 0.0311d0, -0.0538d0/
  data a1, x01, b1, c1/.0621814d0, -.10498d0, 3.72744d0, 12.9352d0/
  data a2, x02, b2, c2/.0310907d0, -.32500d0, 7.06042d0, 18.0578d0/
  data a3, x03, b3, c3/ -.03377373d0, -.0047584d0, 1.13107d0, 13.0045d0/
  data fdd0/1.70992093d0/
  data huges, hugef, dspr/1.d+6, 50.d0, 1.d-4/
  data igl, igh, imj, ibh, ica, icg, ivn, ipg, ivg, ip9, igd, ixlf, iex, &
    xlf/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00d0/
!..
!-----------------------------------------------------------------
!Perdew-zunger parametrization of Ceperley-Alder. g,a,b,c,d in ry.
!for vwn. a1,x01,b1,c1 for para(zta=0), -2 for zta=1, -3 for alfac.
!fdd0: fz''(o).
!-----------------------------------------------------------------

!     PW91   ip9=1  igd=1

!     PW86      imj=1  igd=1


  ip9 = 1
  igd = 1


  pi = acos(-1.d0)
  sml = 1.d-12
  if (zta>1.d0-sml) zta = 1.d0 - sml
!.....
!     vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),(up,dw)
!     vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),(up,dw)
!     xedl,xedg: exchange energy density (local,grad.exp.) in ry.
!     cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
!.....
  vxlu = 0.0d0
  vclu = 0.0d0
  vxld = 0.0d0
  vcld = 0.0d0
  xedl = 0.0d0
  cedl = 0.0d0
  vxgu = 0.0d0
  vcgu = 0.0d0
  vxgd = 0.0d0
  vcgd = 0.0d0
  xedg = 0.0d0
  cedg = 0.0d0
!.....
  if (ro<sml) go to 110
!.....
  c13 = 1.d0/3.d0
  c23 = 2.d0/3.d0
  c32 = 3.d0/2.d0
  c43 = 4.d0/3.d0
  c53 = 5.d0/3.d0
  c76 = 7.d0/6.d0
  c113 = 11.d0/3.d0
!.....ca=2.**(-.33333333)
  ca = 0.793700526d0
!.....alf=-3*(3/4*pai)**(1/3).
  alf = -1.861051473d0
!.....
  ro2 = ro*ro
  ro13 = ro**c13
  ro43 = ro**c43
  ro83 = ro43**2
  ro76 = ro**c76
  ro113 = ro**c113
!.....
  rou = ro*(1.d0+zta)/2.d0
  rou3 = rou**3
  rou13 = rou**c13
  rou23 = rou**c23
  rou43 = rou**c43
!.....
  rod = ro - rou
  rod3 = rod**3
  rod13 = rod**c13
  rod23 = rod**c23
  rod43 = rod**c43
!.....
!     gr2=drr*drr
!     gr2u=drru**2
!     drrd=drr-drru
!     gr2d=drrd**2
!     ddrrd=ddrr-ddrru
!.....
  fz = ffz(zta)
  dfdz = fdfdz(zta)
!.....
!.....gz,gz2,gz3: for Wang-Perdew ssf.
  gz = ((1.d0+zta)**c23+(1.d0-zta)**c23)/2.d0
  gz2 = gz**2
  gz3 = gz**3
!.....
  zta3 = zta**3
  zta4 = zta**4
  zt13p = (1.d0+zta)**c13
  zt13m = (1.d0-zta)**c13
!.....
  rs = 0.620350491d0/ro13
  rs2 = rs*rs
  rs3 = rs*rs2
!.....
!.....xedl: exchange-energy-density in ry.
  xedl = alf*(rou43+rod43)
!.....
!.....exchange-potential, vxp,vxlu,vxld: v-exchange-(para,up,dw).
  vxlu = c43*alf*rou13
  vxld = c43*alf*rod13

!.....
  if (iex==1) go to 100
!.....

!.....xlfa.
  if (ixlf/=0) then
    xlf = 2.d0/3.d0
    vclu = (xlf*c32-1.d0)*vxlu
    vcld = (xlf*c32-1.d0)*vxld
    cedl = (xlf*c32-1.d0)*xedl
    go to 100
  end if
!.....
!.....Gunnarson-Lundqvist.(p.r.b13('76),4274,eqs(54)~(56).) or
!c....  g-l but with beta by hedin-lundq.(j.phys.c.4('71),2064))
  if ((igl/=0) .or. (igh/=0)) then
    xp = rs/11.4d0
    xf = rs/15.9d0
    cep = -0.0666d0*fncf(xp)
    cef = -0.0406d0*fncf(xf)
    ce = cep + (cef-cep)*fz
    cedl = ce*ro
!.....
    if (igl/=0) beta = 1.d0 + 0.0545d0*rs*log(1.d0+11.4d0/rs)
    if (igh/=0) beta = 1.d0 + .03683d0*rs*log(1.d0+21.d0/rs)
    dlta = 1.d0 - 0.036d0*rs + 1.36d0*rs/(1.d0+10.d0*rs)
!.....
    vxp = c43*alf*ca*ro13
    vclu = vxp*(beta+dlta/3.d0*zta/(1.d0+0.297d0*zta)) - vxlu
    vcld = vxp*(beta-dlta/3.d0*zta/(1.d0-0.297d0*zta)) - vxld
!.....
    go to 100
  end if
!.....
!.....Hedin-von Barth. (j.phys.c.5('72),1629) or moruzzi-janak-williams.
  if ((ibh/=0) .or. (imj/=0)) then
    if (ibh/=0) then
      crp = 30.d0
      crf = 75.d0
      ccp = 0.0504d0
      ccf = 0.0254d0
    else if (imj/=0) then
      crp = 21.d0
      crf = 52.916684d0
      ccp = 0.045d0
      ccf = 0.0225d0
!           write(6,*) 'MJW'
    end if
    xp = rs/crp
    xf = rs/crf
    cep = -ccp*fncf(xp)
    cef = -ccf*fncf(xf)
    ce = cep + (cef-cep)*fz
    cedl = ce*ro
!       vclu,vcld: v-correlation-(up,dw). potential.(ry)
    rnc = c43*ca/(1.d0-ca)*(cef-cep)
    vcp = -ccp*log(1.d0+crp/rs)
    brs = vcp - rnc
    vclu = rnc*zt13p + brs
    vcld = rnc*zt13m + brs
!.....
    go to 100
  end if
!.....
!.....Ceperley-Alder.(paramtrzd by Perdew-zunger.(p.r.23('81),5048)).
  if (ica/=0) then
!.....
    if (rs>=1.d0) then
      cep = fncecl(rs, gp, b1p, b2p)
      cef = fncecl(rs, gf, b1f, b2f)
      vcp = fncvcl(cep, rs, b1p, b2p)
      vcf = fncvcl(cef, rs, b1f, b2f)
    else
      cep = fncecs(rs, ap, bp, cp, dp)
      cef = fncecs(rs, af, bf, cf, df)
      vcp = fncvcs(rs, ap, bp, cp, dp)
      vcf = fncvcs(rs, af, bf, cf, df)
    end if
!.....
    ce = cep + (cef-cep)*fz
    cedl = ce*ro
!.....
!.....
    vcl2 = (cef-cep)*dfdz
    vcl1 = vcp + (vcf-vcp)*fz - vcl2*zta
    vclu = vcl1 + vcl2
    vcld = vcl1 - vcl2
!.....
    go to 100
  end if
!.....
!.....Ceperley-Alder.with Wang-Perdew spin-scaling-factor.
  if (icg/=0) then
!.....
    if (rs>=1.d0) then
      cep = fncecl(rs, gp, b1p, b2p)
      vcp = fncvcl(cep, rs, b1p, b2p)
    else
      cep = fncecs(rs, ap, bp, cp, dp)
      vcp = fncvcs(rs, ap, bp, cp, dp)
    end if
!.....
    ce = cep*gz3
    cedl = ce*ro
!.....
    cgz = cep*gz2*(1.d0/zt13p-1.d0/zt13m)
    vcl1 = vcp*gz3 - cgz*zta
    vclu = vcp*gz3 + cgz
    vcld = vcp*gz3 - cgz
!.....
    go to 100
  end if
!.....
!.....Vosko-Wilk-Nusair. Phys.Rev..22,3812,'80.
  if (ivn/=0) then
!.....
!.....xl:x-large. xld:d(xl)/dx. xl0:x-large for x=x0.
    xs = sqrt(rs)
    xl1 = xs**2 + b1*xs + c1
    xl2 = xs**2 + b2*xs + c2
    xl3 = xs**2 + b3*xs + c3
    xld1 = 2.d0*xs + b1
    xld2 = 2.d0*xs + b2
    xld3 = 2.d0*xs + b3
    xl01 = x01**2 + b1*x01 + c1
    xl02 = x02**2 + b2*x02 + c2
    xl03 = x03**2 + b3*x03 + c3
    q1 = fvq(b1, c1)
    q2 = fvq(b2, c2)
    q3 = fvq(b3, c3)
    ecp = fvnec(a1, xs, xl1, x01, xl01, b1, q1)
    ecf = fvnec(a2, xs, xl2, x02, xl02, b2, q2)
    alc = fvnec(a3, xs, xl3, x03, xl03, b3, q3)
    beta = fbet(fdd0, ecf, ecp, alc)
    bz41 = 1.d0 + beta*zta4
!.....
    ce = ecp + alc*fz/fdd0*bz41
    cedl = ce*ro
!.....
!.....alc: alfac.
!.....decdrp,decdrf: d(ec)/dro-para(zta=0), -(zta=1).
!.....dacdr: d(alc)/dro.
!.....dbdr: d(beta)/dro.
    decdrp = fdedr(ro, xs, a1, x01, xl1, xl01, xld1, b1, q1)
    decdrf = fdedr(ro, xs, a2, x02, xl2, xl02, xld2, b2, q2)
    dacdr = fdedr(ro, xs, a3, x03, xl3, xl03, xld3, b3, q3)
!.....
    dbdr = fdd0*((decdrf-decdrp)*alc-(ecf-ecp)*dacdr)/alc**2
    vcl1 = ce + ro*(decdrp+(dacdr*fz*bz41+alc*fz*dbdr*zta4)/fdd0)
    vcl2 = 2.d0*alc/(fdd0*ro)*(dfdz*bz41+4.d0*fz*beta*zta3)
    vclu = vcl1 + vcl2*rod
    vcld = vcl1 + vcl2*(-rou)
!.....
    go to 100
  end if
!.....
  if (ip9==1) then
    go to 100
  end if
!.....
100 continue

!.....gradient expansion.
!.....
  if (igd<=0) go to 110
!       write(6,*)  '  GGA '
!.....
  gr2 = agr**2
  gr2u = agru**2
  gr2d = agrd**2

  c56 = 5.d0/6.d0
  c115 = 1.d0/15.d0
  c1415 = 14.d0/15.d0
  c2915 = 29.d0/15.d0
  c2q23 = 2.d0**c23
  c83 = 8.d0/3.d0
!.....
!.....  dsprs: divergence-suppress-factor.
!       if((log(dspr)+2.*log(agr)-c83*log(ro)).gt.8.0) go to 200
  dsprs = 1.d0
  if (idspr==1) dsprs = exp(-dspr*gr2/ro**c83)
!.....
!     agr,agru,agrd: abs(grad(rho)), for all, up, and down.
!c    gr2,gr2u,gr2d: grad(rho_all)**2, grad(rho_up)**2, grad(rho_d)**2.
!     g2r,g2ru,g2rd: laplacian rho_all, _up and _down.
!     gggru,-d: grad(rho)*grad(abs(grad(rho))) for all,up and down.
!     grgru,-d: grad(rho_all)*grad(rhor_up) and for down.

!       g2r=ddrr+2.*drr/rv
!.....
  rou53 = rou**c53
!.....
!.....  edrru: d(abs(d(rou)/dr))/dr, edrrd for down.
!       edrru=ddrru
!       if(drru.lt.0.) edrru=-ddrru
!.....
!       agr,agbru,-d: abs(grad(rho)),for rou, rod.
!       gggru,-d: grad(rho)*grad(abs(grad(rho))) for up and down.
!.....  su:at ro=2*rou. 1/(2(3*pai**2)**(1/3))*|grad(rou)|/rou**(4/3).
  su = 0.128278244d0*agru/rou43
  if (su>huges) go to 110
!       g2ru=ddrru+2.*drru/rv
  tu = .016455307d0*g2ru/rou53
  uu = 0.002110857d0*gggru/rou3

  if (ip9/=1) then

    su2 = su*su
    su3 = su*su2
    su4 = su2*su2
    su6 = su2*su4
!.....
    f1u = 1.d0 + 1.296d0*su2 + 14.d0*su4 + .2d0*su6
    f2u = 2.592d0 + 56.d0*su2 + 1.2d0*su4
    f3u = 112.d0*su + 4.8d0*su3
!.....
!.....  fu: fgga(su) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
!.....  sidfu: su**(-1)*d(fu)/d(su)).
!.....  dsdfu: d(sidfu)/d(su).
!.....  xedgu; exchange energy density xe at ro=2*rou.(16) of p.w.
!c....      xedgu=ax*rou**(4/3)*(fu-1). ax=2**(4/3)*1.47711(ry).

    fu = f1u**c115
    sidfu = c115*f1u**(-c1415)*f2u
    dsdfu = c115*f1u**(-c2915)*(-c1415*su*f2u**2+f1u*f3u)
!.....
    xedgu = -3.722102942d0*(fu-1.d0)*rou43
!.....
    vxgu = dsprs*alf*rou13*(c43*(fu-1.d0)-tu*sidfu-(uu-c43*su3)*dsdfu)

  else

    dbrou = rou*2.d0

    call exch91(dbrou, su, uu, tu, xedlu, xedgu, vxlu, vxgu)

    xedl = xedlu/2.d0

  end if

!.....
!.....bxu,bxd,bx: grad-coeff. for exchange.

  bxu = xedgu/gr2u*rou43
!.....
  rod53 = rod**c53
!       edrrd=ddrrd
!       if(drrd.lt.0.) edrrd=-ddrrd

  sd = 0.128278244d0*agrd/rod43
  if (sd>huges) go to 110

!       g2rd=ddrrd+2.*drrd/rv

  td = .016455307d0*g2rd/rod53
  ud = 0.002110857d0*gggrd/rod3

  if (ip9/=1) then

    sd2 = sd*sd
    sd3 = sd*sd2
    sd4 = sd2*sd2
    sd6 = sd2*sd4
!.....
    f1d = 1.d0 + 1.296d0*sd2 + 14.d0*sd4 + .2d0*sd6
    f2d = 2.592d0 + 56.d0*sd2 + 1.2d0*sd4
    f3d = 112.d0*sd + 4.8d0*sd3
!.....
!.....  fd: fgga(sd) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
!.....  sidfd: sd**(-1)*d(fd)/d(sd)).
!.....  dsdfd: d(sidfd)/d(sd).
!.....  xedgd; exchange energy density xe at ro=2*rod.(16) of p.w.
!c....      xedgd=ax*rod**(4/3)*(fd-1). ax=2**(4/3)*1.47711(ry).

    fd = f1d**c115
    sidfd = c115*f1d**(-c1415)*f2d
    dsdfd = c115*f1d**(-c2915)*(-c1415*sd*f2d**2+f1d*f3d)
!.....
    xedgd = -3.722102942d0*(fd-1.d0)*rod43
!.....
    vxgd = dsprs*alf*rod13*(c43*(fd-1.d0)-td*sidfd-(ud-c43*sd3)*dsdfd)

  else

    dbrod = rod*2.d0

    call exch91(dbrod, sd, ud, td, xedld, xedgd, vxld, vxgd)

    xedl = xedl + xedld/2.d0

  end if

  bxd = xedgd/gr2d*rod43
!.....

  xedg = dsprs*(xedgu+xedgd)/2.d0

  bx = (bxu+bxd)/2.d0

  if (iex==1) go to 110

!.....
!.... cro: c(n) of (6),Phys.Rev..b33,8822('86). in ry.
!.... dcdr: d(cro)/d(ro).
!.....0.001625816=1.745*f(=0.11)*cro(rs=0).

  if (ip9/=1) then

    crr1 = .005136d0 + .046532d0*rs + 1.4778d-5*rs2
    crr2 = 1.d0 + 8.723d0*rs + .472d0*rs2 + .07389d0*rs3
    cro = .003334d0 + crr1/crr2
    dcdr = ((.046532d0+2.9556d-5*rs)*crr2-crr1*(8.723d0+.944d0*rs+.22167d0*rs2 &
      ))/crr2/crr2*(-rs/ro/3.d0)
!.....
    fai = 0.001625816d0/cro*agr/ro76
    if (fai>hugef) go to 110
    fai2 = fai*fai
    expfai = exp(-fai)
!.....
!.....
    if (ipg==0) then

      dd = 0.707106781d0*sqrt((1.d0+zta)**c53+(1.d0-zta)**c53)
!.....    ssfc: spin-scaling-factor for gradient correlation energy.
      ssfc = 1.d0/dd
      crdc = c56/(ro113*dd**2)*c2q23
      vc45u = -crdc*(rou23-rod23)*((1.d0-fai)*rod*gr2-(2.d0-fai)*ro*grgrd)
      vc45d = -crdc*(rod23-rou23)*((1.d0-fai)*rou*gr2-(2.d0-fai)*ro*grgru)

    else if (ipg==1) then

      ssfc = gz
      crdc = c2q23/(3.d0*gz*ro83)
      vc45u = crdc*(1.d0/rou13-1.d0/rod13)*((1.d0-fai)*rod*gr2-(2.d0-fai)*ro* &
        grgrd)
      vc45d = crdc*(1.d0/rod13-1.d0/rou13)*((1.d0-fai)*rou*gr2-(2.d0-fai)*ro* &
        grgru)

    else if (ivg==1) then

      write (6, fmt='(/'' NON-SPHER MODIFICATION NOT COMPLETED FOR VG'')')
      stop 30

      if (ivn==0) then
        write (6, fmt=120) ivn, ivg
        stop 16
      end if
!.....
      dfdz = fdfdz(zta)
      vz = (1.d0+alc/ecp*fz/fdd0*bz41)**c13
!.....
      ssfc = vz
!.....
!.....    dvdru,dvdrd: d(vz)/drou,-d.
      ef3vi = 1.d0/(ecp*fdd0*3.d0*vz**2)
      dvdr1 = (dacdr*bz41-alc/ecp*decdrp*bz41+alc*dbdr*zta4)*fz*ef3vi
      dvdr2 = 2.d0*(dfdz*bz41+4.d0*fz*beta*zta3)*alc/ro2*ef3vi
      dvdru = dvdr1 + dvdr2*rod
      dvdrd = dvdr1 - dvdr2*rou
!.....
      vc45u = ((1.d0-fai)*gr2*dvdru-(2.d0-fai)*grgrd*(dvdru-dvdrd))/(vz*ro)
      vc45d = ((1.d0-fai)*gr2*dvdrd-(2.d0-fai)*grgru*(dvdrd-dvdru))/(vz*ro)

    end if

!.....  cedg: correlation-energy-density from grad.expansion.
!.....  bcr: grad-coeff. for correlation.
    bcr = ssfc*expfai*cro
    cedg = dsprs*bcr*gr2/ro43
!.....
!.....  vccf:v-correlation-coeff.
    vccf = -ssfc*expfai*cro/ro13
    vc13 = (2.d0-fai)*g2r/ro - (c43-c113*fai+c76*fai2)*gr2/ro2 + &
      fai*(fai-3.d0)*gggr/agr/ro
!    &    fai*(fai-3.)*ddrr/ro
    vc6 = -gr2/ro*(fai2-fai-1.d0)/cro*dcdr
!.....
    vcgu = dsprs*vccf*(vc13+vc6+vc45u)
!.....
    vcgd = dsprs*vccf*(vc13+vc6+vc45d)

  else

!       PW91

    call corlsd(rs, zta, ec, vclu, vcld, ecrs, eczta, alfc)

    vclu = vclu*2.d0
    vcld = vcld*2.d0
    cedl = ec*2.d0*ro

    fk = 1.91915829d0/rs
    sk = sqrt(4.d0*fk/pi)
    tksg = 2.d0*sk*gz
    tc = agr/(ro*tksg)
!           gagr: d(ABS(d(ro)/dr))/dr.
!           gagr=ddrr
!           if(drr.lt.0.) gagr=-ddrr
    uc = gggr/(ro2*tksg**3)
!         uc=drr*gagr/(ro2*tksg**3)
    vc = g2r/(ro*tksg**2)
    wc = gzgr/(ro*tksg**2)
!         wc=drr*dzr/(ro*tksg**2)

    call cpw91(fk, sk, gz, ec, ecrs, eczta, rs, zta, tc, uc, vc, wc, cedg, &
      vcgu, vcgd)

    vcgu = vcgu*2.d0
    vcgd = vcgd*2.d0
    cedg = cedg*ro*2.d0*dsprs

    bcr = cedg/gr2*ro43

  end if
!.....
110 continue

  xcptu = vxlu + vclu + vxgu + vcgu
  xcptd = vxld + vcld + vxgd + vcgd
!heck
!     ro is small

  xced = 0.0d0
  if (ro>sml) xced = (xedl+cedl+xedg+cedg)/ro

!     write(6,'(/'' vxlu,vxld,vclu,vcld,xedl,cedl ro='',7f11.5)') vxlu,
!    &    vxld,vclu,vcld,xedl,cedl,ro
!       write(6,'(/'' vxgu,vxgd,vcgu,vcgd,xedg,cedg='',6f12.7)') vxgu,
!    &    vxgd,vcgu,vcgd,xedg,cedg

  return
120 format (/, ' ivn should be 1 for ivg=1. ivn,ivg=', 2i5, /)
end subroutine
