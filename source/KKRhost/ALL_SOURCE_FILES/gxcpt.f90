module mod_gxcpt
  
  private
  public :: gxcpt

contains

  !-------------------------------------------------------------------------------
  !> Summary: Exchange-correlation potential and total energy for PW91 (GGA)
  !> Author: 
  !> Category: KKRhost, xc-potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Energies in Ry
  !-------------------------------------------------------------------------------
  subroutine gxcpt(idspr, ro, zta, agr, agru, agrd, g2r, g2ru, g2rd, gggr, gggru, gggrd, grgru, grgrd, gzgr, xcptu, xcptd, xced, vxlu, vxld, vclu, vcld, xedl, cedl, vxgu, vxgd, &
    vcgu, vcgd, xedg, cedg)

    use :: mod_datatypes, only: dp
    use :: mod_corlsd, only: corlsd
    use :: mod_cpw91, only: cpw91
    use :: mod_exch91, only: exch91
    use :: mod_constants, only: pi
    implicit none
    real (kind=dp), parameter :: sml = 1.e-12_dp
    ! .. Scalar Arguments ..
    real (kind=dp) :: agr, agrd, agru, cedg, cedl, g2r, g2rd, g2ru, gggr, gggrd, gggru, grgrd, grgru, gzgr, ro, vcgd, vcgu, vcld, vclu, vxgd, vxgu, vxld, vxlu, xced, xcptd, xcptu, &
      xedg, xedl, zta
    integer :: idspr
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: a1, a2, a3, af, alc, alf, alfc, ap, b1, b1f, b1p, b2, b2f, b2p, b3, bcr, beta, bf, bp, brs, bx, bxd, bxu, bz41, c1, c113, c115, c13, c1415, c2, c23, c2915, &
      c2q23, c3, c32, c43, c53, c56, c76, c83, ca, ccf, ccp, ce, cef, cep, cf, cgz, cp, crdc, crf, cro, crp, crr1, crr2, dacdr, dbdr, dbrod, dbrou, dcdr, dd, decdrf, decdrp, df, &
      dfdz, dlta, d_p, dsdfd, dsdfu, dspr, dsprs, dvdr1, dvdr2, dvdrd, dvdru, ec, ecf, ecp, ecrs, eczta, ef3vi, expfai, f1d, f1u, f2d, f2u, f3d, f3u, fai, fai2, fd, fdd0, fk, fu, &
      fz, gf, gp, gr2, gr2d, gr2u, gz, gz2, gz3, hugef, huges, q1, q2, q3, rnc, ro113, ro13, ro2, ro43, ro76, ro83, rod, rod13, rod23, rod3, rod43, rod53, rou, rou13, rou23, &
      rou3, rou43, rou53, rs, rs2, rs3, sd, sd2, sd3, sd4, sd6, sidfd, sidfu, sk, ssfc, su, su2, su3, su4, su6, tc, td, tksg, tu, uc, ud, uu, vc, vc13, vc45d, vc45u, vc6, &
      vccf, vcf, vcl1, vcl2, vcp, vxp, vz, wc, x01, x02, x03, xedgd, xedgu, xedld, xedlu, xf, xl01, xl02, xl03, xl1, xl2, xl3, xld1, xld2, xld3, xlf, xp, xs, zt13m, zt13p, zta3, &
      zta4
    integer :: ibh, ica, icg, iex, igd, igh, igl, imj, ip9, ipg, ivg, ivn, ixlf
    ! ..
    ! .. Save statement ..
    save :: gp, gf, b1p, b1f, b2p, b2f, cp, cf, d_p, df, ap, bp, af, bf, a1, x01, b1, c1, a2, x02, b2, c2, a3, x03, b3, c3, fdd0, huges, hugef, dspr, igl, igh, imj, ibh, ica, icg, &
      ivn, ipg, ivg, ip9, igd, ixlf, iex, xlf
    ! ..
    ! .. Data statements ..
    data gp, gf, b1p, b1f, b2p, b2f, cp, cf, d_p, df/ -.2846e0_dp, -.1686e0_dp, 1.0529e0_dp, 1.3981e0_dp, 0.3334e0_dp, 0.2611e0_dp, 0.0040e0_dp, 0.0014e0_dp, -.0232e0_dp, &
      -.0096e0_dp/
    data ap, bp, af, bf/0.0622e0_dp, -.096e0_dp, 0.0311e0_dp, -0.0538e0_dp/
    data a1, x01, b1, c1/.0621814e0_dp, -.10498e0_dp, 3.72744e0_dp, 12.9352e0_dp/
    data a2, x02, b2, c2/.0310907e0_dp, -.32500e0_dp, 7.06042e0_dp, 18.0578e0_dp/
    data a3, x03, b3, c3/ -.03377373e0_dp, -.0047584e0_dp, 1.13107e0_dp, 13.0045e0_dp/
    data fdd0/1.70992093e0_dp/
    data huges, hugef, dspr/1.e+6_dp, 50.e0_dp, 1.e-4_dp/
    data igl, igh, imj, ibh, ica, icg, ivn, ipg, ivg, ip9, igd, ixlf, iex, xlf/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00e0_dp/
    ! ..
    ! -----------------------------------------------------------------
    ! Perdew-zunger parametrization of Ceperley-Alder. g,a,b,c,d in ry.
    ! for vwn. a1,x01,b1,c1 for para(zta=0), -2 for zta=1, -3 for alfac.
    ! fdd0: fz''(o).
    ! -----------------------------------------------------------------

    ! PW91   ip9=1  igd=1
    ! PW86      imj=1  igd=1

    ip9 = 1
    igd = 1

    if (zta>1.e0_dp-sml) zta = 1.e0_dp - sml
    ! .....
    ! vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),(up,dw)
    ! vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),(up,dw)
    ! xedl,xedg: exchange energy density (local,grad.exp.) in ry.
    ! cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
    ! .....
    vxlu = 0.0e0_dp
    vclu = 0.0e0_dp
    vxld = 0.0e0_dp
    vcld = 0.0e0_dp
    xedl = 0.0e0_dp
    cedl = 0.0e0_dp
    vxgu = 0.0e0_dp
    vcgu = 0.0e0_dp
    vxgd = 0.0e0_dp
    vcgd = 0.0e0_dp
    xedg = 0.0e0_dp
    cedg = 0.0e0_dp
    ! .....
    if (ro<sml) go to 110
    ! .....
    c13 = 1.e0_dp/3.e0_dp
    c23 = 2.e0_dp/3.e0_dp
    c32 = 3.e0_dp/2.e0_dp
    c43 = 4.e0_dp/3.e0_dp
    c53 = 5.e0_dp/3.e0_dp
    c76 = 7.e0_dp/6.e0_dp
    c113 = 11.e0_dp/3.e0_dp
    ! .....ca=2.**(-.33333333)
    ca = 0.793700526e0_dp
    ! .....alf=-3*(3/4*pai)**(1/3).
    alf = -1.861051473e0_dp
    ! .....
    ro2 = ro*ro
    ro13 = ro**c13
    ro43 = ro**c43
    ro83 = ro43**2
    ro76 = ro**c76
    ro113 = ro**c113
    ! .....
    rou = ro*(1.e0_dp+zta)/2.e0_dp
    rou3 = rou**3
    rou13 = rou**c13
    rou23 = rou**c23
    rou43 = rou**c43
    ! .....
    rod = ro - rou
    rod3 = rod**3
    rod13 = rod**c13
    rod23 = rod**c23
    rod43 = rod**c43
    ! .....
    ! gr2=drr*drr
    ! gr2u=drru**2
    ! drrd=drr-drru
    ! gr2d=drrd**2
    ! ddrrd=ddrr-ddrru
    ! .....
    fz = ffz(zta)
    dfdz = fdfdz(zta)
    ! .....
    ! .....gz,gz2,gz3: for Wang-Perdew ssf.
    gz = ((1.e0_dp+zta)**c23+(1.e0_dp-zta)**c23)/2.e0_dp
    gz2 = gz**2
    gz3 = gz**3
    ! .....
    zta3 = zta**3
    zta4 = zta**4
    zt13p = (1.e0_dp+zta)**c13
    zt13m = (1.e0_dp-zta)**c13
    ! .....
    rs = 0.620350491e0_dp/ro13
    rs2 = rs*rs
    rs3 = rs*rs2
    ! .....
    ! .....xedl: exchange-energy-density in ry.
    xedl = alf*(rou43+rod43)
    ! .....
    ! .....exchange-potential, vxp,vxlu,vxld: v-exchange-(para,up,dw).
    vxlu = c43*alf*rou13
    vxld = c43*alf*rod13

    ! .....
    if (iex==1) go to 100
    ! .....

    ! .....xlfa.
    if (ixlf/=0) then
      xlf = 2.e0_dp/3.e0_dp
      vclu = (xlf*c32-1.e0_dp)*vxlu
      vcld = (xlf*c32-1.e0_dp)*vxld
      cedl = (xlf*c32-1.e0_dp)*xedl
      go to 100
    end if
    ! .....
    ! .....Gunnarson-Lundqvist.(p.r.b13('76),4274,eqs(54)~(56).) or
    ! c....  g-l but with beta by hedin-lundq.(j.phys.c.4('71),2064))
    if ((igl/=0) .or. (igh/=0)) then
      xp = rs/11.4e0_dp
      xf = rs/15.9e0_dp
      cep = -0.0666e0_dp*fncf(xp)
      cef = -0.0406e0_dp*fncf(xf)
      ce = cep + (cef-cep)*fz
      cedl = ce*ro
      ! .....
      if (igl/=0) beta = 1.e0_dp + 0.0545e0_dp*rs*log(1.e0_dp+11.4e0_dp/rs)
      if (igh/=0) beta = 1.e0_dp + .03683e0_dp*rs*log(1.e0_dp+21.e0_dp/rs)
      dlta = 1.e0_dp - 0.036e0_dp*rs + 1.36e0_dp*rs/(1.e0_dp+10.e0_dp*rs)
      ! .....
      vxp = c43*alf*ca*ro13
      vclu = vxp*(beta+dlta/3.e0_dp*zta/(1.e0_dp+0.297e0_dp*zta)) - vxlu
      vcld = vxp*(beta-dlta/3.e0_dp*zta/(1.e0_dp-0.297e0_dp*zta)) - vxld
      ! .....
      go to 100
    end if
    ! .....
    ! .....Hedin-von Barth. (j.phys.c.5('72),1629) or moruzzi-janak-williams.
    if ((ibh/=0) .or. (imj/=0)) then
      if (ibh/=0) then
        crp = 30.e0_dp
        crf = 75.e0_dp
        ccp = 0.0504e0_dp
        ccf = 0.0254e0_dp
      else if (imj/=0) then
        crp = 21.e0_dp
        crf = 52.916684e0_dp
        ccp = 0.045e0_dp
        ccf = 0.0225e0_dp
        ! write(6,*) 'MJW'
      end if
      xp = rs/crp
      xf = rs/crf
      cep = -ccp*fncf(xp)
      cef = -ccf*fncf(xf)
      ce = cep + (cef-cep)*fz
      cedl = ce*ro
      ! vclu,vcld: v-correlation-(up,dw). potential.(ry)
      rnc = c43*ca/(1.e0_dp-ca)*(cef-cep)
      vcp = -ccp*log(1.e0_dp+crp/rs)
      brs = vcp - rnc
      vclu = rnc*zt13p + brs
      vcld = rnc*zt13m + brs
      ! .....
      go to 100
    end if
    ! .....
    ! .....Ceperley-Alder.(paramtrzd by Perdew-zunger.(p.r.23('81),5048)).
    if (ica/=0) then
      ! .....
      if (rs>=1.e0_dp) then
        cep = fncecl(rs, gp, b1p, b2p)
        cef = fncecl(rs, gf, b1f, b2f)
        vcp = fncvcl(cep, rs, b1p, b2p)
        vcf = fncvcl(cef, rs, b1f, b2f)
      else
        cep = fncecs(rs, ap, bp, cp, d_p)
        cef = fncecs(rs, af, bf, cf, df)
        vcp = fncvcs(rs, ap, bp, cp, d_p)
        vcf = fncvcs(rs, af, bf, cf, df)
      end if
      ! .....
      ce = cep + (cef-cep)*fz
      cedl = ce*ro
      ! .....
      ! .....
      vcl2 = (cef-cep)*dfdz
      vcl1 = vcp + (vcf-vcp)*fz - vcl2*zta
      vclu = vcl1 + vcl2
      vcld = vcl1 - vcl2
      ! .....
      go to 100
    end if
    ! .....
    ! .....Ceperley-Alder.with Wang-Perdew spin-scaling-factor.
    if (icg/=0) then
      ! .....
      if (rs>=1.e0_dp) then
        cep = fncecl(rs, gp, b1p, b2p)
        vcp = fncvcl(cep, rs, b1p, b2p)
      else
        cep = fncecs(rs, ap, bp, cp, d_p)
        vcp = fncvcs(rs, ap, bp, cp, d_p)
      end if
      ! .....
      ce = cep*gz3
      cedl = ce*ro
      ! .....
      cgz = cep*gz2*(1.e0_dp/zt13p-1.e0_dp/zt13m)
      vcl1 = vcp*gz3 - cgz*zta
      vclu = vcp*gz3 + cgz
      vcld = vcp*gz3 - cgz
      ! .....
      go to 100
    end if
    ! .....
    ! .....Vosko-Wilk-Nusair. Phys.Rev..22,3812,'80.
    if (ivn/=0) then
      ! .....
      ! .....xl:x-large. xld:d(xl)/dx. xl0:x-large for x=x0.
      xs = sqrt(rs)
      xl1 = xs**2 + b1*xs + c1
      xl2 = xs**2 + b2*xs + c2
      xl3 = xs**2 + b3*xs + c3
      xld1 = 2.e0_dp*xs + b1
      xld2 = 2.e0_dp*xs + b2
      xld3 = 2.e0_dp*xs + b3
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
      bz41 = 1.e0_dp + beta*zta4
      ! .....
      ce = ecp + alc*fz/fdd0*bz41
      cedl = ce*ro
      ! .....
      ! .....alc: alfac.
      ! .....decdrp,decdrf: d(ec)/dro-para(zta=0), -(zta=1).
      ! .....dacdr: d(alc)/dro.
      ! .....dbdr: d(beta)/dro.
      decdrp = fdedr(ro, xs, a1, x01, xl1, xl01, xld1, b1, q1)
      decdrf = fdedr(ro, xs, a2, x02, xl2, xl02, xld2, b2, q2)
      dacdr = fdedr(ro, xs, a3, x03, xl3, xl03, xld3, b3, q3)
      ! .....
      dbdr = fdd0*((decdrf-decdrp)*alc-(ecf-ecp)*dacdr)/alc**2
      vcl1 = ce + ro*(decdrp+(dacdr*fz*bz41+alc*fz*dbdr*zta4)/fdd0)
      vcl2 = 2.e0_dp*alc/(fdd0*ro)*(dfdz*bz41+4.e0_dp*fz*beta*zta3)
      vclu = vcl1 + vcl2*rod
      vcld = vcl1 + vcl2*(-rou)
      ! .....
      go to 100
    end if
    ! .....
    if (ip9==1) then
      go to 100
    end if
    ! .....
100 continue

    ! .....gradient expansion.
    ! .....
    if (igd<=0) go to 110
    ! write(6,*)  '  GGA '
    ! .....
    gr2 = agr**2
    gr2u = agru**2
    gr2d = agrd**2

    c56 = 5.e0_dp/6.e0_dp
    c115 = 1.e0_dp/15.e0_dp
    c1415 = 14.e0_dp/15.e0_dp
    c2915 = 29.e0_dp/15.e0_dp
    c2q23 = 2.e0_dp**c23
    c83 = 8.e0_dp/3.e0_dp
    ! .....
    ! .....  dsprs: divergence-suppress-factor.
    ! if((log(dspr)+2.*log(agr)-c83*log(ro)).gt.8.0) go to 200
    dsprs = 1.e0_dp
    if (idspr==1) dsprs = exp(-dspr*gr2/ro**c83)
    ! .....
    ! agr,agru,agrd: abs(grad(rho)), for all, up, and down.
    ! c    gr2,gr2u,gr2d: grad(rho_all)**2, grad(rho_up)**2, grad(rho_d)**2.
    ! g2r,g2ru,g2rd: laplacian rho_all, _up and _down.
    ! gggru,-d: grad(rho)*grad(abs(grad(rho))) for all,up and down.
    ! grgru,-d: grad(rho_all)*grad(rhor_up) and for down.

    ! g2r=ddrr+2.*drr/rv
    ! .....
    rou53 = rou**c53
    ! .....
    ! .....  edrru: d(abs(d(rou)/dr))/dr, edrrd for down.
    ! edrru=ddrru
    ! if(drru.lt.0.) edrru=-ddrru
    ! .....
    ! agr,agbru,-d: abs(grad(rho)),for rou, rod.
    ! gggru,-d: grad(rho)*grad(abs(grad(rho))) for up and down.
    ! .....  su:at ro=2*rou. 1/(2(3*pai**2)**(1/3))*|grad(rou)|/rou**(4/3).
    su = 0.128278244e0_dp*agru/rou43
    if (su>huges) go to 110
    ! g2ru=ddrru+2.*drru/rv
    tu = .016455307e0_dp*g2ru/rou53
    uu = 0.002110857e0_dp*gggru/rou3

    if (ip9/=1) then

      su2 = su*su
      su3 = su*su2
      su4 = su2*su2
      su6 = su2*su4
      ! .....
      f1u = 1.e0_dp + 1.296e0_dp*su2 + 14.e0_dp*su4 + .2e0_dp*su6
      f2u = 2.592e0_dp + 56.e0_dp*su2 + 1.2e0_dp*su4
      f3u = 112.e0_dp*su + 4.8e0_dp*su3
      ! .....
      ! .....  fu: fgga(su) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
      ! .....  sidfu: su**(-1)*d(fu)/d(su)).
      ! .....  dsdfu: d(sidfu)/d(su).
      ! .....  xedgu; exchange energy density xe at ro=2*rou.(16) of p.w.
      ! c....      xedgu=ax*rou**(4/3)*(fu-1). ax=2**(4/3)*1.47711(ry).

      fu = f1u**c115
      sidfu = c115*f1u**(-c1415)*f2u
      dsdfu = c115*f1u**(-c2915)*(-c1415*su*f2u**2+f1u*f3u)
      ! .....
      xedgu = -3.722102942e0_dp*(fu-1.e0_dp)*rou43
      ! .....
      vxgu = dsprs*alf*rou13*(c43*(fu-1.e0_dp)-tu*sidfu-(uu-c43*su3)*dsdfu)

    else

      dbrou = rou*2.e0_dp

      call exch91(dbrou, su, uu, tu, xedlu, xedgu, vxlu, vxgu)

      xedl = xedlu/2.e0_dp

    end if

    ! .....
    ! .....bxu,bxd,bx: grad-coeff. for exchange.

    bxu = xedgu/gr2u*rou43
    ! .....
    rod53 = rod**c53
    ! edrrd=ddrrd
    ! if(drrd.lt.0.) edrrd=-ddrrd

    sd = 0.128278244e0_dp*agrd/rod43
    if (sd>huges) go to 110

    ! g2rd=ddrrd+2.*drrd/rv

    td = .016455307e0_dp*g2rd/rod53
    ud = 0.002110857e0_dp*gggrd/rod3

    if (ip9/=1) then

      sd2 = sd*sd
      sd3 = sd*sd2
      sd4 = sd2*sd2
      sd6 = sd2*sd4
      ! .....
      f1d = 1.e0_dp + 1.296e0_dp*sd2 + 14.e0_dp*sd4 + .2e0_dp*sd6
      f2d = 2.592e0_dp + 56.e0_dp*sd2 + 1.2e0_dp*sd4
      f3d = 112.e0_dp*sd + 4.8e0_dp*sd3
      ! .....
      ! .....  fd: fgga(sd) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
      ! .....  sidfd: sd**(-1)*d(fd)/d(sd)).
      ! .....  dsdfd: d(sidfd)/d(sd).
      ! .....  xedgd; exchange energy density xe at ro=2*rod.(16) of p.w.
      ! c....      xedgd=ax*rod**(4/3)*(fd-1). ax=2**(4/3)*1.47711(ry).

      fd = f1d**c115
      sidfd = c115*f1d**(-c1415)*f2d
      dsdfd = c115*f1d**(-c2915)*(-c1415*sd*f2d**2+f1d*f3d)
      ! .....
      xedgd = -3.722102942e0_dp*(fd-1.e0_dp)*rod43
      ! .....
      vxgd = dsprs*alf*rod13*(c43*(fd-1.e0_dp)-td*sidfd-(ud-c43*sd3)*dsdfd)

    else

      dbrod = rod*2.e0_dp

      call exch91(dbrod, sd, ud, td, xedld, xedgd, vxld, vxgd)

      xedl = xedl + xedld/2.e0_dp

    end if

    bxd = xedgd/gr2d*rod43
    ! .....

    xedg = dsprs*(xedgu+xedgd)/2.e0_dp

    bx = (bxu+bxd)/2.e0_dp

    if (iex==1) go to 110

    ! .....
    ! .... cro: c(n) of (6),Phys.Rev..b33,8822('86). in ry.
    ! .... dcdr: d(cro)/d(ro).
    ! .....0.001625816=1.745*f(=0.11)*cro(rs=0).

    if (ip9/=1) then

      crr1 = .005136e0_dp + .046532e0_dp*rs + 1.4778e-5_dp*rs2
      crr2 = 1.e0_dp + 8.723e0_dp*rs + .472e0_dp*rs2 + .07389e0_dp*rs3
      cro = .003334e0_dp + crr1/crr2
      dcdr = ((.046532e0_dp+2.9556e-5_dp*rs)*crr2-crr1*(8.723e0_dp+.944e0_dp*rs+.22167e0_dp*rs2))/crr2/crr2*(-rs/ro/3.e0_dp)
      ! .....
      fai = 0.001625816e0_dp/cro*agr/ro76
      if (fai>hugef) go to 110
      fai2 = fai*fai
      expfai = exp(-fai)
      ! .....
      ! .....
      if (ipg==0) then

        dd = 0.707106781e0_dp*sqrt((1.e0_dp+zta)**c53+(1.e0_dp-zta)**c53)
        ! .....    ssfc: spin-scaling-factor for gradient correlation energy.
        ssfc = 1.e0_dp/dd
        crdc = c56/(ro113*dd**2)*c2q23
        vc45u = -crdc*(rou23-rod23)*((1.e0_dp-fai)*rod*gr2-(2.e0_dp-fai)*ro*grgrd)
        vc45d = -crdc*(rod23-rou23)*((1.e0_dp-fai)*rou*gr2-(2.e0_dp-fai)*ro*grgru)

      else if (ipg==1) then

        ssfc = gz
        crdc = c2q23/(3.e0_dp*gz*ro83)
        vc45u = crdc*(1.e0_dp/rou13-1.e0_dp/rod13)*((1.e0_dp-fai)*rod*gr2-(2.e0_dp-fai)*ro*grgrd)
        vc45d = crdc*(1.e0_dp/rod13-1.e0_dp/rou13)*((1.e0_dp-fai)*rou*gr2-(2.e0_dp-fai)*ro*grgru)

      else if (ivg==1) then

        write (6, fmt='(/'' NON-SPHER MODIFICATION NOT COMPLETED FOR VG'')')
        stop 30

        if (ivn==0) then
          write (6, fmt=120) ivn, ivg
          stop 16
        end if
        ! .....
        dfdz = fdfdz(zta)
        vz = (1.e0_dp+alc/ecp*fz/fdd0*bz41)**c13
        ! .....
        ssfc = vz
        ! .....
        ! .....    dvdru,dvdrd: d(vz)/drou,-d.
        ef3vi = 1.e0_dp/(ecp*fdd0*3.e0_dp*vz**2)
        dvdr1 = (dacdr*bz41-alc/ecp*decdrp*bz41+alc*dbdr*zta4)*fz*ef3vi
        dvdr2 = 2.e0_dp*(dfdz*bz41+4.e0_dp*fz*beta*zta3)*alc/ro2*ef3vi
        dvdru = dvdr1 + dvdr2*rod
        dvdrd = dvdr1 - dvdr2*rou
        ! .....
        vc45u = ((1.e0_dp-fai)*gr2*dvdru-(2.e0_dp-fai)*grgrd*(dvdru-dvdrd))/(vz*ro)
        vc45d = ((1.e0_dp-fai)*gr2*dvdrd-(2.e0_dp-fai)*grgru*(dvdrd-dvdru))/(vz*ro)

      end if

      ! .....  cedg: correlation-energy-density from grad.expansion.
      ! .....  bcr: grad-coeff. for correlation.
      bcr = ssfc*expfai*cro
      cedg = dsprs*bcr*gr2/ro43
      ! .....
      ! .....  vccf:v-correlation-coeff.
      vccf = -ssfc*expfai*cro/ro13
      vc13 = (2.e0_dp-fai)*g2r/ro - (c43-c113*fai+c76*fai2)*gr2/ro2 + fai*(fai-3.e0_dp)*gggr/agr/ro
      ! &    fai*(fai-3.)*ddrr/ro
      vc6 = -gr2/ro*(fai2-fai-1.e0_dp)/cro*dcdr
      ! .....
      vcgu = dsprs*vccf*(vc13+vc6+vc45u)
      ! .....
      vcgd = dsprs*vccf*(vc13+vc6+vc45d)

    else

      ! PW91

      call corlsd(rs, zta, ec, vclu, vcld, ecrs, eczta, alfc)

      vclu = vclu*2.e0_dp
      vcld = vcld*2.e0_dp
      cedl = ec*2.e0_dp*ro

      fk = 1.91915829e0_dp/rs
      sk = sqrt(4.e0_dp*fk/pi)
      tksg = 2.e0_dp*sk*gz
      tc = agr/(ro*tksg)
      ! gagr: d(ABS(d(ro)/dr))/dr.
      ! gagr=ddrr
      ! if(drr.lt.0.) gagr=-ddrr
      uc = gggr/(ro2*tksg**3)
      ! uc=drr*gagr/(ro2*tksg**3)
      vc = g2r/(ro*tksg**2)
      wc = gzgr/(ro*tksg**2)
      ! wc=drr*dzr/(ro*tksg**2)

      call cpw91(fk, sk, gz, ec, ecrs, eczta, rs, zta, tc, uc, vc, wc, cedg, vcgu, vcgd)

      vcgu = vcgu*2.e0_dp
      vcgd = vcgd*2.e0_dp
      cedg = cedg*ro*2.e0_dp*dsprs

      bcr = cedg/gr2*ro43

    end if
    ! .....
110 continue

    xcptu = vxlu + vclu + vxgu + vcgu
    xcptd = vxld + vcld + vxgd + vcgd
    ! heck
    ! ro is small

    xced = 0.0e0_dp
    if (ro>sml) xced = (xedl+cedl+xedg+cedg)/ro

    ! write(6,'(/'' vxlu,vxld,vclu,vcld,xedl,cedl ro='',7f11.5)') vxlu,
    ! &    vxld,vclu,vcld,xedl,cedl,ro
    ! write(6,'(/'' vxgu,vxgd,vcgu,vcgd,xedg,cedg='',6f12.7)') vxgu,
    ! &    vxgd,vcgu,vcgd,xedg,cedg

    return
120 format (/, ' ivn should be 1 for ivg=1. ivn,ivg=', 2i5, /)
  end subroutine gxcpt

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fncf(x)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fncf
    real (kind=dp), intent (in) :: x

    fncf = (1.e0_dp+x*x*x)*log(1.e0_dp+1.e0_dp/x) + x/2.e0_dp - x*x - 0.333333333e0_dp
  end function fncf

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fncecl(r, g, b1, b2)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fncecl
    real (kind=dp), intent (in) :: r, g, b1, b2

    fncecl = g/(1.e0_dp+b1*sqrt(r)+b2*r)
  end function fncecl

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fncvcl(ce, r, b1, b2)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fncvcl
    real (kind=dp), intent (in) :: ce, r, b1, b2

    fncvcl = ce*(1.e0_dp+1.16666667e0_dp*b1*sqrt(r)+1.33333333e0_dp*b2*r)/(1.e0_dp+b1*sqrt(r)+b2*r)
  end function fncvcl

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fncecs(r, a, b, c, d)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fncecs
    real (kind=dp), intent (in) :: r, a, b, c, d

    fncecs = a*log(r) + b + c*r*log(r) + d*r
  end function fncecs

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fncvcs(r, a, b, c, d)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fncvcs
    real (kind=dp), intent (in) :: r, a, b, c, d

    fncvcs = a*log(r) + (b-a/3.e0_dp) + 0.666666667e0_dp*c*r*log(r) + (2.e0_dp*d-c)*r/3.e0_dp
  end function fncvcs

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function ffz(zta)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: ffz
    real (kind=dp), intent (in) :: zta

    ffz = 1.923661051e0_dp*((1.e0_dp+zta)**1.3333333333e0_dp+(1.e0_dp-zta)**1.3333333333e0_dp-2.e0_dp)
  end function ffz

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fdfdz(zta)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fdfdz
    real (kind=dp), intent (in) :: zta

    fdfdz = 2.564881401e0_dp*((1.e0_dp+zta)**.333333333333e0_dp-(1.e0_dp-zta)**.333333333333e0_dp)
  end function fdfdz

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fvq(b, c)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fvq
    real (kind=dp), intent (in) :: b, c

    fvq = sqrt(4.e0_dp*c-b**2)
  end function fvq

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fvnec(a, x, xl, x0, xl0, b, q)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fvnec
    real (kind=dp), intent (in) :: a, x, xl, x0, xl0, b, q

    fvnec = a*(log(x*x/xl)+2.e0_dp*b/q*atan(q/(2.e0_dp*x+b))-b*x0/xl0*(log((x-x0)**2/xl)+2.e0_dp*(b+2.e0_dp*x0)/q*atan(q/(2.e0_dp*x+b))))
  end function fvnec

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fbet(fdd0, ecf, ecp, alc)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fbet
    real (kind=dp), intent (in) :: fdd0, ecf, ecp, alc

    fbet = fdd0*(ecf-ecp)/alc - 1.e0_dp
  end function fbet

!-------------------------------------------------------------------------------
!> Summary: 
!> Author:
!> Date: 
!> Category: KKRhost, numerical-tools
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
  function fdedr(ro, x, a, x0, xl, xl0, xld, b, q)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp) :: fdedr
    real (kind=dp), intent (in) :: ro, x, a, x0, xl, xl0, xld, b, q

    fdedr = -x/(6.e0_dp*ro)*a*((2.e0_dp*xl-x*xld)/(x*xl)-b*(4.e0_dp/(xld**2+q**2)+x0/xl0*((2.e0_dp*xl-(x-x0)*xld)/((x-x0)*xl)-4.e0_dp*(b+2.e0_dp*x0)/(xld**2+q**2))))
  end function fdedr


end module mod_gxcpt
