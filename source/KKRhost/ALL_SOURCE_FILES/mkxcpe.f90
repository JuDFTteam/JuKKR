module mod_mkxcpe

contains

subroutine mkxcpe(nspin, ir, np, l1max, rv, rholm, vxcp, excp, thet, ylm, &
  dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, &
  lmpotd)
  use :: mod_datatypes, only: dp
   use mod_gxcpt
  ! ..
  implicit none
  ! .. Parameters ..
  real (kind=dp), parameter :: eps=1.0D-12
  integer :: ijd
  parameter (ijd=434)
  ! ..
  ! .. Scalar Arguments ..
  real (kind=dp) :: rv
  integer :: ir, irmd, l1max, lmpotd, np, nspin
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: ddrrl(irmd, lmpotd), ddrrul(irmd, lmpotd), &
    drrl(irmd, lmpotd), drrul(irmd, lmpotd), dylmf1(ijd, lmpotd), &
    dylmf2(ijd, lmpotd), dylmt1(ijd, lmpotd), dylmt2(ijd, lmpotd), &
    dylmtf(ijd, lmpotd), excp(ijd), rholm(lmpotd, 2), thet(ijd), vxcp(ijd, 2), &
    ylm(ijd, lmpotd)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: agr, agrd, agru, cedg, cedl, chg, cosx, dagrf, dagrfd, &
    dagrfu, dagrr, dagrrd, dagrru, dagrt, dagrtd, dagrtu, ddrr, ddrrd, ddrru, &
    df1, df2, drr, drrd, drru, dt1, dt2, dtf, dzdfs, dzdr, dzdtr, etot0, &
    etota0, g2r, g2rd, g2ru, gggr, gggrd, gggru, grf, grfd, grfu, grgrd, &
    grgru, grr, grrd, grru, grt, grtd, grtu, gzgr, rdspr, ro, rod, rou, rv2, &
    rv3, rvsin1, rvsin2, ry2, rylm, sint1, sint2, smag, spi, tant1, vcg1, &
    vcg2, vcl1, vcl2, vtot1, vtot2, vtota1, vtota2, vxg1, vxg2, vxl1, vxl2, &
    xedg, xedl, zta
  integer :: idspr, im, ip, l1, ll, llmax, lm, lmax, nn, nn1
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: ddry(ijd), ddryd(ijd), ddryu(ijd), drdf(ijd), drdfd(ijd), &
    drdfu(ijd), drdt(ijd), drdtd(ijd), drdtu(ijd), dry(ijd), dryd(ijd), &
    dryu(ijd), rdf1(ijd), rdf1d(ijd), rdf1u(ijd), rdf2(ijd), rdf2d(ijd), &
    rdf2u(ijd), rdt1(ijd), rdt1d(ijd), rdt1u(ijd), rdt2(ijd), rdt2d(ijd), &
    rdt2u(ijd), rdtf(ijd), rdtfd(ijd), rdtfu(ijd), ry(ijd), ryd(ijd), ryu(ijd)
  ! ..
  ! .. External Subroutines ..
  external :: gxcpt
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: abs, cos, max, min, sign, sin, sqrt, tan
  ! ..
  ! .. Data statements ..
  data rdspr/9.0e0_dp/
  ! ..
  ! .....------------------------------------------------------------------
  ! rl: charge=sumlm(rl*ylm)
  ! ry=sumlm(ro*ylm), dry=sumlm(drr*ylm), ddry=sumlm(ddrr*ylm),
  ! c    rdt1=sumlm(ro*dylmt1), rdt2=sumlm(ro*dylmt2), ...
  ! c    rdf1=sumlm(ro*dylmf1), rdf2=sumlm(ro*dylmf2), ...
  ! c    rdtf=sumlm(ro*dylmtf), rdf2=sumlm(ro*dylmf2), ...
  ! c    drdt=sumlm(drr*dylmt1),drdf=sumlm(drr*dylmf1),

  ! agr: abs(grad(ro)), g2r: laplacian(ro),
  ! c    gggr: grad(ro)*grad(agr),
  ! c    grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro).

  ! dagrr,-t,-f: d(agr)/dr, d(agr)/dth, d(agr)/dfi.
  ! .....------------------------------------------------------------------
  ! if(meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or.
  ! &   mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD) then
  ! write(6,'(/'' meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or. '',
  ! &    ''mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD.''/
  ! &    '' meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD='',
  ! & 10i4)') meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD
  ! stop14
  ! endif
  ! heck  ist=mesh
  llmax = l1max*l1max
  lmax = l1max - 1
  ! lmax2=lmax*2
  ! llmax2=(lmax2+1)**2
  ! lmax3=lmax*1
  ! llmax3=(lmax3+1)**2





  ! write(6,9030) (ii,drrs(ii,1),ddrrs(ii,1),drrus(ii,1),
  ! &               ddrrus(ii,1),ii=ir,ir)
  ! 9030 format(1x,' ist drrs  ddrrs drrus ddrrus',i5,4e12.5)


  do ip = 1, np

    ry(ip) = 0.e0_dp
    dry(ip) = 0.e0_dp
    ddry(ip) = 0.e0_dp
    rdt1(ip) = 0.e0_dp
    rdt2(ip) = 0.e0_dp
    rdf1(ip) = 0.e0_dp
    rdf2(ip) = 0.e0_dp
    rdtf(ip) = 0.e0_dp
    drdt(ip) = 0.e0_dp
    drdf(ip) = 0.e0_dp

    ryu(ip) = 0.e0_dp
    dryu(ip) = 0.e0_dp
    ddryu(ip) = 0.e0_dp
    rdt1u(ip) = 0.e0_dp
    rdt2u(ip) = 0.e0_dp
    rdf1u(ip) = 0.e0_dp
    rdf2u(ip) = 0.e0_dp
    rdtfu(ip) = 0.e0_dp
    drdtu(ip) = 0.e0_dp
    drdfu(ip) = 0.e0_dp

    ryd(ip) = 0.e0_dp
    dryd(ip) = 0.e0_dp
    ddryd(ip) = 0.e0_dp
    rdt1d(ip) = 0.e0_dp
    rdt2d(ip) = 0.e0_dp
    rdf1d(ip) = 0.e0_dp
    rdf2d(ip) = 0.e0_dp
    rdtfd(ip) = 0.e0_dp
    drdtd(ip) = 0.e0_dp
    drdfd(ip) = 0.e0_dp

  end do

  ! write(6,'(/'' nspin,mesh,np,l1max='',4i5)') nspin,mesh,np,l1max

  lm = 0

  do l1 = 1, l1max

    ll = l1 - 1



    do im = -ll, ll

      lm = lm + 1


      ro = rholm(lm, 1)*2.e0_dp
      rou = ro/2.e0_dp
      rod = rou

      if (nspin/=1) then

        ro = rholm(lm, 1) + rholm(lm, 2)
        rou = rholm(lm, 2)
        rod = rholm(lm, 1)
        ! write(6,9001) ro,rou,rod
      end if
      drr = drrl(ir, lm)
      ddrr = ddrrl(ir, lm)
      drru = drrul(ir, lm)
      ddrru = ddrrul(ir, lm)
      drrd = drr - drru
      ddrrd = ddrr - ddrru


      do ip = 1, np

        rylm = ylm(ip, lm)
        dt1 = dylmt1(ip, lm)
        dt2 = dylmt2(ip, lm)
        df1 = dylmf1(ip, lm)
        df2 = dylmf2(ip, lm)
        dtf = dylmtf(ip, lm)

        ry(ip) = ry(ip) + ro*rylm
        dry(ip) = dry(ip) + drr*rylm
        ddry(ip) = ddry(ip) + ddrr*rylm

        ryu(ip) = ryu(ip) + rou*rylm
        dryu(ip) = dryu(ip) + drru*rylm
        ddryu(ip) = ddryu(ip) + ddrru*rylm

        ryd(ip) = ryd(ip) + rod*rylm
        dryd(ip) = dryd(ip) + drrd*rylm
        ddryd(ip) = ddryd(ip) + ddrrd*rylm

        rdt1(ip) = rdt1(ip) + ro*dt1
        rdt2(ip) = rdt2(ip) + ro*dt2
        rdf1(ip) = rdf1(ip) + ro*df1
        rdf2(ip) = rdf2(ip) + ro*df2
        rdtf(ip) = rdtf(ip) + ro*dtf
        drdt(ip) = drdt(ip) + drr*dt1
        drdf(ip) = drdf(ip) + drr*df1

        rdt1u(ip) = rdt1u(ip) + rou*dt1
        rdt2u(ip) = rdt2u(ip) + rou*dt2
        rdf1u(ip) = rdf1u(ip) + rou*df1
        rdf2u(ip) = rdf2u(ip) + rou*df2
        rdtfu(ip) = rdtfu(ip) + rou*dtf
        drdtu(ip) = drdtu(ip) + drru*dt1
        drdfu(ip) = drdfu(ip) + drru*df1

        rdt1d(ip) = rdt1d(ip) + rod*dt1
        rdt2d(ip) = rdt2d(ip) + rod*dt2
        rdf1d(ip) = rdf1d(ip) + rod*df1
        rdf2d(ip) = rdf2d(ip) + rod*df2
        rdtfd(ip) = rdtfd(ip) + rod*dtf
        drdtd(ip) = drdtd(ip) + drrd*dt1
        drdfd(ip) = drdfd(ip) + drrd*df1
        ! rc             if (ip.eq.5.or.ip.eq.6) then
        ! write(6,9907) lm,rylm,dt1,dt2,df1,df2,dtf
        ! 9907         format(1x,' lmt ',i3,6d12.6)
        ! write(6,*) 'nikos',dry(ip),ddry(ip)
        ! end if


      end do
    end do
  end do



  do ip = 1, np
    sint1 = sin(thet(ip))
    sint2 = sint1**2
    tant1 = tan(thet(ip))
    if (abs(sint1)<eps) then
      vxcp(ip, 1) = 0.e0_dp
      vxcp(ip, 2) = 0.e0_dp
      excp(ip) = 0.e0_dp
      ! WRITE (6,FMT=*) 'interpolate'

      ! set values later
    else
      rv2 = rv**2
      rv3 = rv**3


      rvsin1 = rv*sint1
      rvsin2 = rv2*sint2

      grr = dry(ip)
      grt = rdt1(ip)/rv
      grf = rdf1(ip)/rvsin1
      ry2 = ry(ip)**2

      agr = sqrt(grr**2+grt**2+grf**2)

      dagrr = (dry(ip)*ddry(ip)*rv3+rdt1(ip)*(drdt(ip)*rv-rdt1( &
        ip))+rdf1(ip)*(drdf(ip)*rv-rdf1(ip))/sint2)/agr/rv3

      dagrt = (dry(ip)*drdt(ip)*rv2+rdt1(ip)*rdt2(ip)+rdf1(ip)*(-rdf1( &
        ip)/tant1+rdtf(ip))/sint2)/(agr*rv3)

      dagrf = (dry(ip)*drdf(ip)*rv2+rdt1(ip)*rdtf(ip)+rdf1(ip)*rdf2(ip)/sint2) &
        /(agr*rv3*sint1)


      dzdr = ((dryu(ip)-dryd(ip))*ry(ip)-(ryu(ip)-ryd(ip))*dry(ip))/ry2

      dzdtr = ((rdt1u(ip)-rdt1d(ip))*ry(ip)-(ryu(ip)-ryd(ip))*rdt1(ip))/ry2/rv

      dzdfs = ((rdf1u(ip)-rdf1d(ip))*ry(ip)-(ryu(ip)-ryd(ip))*rdf1(ip))/ry2/ &
        rvsin1

      g2r = ddry(ip) + 2.e0_dp*dry(ip)/rv + (rdt2(ip)+rdt1(ip)/tant1+rdf2(ip)/ &
        sint2)/rv2

      gggr = grr*dagrr + grt*dagrt + grf*dagrf

      gzgr = dzdr*grr + dzdtr*grt + dzdfs*grf


      chg = ry(ip)
      spi = ryu(ip) - ryd(ip)
      chg = max(1.0e-12_dp, chg)
      smag = sign(1.0e0_dp, spi)
      spi = smag*min(chg-1.0e-12_dp, abs(spi))
      zta = spi/chg


      grru = dryu(ip)
      grtu = rdt1u(ip)/rv
      grfu = rdf1u(ip)/rvsin1

      agru = sqrt(grru**2+grtu**2+grfu**2)

      dagrru = (dryu(ip)*ddryu(ip)*rv3+rdt1u(ip)*(drdtu(ip)*rv-rdt1u( &
        ip))+rdf1u(ip)*(drdfu(ip)*rv-rdf1u(ip))/sint2)/agru/rv3

      dagrtu = (dryu(ip)*drdtu(ip)*rv2+rdt1u(ip)*rdt2u(ip)+rdf1u(ip)*(-rdf1u( &
        ip)/tant1+rdtfu(ip))/sint2)/(agru*rv3)

      dagrfu = (dryu(ip)*drdfu(ip)*rv2+rdt1u(ip)*rdtfu(ip)+ &
        rdf1u(ip)*rdf2u(ip)/sint2)/(agru*rv3*sint1)



      g2ru = ddryu(ip) + 2.e0_dp*dryu(ip)/rv + (rdt2u(ip)+rdt1u(ip)/tant1+ &
        rdf2u(ip)/sint2)/rv2

      gggru = grru*dagrru + grtu*dagrtu + grfu*dagrfu

      grgru = grr*grru + grt*grtu + grf*grfu


      grrd = dryd(ip)
      grtd = rdt1d(ip)/rv
      grfd = rdf1d(ip)/rvsin1

      agrd = sqrt(grrd**2+grtd**2+grfd**2)

      dagrrd = (dryd(ip)*ddryd(ip)*rv3+rdt1d(ip)*(drdtd(ip)*rv-rdt1d( &
        ip))+rdf1d(ip)*(drdfd(ip)*rv-rdf1d(ip))/sint2)/agrd/rv3

      dagrtd = (dryd(ip)*drdtd(ip)*rv2+rdt1d(ip)*rdt2d(ip)+rdf1d(ip)*(-rdf1d( &
        ip)/tant1+rdtfd(ip))/sint2)/(agrd*rv3)

      dagrfd = (dryd(ip)*drdfd(ip)*rv2+rdt1d(ip)*rdtfd(ip)+ &
        rdf1d(ip)*rdf2d(ip)/sint2)/(agrd*rv3*sint1)



      g2rd = ddryd(ip) + 2.e0_dp*dryd(ip)/rv + (rdt2d(ip)+rdt1d(ip)/tant1+ &
        rdf2d(ip)/sint2)/rv2

      gggrd = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd

      grgrd = grr*grrd + grt*grtd + grf*grfd


      idspr = 0
      if (rv>rdspr) idspr = 1





      ! for debug
      call gxcpt(idspr, chg, zta, agr, agru, agrd, g2r, g2ru, g2rd, gggr, &
        gggru, gggrd, grgru, grgrd, gzgr, vxcp(ip,2), vxcp(ip,1), excp(ip), &
        vxl1, vxl2, vcl1, vcl2, xedl, cedl, vxg1, vxg2, vcg1, vcg2, xedg, &
        cedg)


      ! if(ip.eq.202) then
      ! write(6,9912) ir,ip,ry(ip),zta,vxcp(ip,2),vxcp(ip,1)
      ! 9912 format(1x,' ir ip ry zta',2i5,5e15.6)
      ! write(6,*) 'mkxcpe',sint1,sint2,tant1,thet(ip)

      ! write(6,9911) agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,
      ! &              gggrd,grgru,grgrd,gzgr
      ! 9911 format(1x,' agr  ',6E15.6)

      ! write(6,7777) vxcp(ip,1),vxcp(ip,2),
      ! &              vxl1,vxl2,vcl1,vcl2
      ! write(6,7778) vxg1,vxg2,vcg1,vcg2
      ! 7777 format(1x,'vxcp(1,2) vxl(1,2) vcl(1,2) ',6D15.6)
      ! 7778 format(1x,'vxg(1,2) vcg(1,2)  (asada)  ',4D15.6)

      ! end if
    end if

  end do

  ! This is expected to work only for the Lebedev points
  nn = 0
  nn1 = 0
  vtot1 = 0.e0_dp
  vtot2 = 0.e0_dp
  etot0 = 0.e0_dp
  vtota1 = 0.e0_dp
  vtota2 = 0.e0_dp
  etota0 = 0.e0_dp
  do ip = 1, np
    cosx = cos(thet(ip))
    if (cosx>0.99e0_dp .and. abs(cosx-1.e0_dp)>eps) then
      nn = nn + 1
      vtot1 = vtot1 + vxcp(ip, 1)
      vtot2 = vtot2 + vxcp(ip, 2)
      etot0 = etot0 + excp(ip)
      ! write(6,*) 'more',ip,vxcp(ip,1),nn
    end if
    if (cosx<-0.99e0_dp .and. abs(cosx+1.e0_dp)>eps) then
      nn1 = nn1 + 1
      vtota1 = vtota1 + vxcp(ip, 1)
      vtota2 = vtota2 + vxcp(ip, 2)
      etota0 = etota0 + excp(ip)
      ! write(6,*) 'less',ip,vxcp(ip,1),nn1
    end if
  end do
  do ip = 1, np
    cosx = cos(thet(ip))
    if (abs(cosx-1.e0_dp)<eps) then
      vxcp(ip, 1) = vtot1/nn
      vxcp(ip, 2) = vtot2/nn
      excp(ip) = etot0/nn
      ! write(6,*) 'averaging ',ip,vxcp(ip,1),vxcp(ip,2),excp(ip)
      ! write(6,*) 'averaging1 ',vtot1,vtot2,etot0,nn
    end if
    if (abs(cosx+1.e0_dp)<eps) then
      vxcp(ip, 1) = vtota1/nn1
      vxcp(ip, 2) = vtota2/nn1
      excp(ip) = etota0/nn1
      ! write(6,*) 'averaging ',ip,cosx,vxcp(ip,1),vxcp(ip,2),excp(ip)
      ! write(6,*)'averaging2 ',vtota1,vtota2,etota0,nn1
    end if
  end do
  return
end subroutine mkxcpe

end module mod_mkxcpe
