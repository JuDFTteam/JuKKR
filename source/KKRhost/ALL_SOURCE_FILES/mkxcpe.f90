SUBROUTINE mkxcpe(nspin,ir,np,l1max,rv,rholm,vxcp,excp,thet,ylm,  &
        dylmt1,dylmt2,dylmf1,dylmf2,dylmtf,drrl,ddrrl,  &
        drrul,ddrrul,irmd,lmpotd)
!     ..
implicit none
!.. Parameters ..
      INTEGER IJD
      PARAMETER (IJD=434)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION RV
      INTEGER IR,IRMD,L1MAX,LMPOTD,NP,NSPIN
!..
!.. Array Arguments ..
DOUBLE PRECISION DDRRL(IRMD,LMPOTD),DDRRUL(IRMD,LMPOTD), &
                 DRRL(IRMD,LMPOTD),DRRUL(IRMD,LMPOTD), &
                 DYLMF1(IJD,LMPOTD),DYLMF2(IJD,LMPOTD), &
                 DYLMT1(IJD,LMPOTD),DYLMT2(IJD,LMPOTD), &
                 DYLMTF(IJD,LMPOTD),EXCP(IJD),RHOLM(LMPOTD,2), &
                 THET(IJD),VXCP(IJD,2),YLM(IJD,LMPOTD)
!..
!.. Local Scalars ..
DOUBLE PRECISION AGR,AGRD,AGRU,CEDG,CEDL,CHG,COSX,DAGRF,DAGRFD, &
                 DAGRFU,DAGRR,DAGRRD,DAGRRU,DAGRT,DAGRTD,DAGRTU, &
                 DDRR,DDRRD,DDRRU,DF1,DF2,DRR,DRRD,DRRU,DT1,DT2, &
                 DTF,DZDFS,DZDR,DZDTR,ETOT0,ETOTA0,G2R,G2RD,G2RU, &
                 GGGR,GGGRD,GGGRU,GRF,GRFD,GRFU,GRGRD,GRGRU,GRR, &
                 GRRD,GRRU,GRT,GRTD,GRTU,GZGR,RDSPR,RO,ROD,ROU, &
                 RV2,RV3,RVSIN1,RVSIN2,RY2,RYLM,SINT1,SINT2,SMAG, &
                 SPI,TANT1,VCG1,VCG2,VCL1,VCL2,VTOT1,VTOT2,VTOTA1, &
                 VTOTA2,VXG1,VXG2,VXL1,VXL2,XEDG,XEDL,ZTA
INTEGER IDSPR,IM,IP,L1,LL,LLMAX,LM,LMAX,NN,NN1
!..
!.. Local Arrays ..
DOUBLE PRECISION DDRY(IJD),DDRYD(IJD),DDRYU(IJD),DRDF(IJD), &
                 DRDFD(IJD),DRDFU(IJD),DRDT(IJD),DRDTD(IJD), &
                 DRDTU(IJD),DRY(IJD),DRYD(IJD),DRYU(IJD), &
                 RDF1(IJD),RDF1D(IJD),RDF1U(IJD),RDF2(IJD), &
                 RDF2D(IJD),RDF2U(IJD),RDT1(IJD),RDT1D(IJD), &
                 RDT1U(IJD),RDT2(IJD),RDT2D(IJD),RDT2U(IJD), &
                 RDTF(IJD),RDTFD(IJD),RDTFU(IJD),RY(IJD),RYD(IJD), &
                 RYU(IJD)
!..
!.. External Subroutines ..
      EXTERNAL GXCPT
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,COS,MAX,MIN,SIGN,SIN,SQRT,TAN
!..
!.. Data statements ..
      DATA RDSPR/9.0d0/
!..
!.....------------------------------------------------------------------
!     rl: charge=sumlm(rl*ylm)
!     ry=sumlm(ro*ylm), dry=sumlm(drr*ylm), ddry=sumlm(ddrr*ylm),
!c    rdt1=sumlm(ro*dylmt1), rdt2=sumlm(ro*dylmt2), ...
!c    rdf1=sumlm(ro*dylmf1), rdf2=sumlm(ro*dylmf2), ...
!c    rdtf=sumlm(ro*dylmtf), rdf2=sumlm(ro*dylmf2), ...
!c    drdt=sumlm(drr*dylmt1),drdf=sumlm(drr*dylmf1),

!     agr: abs(grad(ro)), g2r: laplacian(ro),
!c    gggr: grad(ro)*grad(agr),
!c    grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro).

!     dagrr,-t,-f: d(agr)/dr, d(agr)/dth, d(agr)/dfi.
!.....------------------------------------------------------------------
!     if(meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or.
!    &   mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD) then
!       write(6,'(/'' meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or. '',
!    &    ''mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD.''/
!    &    '' meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD='',
!    & 10i4)') meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD
!       stop14
!     endif
!heck  ist=mesh
llmax = l1max*l1max
lmax = l1max - 1
!     lmax2=lmax*2
!     llmax2=(lmax2+1)**2
!     lmax3=lmax*1
!     llmax3=(lmax3+1)**2





!     write(6,9030) (ii,drrs(ii,1),ddrrs(ii,1),drrus(ii,1),
!    &               ddrrus(ii,1),ii=ir,ir)
!9030 format(1x,' ist drrs  ddrrs drrus ddrrus',i5,4e12.5)


DO  ip = 1,np
  
  ry(ip) = 0.d0
  dry(ip) = 0.d0
  ddry(ip) = 0.d0
  rdt1(ip) = 0.d0
  rdt2(ip) = 0.d0
  rdf1(ip) = 0.d0
  rdf2(ip) = 0.d0
  rdtf(ip) = 0.d0
  drdt(ip) = 0.d0
  drdf(ip) = 0.d0
  
  ryu(ip) = 0.d0
  dryu(ip) = 0.d0
  ddryu(ip) = 0.d0
  rdt1u(ip) = 0.d0
  rdt2u(ip) = 0.d0
  rdf1u(ip) = 0.d0
  rdf2u(ip) = 0.d0
  rdtfu(ip) = 0.d0
  drdtu(ip) = 0.d0
  drdfu(ip) = 0.d0
  
  ryd(ip) = 0.d0
  dryd(ip) = 0.d0
  ddryd(ip) = 0.d0
  rdt1d(ip) = 0.d0
  rdt2d(ip) = 0.d0
  rdf1d(ip) = 0.d0
  rdf2d(ip) = 0.d0
  rdtfd(ip) = 0.d0
  drdtd(ip) = 0.d0
  drdfd(ip) = 0.d0
  
END DO

!     write(6,'(/'' nspin,mesh,np,l1max='',4i5)') nspin,mesh,np,l1max

lm = 0

DO  l1 = 1,l1max
  
  ll = l1 - 1
  
  
  
  DO  im = -ll,ll
    
    lm = lm + 1
    
    
    ro = rholm(lm,1)*2.d0
    rou = ro/2.d0
    rod = rou
    
    IF (nspin /= 1) THEN
      
      ro = rholm(lm,1) + rholm(lm,2)
      rou = rholm(lm,2)
      rod = rholm(lm,1)
!        write(6,9001) ro,rou,rod
    END IF
    drr = drrl(ir,lm)
    ddrr = ddrrl(ir,lm)
    drru = drrul(ir,lm)
    ddrru = ddrrul(ir,lm)
    drrd = drr - drru
    ddrrd = ddrr - ddrru
    
    
    DO  ip = 1,np
      
      rylm = ylm(ip,lm)
      dt1 = dylmt1(ip,lm)
      dt2 = dylmt2(ip,lm)
      df1 = dylmf1(ip,lm)
      df2 = dylmf2(ip,lm)
      dtf = dylmtf(ip,lm)
      
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
!rc             if (ip.eq.5.or.ip.eq.6) then
!             write(6,9907) lm,rylm,dt1,dt2,df1,df2,dtf
!9907         format(1x,' lmt ',i3,6d12.6)
!             write(6,*) 'nikos',dry(ip),ddry(ip)
!             end if
      
      
    END DO
  END DO
END DO



DO  ip = 1,np
  sint1 = SIN(thet(ip))
  sint2 = sint1**2
  tant1 = TAN(thet(ip))
  IF (sint1 == 0.d0) THEN
    vxcp(ip,1) = 0.d0
    vxcp(ip,2) = 0.d0
    excp(ip) = 0.d0
!          WRITE (6,FMT=*) 'interpolate'
    
! set values later
  ELSE
    rv2 = rv**2
    rv3 = rv**3
    
    
    rvsin1 = rv*sint1
    rvsin2 = rv2*sint2
    
    grr = dry(ip)
    grt = rdt1(ip)/rv
    grf = rdf1(ip)/rvsin1
    ry2 = ry(ip)**2
    
    agr = SQRT(grr**2+grt**2+grf**2)
    
    dagrr = (dry(ip)*ddry(ip)*rv3+ rdt1(ip)* (drdt(ip)*rv-rdt1(ip))+  &
        rdf1(ip)* (drdf(ip)*rv-rdf1(ip))/sint2)/agr/rv3
    
    dagrt = (dry(ip)*drdt(ip)*rv2+rdt1(ip)*rdt2(ip)+  &
        rdf1(ip)* (-rdf1(ip)/tant1+rdtf(ip))/sint2)/ (agr*rv3)
    
    dagrf = (dry(ip)*drdf(ip)*rv2+rdt1(ip)*rdtf(ip)+  &
        rdf1(ip)*rdf2(ip)/sint2)/ (agr*rv3*sint1)
    
    
    dzdr = ((dryu(ip)-dryd(ip))*ry(ip)- (ryu(ip)-ryd(ip))*dry(ip))/ry2
    
    dzdtr = ((rdt1u(ip)-rdt1d(ip))*ry(ip)- (ryu(ip)-ryd(ip))*rdt1(ip))/ry2/rv
    
    dzdfs = ((rdf1u(ip)-rdf1d(ip))*ry(ip)-  &
        (ryu(ip)-ryd(ip))*rdf1(ip))/ry2/rvsin1
    
    g2r = ddry(ip) + 2.d0*dry(ip)/rv +  &
        (rdt2(ip)+rdt1(ip)/tant1+rdf2(ip)/sint2)/rv2
    
    gggr = grr*dagrr + grt*dagrt + grf*dagrf
    
    gzgr = dzdr*grr + dzdtr*grt + dzdfs*grf
    
    
    chg = ry(ip)
    spi = ryu(ip) - ryd(ip)
    chg = MAX(1.0D-12,chg)
    smag = SIGN(1.0D0,spi)
    spi = smag*MIN(chg-1.0D-12,ABS(spi))
    zta = spi/chg
    
    
    grru = dryu(ip)
    grtu = rdt1u(ip)/rv
    grfu = rdf1u(ip)/rvsin1
    
    agru = SQRT(grru**2+grtu**2+grfu**2)
    
    dagrru = (dryu(ip)*ddryu(ip)*rv3+ rdt1u(ip)* (drdtu(ip)*rv-rdt1u(ip))+  &
        rdf1u(ip)* (drdfu(ip)*rv-rdf1u(ip))/sint2)/agru/rv3
    
    dagrtu = (dryu(ip)*drdtu(ip)*rv2+rdt1u(ip)*rdt2u(ip)+  &
        rdf1u(ip)* (-rdf1u(ip)/tant1+rdtfu(ip))/sint2)/ (agru*rv3)
    
    dagrfu = (dryu(ip)*drdfu(ip)*rv2+rdt1u(ip)*rdtfu(ip)+  &
        rdf1u(ip)*rdf2u(ip)/sint2)/ (agru*rv3*sint1)
    
    
    
    g2ru = ddryu(ip) + 2.d0*dryu(ip)/rv +  &
        (rdt2u(ip)+rdt1u(ip)/tant1+rdf2u(ip)/sint2)/rv2
    
    gggru = grru*dagrru + grtu*dagrtu + grfu*dagrfu
    
    grgru = grr*grru + grt*grtu + grf*grfu
    
    
    grrd = dryd(ip)
    grtd = rdt1d(ip)/rv
    grfd = rdf1d(ip)/rvsin1
    
    agrd = SQRT(grrd**2+grtd**2+grfd**2)
    
    dagrrd = (dryd(ip)*ddryd(ip)*rv3+ rdt1d(ip)* (drdtd(ip)*rv-rdt1d(ip))+  &
        rdf1d(ip)* (drdfd(ip)*rv-rdf1d(ip))/sint2)/agrd/rv3
    
    dagrtd = (dryd(ip)*drdtd(ip)*rv2+rdt1d(ip)*rdt2d(ip)+  &
        rdf1d(ip)* (-rdf1d(ip)/tant1+rdtfd(ip))/sint2)/ (agrd*rv3)
    
    dagrfd = (dryd(ip)*drdfd(ip)*rv2+rdt1d(ip)*rdtfd(ip)+  &
        rdf1d(ip)*rdf2d(ip)/sint2)/ (agrd*rv3*sint1)
    
    
    
    g2rd = ddryd(ip) + 2.d0*dryd(ip)/rv +  &
        (rdt2d(ip)+rdt1d(ip)/tant1+rdf2d(ip)/sint2)/rv2
    
    gggrd = grrd*dagrrd + grtd*dagrtd + grfd*dagrfd
    
    grgrd = grr*grrd + grt*grtd + grf*grfd
    
    
    idspr = 0
    IF (rv > rdspr) idspr = 1
    
    
    
    
    
! for debug
    CALL gxcpt(idspr,chg,zta,agr,agru,agrd,g2r,g2ru,g2rd,gggr,  &
        gggru,gggrd,grgru,grgrd,gzgr,vxcp(ip,2),vxcp(ip,1),  &
        excp(ip),vxl1,vxl2,vcl1,vcl2,xedl,cedl,vxg1,vxg2, vcg1,vcg2,xedg,cedg)
    
    
!     if(ip.eq.202) then
!     write(6,9912) ir,ip,ry(ip),zta,vxcp(ip,2),vxcp(ip,1)
!9912 format(1x,' ir ip ry zta',2i5,5e15.6)
!     write(6,*) 'mkxcpe',sint1,sint2,tant1,thet(ip)
    
!     write(6,9911) agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,
!    &              gggrd,grgru,grgrd,gzgr
!9911 format(1x,' agr  ',6E15.6)
    
!     write(6,7777) vxcp(ip,1),vxcp(ip,2),
!    &              vxl1,vxl2,vcl1,vcl2
!     write(6,7778) vxg1,vxg2,vcg1,vcg2
!7777 format(1x,'vxcp(1,2) vxl(1,2) vcl(1,2) ',6D15.6)
!7778 format(1x,'vxg(1,2) vcg(1,2)  (asada)  ',4D15.6)
    
!     end if
  END IF
  
END DO

! This is expected to work only for the Lebedev points
nn = 0
nn1 = 0
vtot1 = 0.d0
vtot2 = 0.d0
etot0 = 0.d0
vtota1 = 0.d0
vtota2 = 0.d0
etota0 = 0.d0
DO  ip = 1,np
  cosx = COS(thet(ip))
  IF (cosx > 0.99D0 .AND. cosx /= 1.d0) THEN
    nn = nn + 1
    vtot1 = vtot1 + vxcp(ip,1)
    vtot2 = vtot2 + vxcp(ip,2)
    etot0 = etot0 + excp(ip)
!        write(6,*) 'more',ip,vxcp(ip,1),nn
  END IF
  IF (cosx < -0.99D0 .AND. cosx /= -1.d0) THEN
    nn1 = nn1 + 1
    vtota1 = vtota1 + vxcp(ip,1)
    vtota2 = vtota2 + vxcp(ip,2)
    etota0 = etota0 + excp(ip)
!           write(6,*) 'less',ip,vxcp(ip,1),nn1
  END IF
END DO
DO  ip = 1,np
  cosx = COS(thet(ip))
  IF (cosx == 1.d0) THEN
    vxcp(ip,1) = vtot1/nn
    vxcp(ip,2) = vtot2/nn
    excp(ip) = etot0/nn
!     write(6,*) 'averaging ',ip,vxcp(ip,1),vxcp(ip,2),excp(ip)
!     write(6,*) 'averaging1 ',vtot1,vtot2,etot0,nn
  END IF
  IF (cosx == -1.d0) THEN
    vxcp(ip,1) = vtota1/nn1
    vxcp(ip,2) = vtota2/nn1
    excp(ip) = etota0/nn1
!     write(6,*) 'averaging ',ip,cosx,vxcp(ip,1),vxcp(ip,2),excp(ip)
!     write(6,*)'averaging2 ',vtota1,vtota2,etota0,nn1
  END IF
END DO
RETURN
END SUBROUTINE mkxcpe
