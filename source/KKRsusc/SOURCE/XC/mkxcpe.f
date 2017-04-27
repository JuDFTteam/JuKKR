      SUBROUTINE MKXCPE(NSPIN,IR,NP,L1MAX,RV,RHOLM,VXCP,EXCP,THET,YLM,
     +                  DYLMT1,DYLMT2,DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,
     +                  DRRUL,DDRRUL,IRMD,LMPOTD)
C     .. Parameters ..
      INTEGER IJD
      PARAMETER (IJD=434)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION RV
      INTEGER IR,IRMD,L1MAX,LMPOTD,NP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DDRRL(IRMD,LMPOTD),DDRRUL(IRMD,LMPOTD),
     +                 DRRL(IRMD,LMPOTD),DRRUL(IRMD,LMPOTD),
     +                 DYLMF1(IJD,LMPOTD),DYLMF2(IJD,LMPOTD),
     +                 DYLMT1(IJD,LMPOTD),DYLMT2(IJD,LMPOTD),
     +                 DYLMTF(IJD,LMPOTD),EXCP(IJD),RHOLM(LMPOTD,2),
     +                 THET(IJD),VXCP(IJD,2),YLM(IJD,LMPOTD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AGR,AGRD,AGRU,CEDG,CEDL,CHG,COSX,DAGRF,DAGRFD,
     +                 DAGRFU,DAGRR,DAGRRD,DAGRRU,DAGRT,DAGRTD,DAGRTU,
     +                 DDRR,DDRRD,DDRRU,DF1,DF2,DRR,DRRD,DRRU,DT1,DT2,
     +                 DTF,DZDFS,DZDR,DZDTR,ETOT0,ETOTA0,G2R,G2RD,G2RU,
     +                 GGGR,GGGRD,GGGRU,GRF,GRFD,GRFU,GRGRD,GRGRU,GRR,
     +                 GRRD,GRRU,GRT,GRTD,GRTU,GZGR,RDSPR,RO,ROD,ROU,
     +                 RV2,RV3,RVSIN1,RVSIN2,RY2,RYLM,SINT1,SINT2,SMAG,
     +                 SPI,TANT1,VCG1,VCG2,VCL1,VCL2,VTOT1,VTOT2,VTOTA1,
     +                 VTOTA2,VXG1,VXG2,VXL1,VXL2,XEDG,XEDL,ZTA
      INTEGER IDSPR,IM,IP,L1,LL,LLMAX,LM,LMAX,NN,NN1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DDRY(IJD),DDRYD(IJD),DDRYU(IJD),DRDF(IJD),
     +                 DRDFD(IJD),DRDFU(IJD),DRDT(IJD),DRDTD(IJD),
     +                 DRDTU(IJD),DRY(IJD),DRYD(IJD),DRYU(IJD),
     +                 RDF1(IJD),RDF1D(IJD),RDF1U(IJD),RDF2(IJD),
     +                 RDF2D(IJD),RDF2U(IJD),RDT1(IJD),RDT1D(IJD),
     +                 RDT1U(IJD),RDT2(IJD),RDT2D(IJD),RDT2U(IJD),
     +                 RDTF(IJD),RDTFD(IJD),RDTFU(IJD),RY(IJD),RYD(IJD),
     +                 RYU(IJD)
C     ..
C     .. External Subroutines ..
!       EXTERNAL GXCPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,MAX,MIN,SIGN,SIN,SQRT,TAN
C     ..
C     .. Data statements ..
      DATA RDSPR/9.0d0/
C     ..
c.....------------------------------------------------------------------
c     rl: charge=sumlm(rl*ylm)
c     ry=sumlm(ro*ylm), dry=sumlm(drr*ylm), ddry=sumlm(ddrr*ylm),
cc    rdt1=sumlm(ro*dylmt1), rdt2=sumlm(ro*dylmt2), ...
cc    rdf1=sumlm(ro*dylmf1), rdf2=sumlm(ro*dylmf2), ...
cc    rdtf=sumlm(ro*dylmtf), rdf2=sumlm(ro*dylmf2), ...
cc    drdt=sumlm(drr*dylmt1),drdf=sumlm(drr*dylmf1),

c     agr: abs(grad(ro)), g2r: laplacian(ro),
cc    gggr: grad(ro)*grad(agr),
cc    grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro).

c     dagrr,-t,-f: d(agr)/dr, d(agr)/dth, d(agr)/dfi.
c.....------------------------------------------------------------------
c     if(meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or.
c    &   mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD) then
c       write(6,'(/'' meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or. '',
c    &    ''mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD.''/
c    &    '' meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD='',
c    & 10i4)') meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD
c       stop14
c     endif
check  ist=mesh
      LLMAX = L1MAX*L1MAX
      LMAX = L1MAX - 1
c     lmax2=lmax*2
c     llmax2=(lmax2+1)**2
c     lmax3=lmax*1
c     llmax3=(lmax3+1)**2
c
c
c
c
c
c     write(6,9030) (ii,drrs(ii,1),ddrrs(ii,1),drrus(ii,1),
c    &               ddrrus(ii,1),ii=ir,ir)
c9030 format(1x,' ist drrs  ddrrs drrus ddrrus',i5,4e12.5)
c

      DO 10 IP = 1,NP

        RY(IP) = 0.D0
        DRY(IP) = 0.D0
        DDRY(IP) = 0.D0
        RDT1(IP) = 0.D0
        RDT2(IP) = 0.D0
        RDF1(IP) = 0.D0
        RDF2(IP) = 0.D0
        RDTF(IP) = 0.D0
        DRDT(IP) = 0.D0
        DRDF(IP) = 0.D0

        RYU(IP) = 0.D0
        DRYU(IP) = 0.D0
        DDRYU(IP) = 0.D0
        RDT1U(IP) = 0.D0
        RDT2U(IP) = 0.D0
        RDF1U(IP) = 0.D0
        RDF2U(IP) = 0.D0
        RDTFU(IP) = 0.D0
        DRDTU(IP) = 0.D0
        DRDFU(IP) = 0.D0

        RYD(IP) = 0.D0
        DRYD(IP) = 0.D0
        DDRYD(IP) = 0.D0
        RDT1D(IP) = 0.D0
        RDT2D(IP) = 0.D0
        RDF1D(IP) = 0.D0
        RDF2D(IP) = 0.D0
        RDTFD(IP) = 0.D0
        DRDTD(IP) = 0.D0
        DRDFD(IP) = 0.D0

   10 CONTINUE

c     write(6,'(/'' nspin,mesh,np,l1max='',4i5)') nspin,mesh,np,l1max

      LM = 0

      DO 40 L1 = 1,L1MAX

        LL = L1 - 1



        DO 30 IM = -LL,LL

          LM = LM + 1


          RO = RHOLM(LM,1)*2.d0
          ROU = RO/2.D0
          ROD = ROU
c
          IF (NSPIN.NE.1) THEN
c
            RO = RHOLM(LM,1) + RHOLM(LM,2)
            ROU = RHOLM(LM,2)
            ROD = RHOLM(LM,1)
c        write(6,9001) ro,rou,rod
          END IF
          DRR = DRRL(IR,LM)
          DDRR = DDRRL(IR,LM)
          DRRU = DRRUL(IR,LM)
          DDRRU = DDRRUL(IR,LM)
          DRRD = DRR - DRRU
          DDRRD = DDRR - DDRRU


          DO 20 IP = 1,NP

            RYLM = YLM(IP,LM)
            DT1 = DYLMT1(IP,LM)
            DT2 = DYLMT2(IP,LM)
            DF1 = DYLMF1(IP,LM)
            DF2 = DYLMF2(IP,LM)
            DTF = DYLMTF(IP,LM)

            RY(IP) = RY(IP) + RO*RYLM
            DRY(IP) = DRY(IP) + DRR*RYLM
            DDRY(IP) = DDRY(IP) + DDRR*RYLM

            RYU(IP) = RYU(IP) + ROU*RYLM
            DRYU(IP) = DRYU(IP) + DRRU*RYLM
            DDRYU(IP) = DDRYU(IP) + DDRRU*RYLM

            RYD(IP) = RYD(IP) + ROD*RYLM
            DRYD(IP) = DRYD(IP) + DRRD*RYLM
            DDRYD(IP) = DDRYD(IP) + DDRRD*RYLM

            RDT1(IP) = RDT1(IP) + RO*DT1
            RDT2(IP) = RDT2(IP) + RO*DT2
            RDF1(IP) = RDF1(IP) + RO*DF1
            RDF2(IP) = RDF2(IP) + RO*DF2
            RDTF(IP) = RDTF(IP) + RO*DTF
            DRDT(IP) = DRDT(IP) + DRR*DT1
            DRDF(IP) = DRDF(IP) + DRR*DF1

            RDT1U(IP) = RDT1U(IP) + ROU*DT1
            RDT2U(IP) = RDT2U(IP) + ROU*DT2
            RDF1U(IP) = RDF1U(IP) + ROU*DF1
            RDF2U(IP) = RDF2U(IP) + ROU*DF2
            RDTFU(IP) = RDTFU(IP) + ROU*DTF
            DRDTU(IP) = DRDTU(IP) + DRRU*DT1
            DRDFU(IP) = DRDFU(IP) + DRRU*DF1

            RDT1D(IP) = RDT1D(IP) + ROD*DT1
            RDT2D(IP) = RDT2D(IP) + ROD*DT2
            RDF1D(IP) = RDF1D(IP) + ROD*DF1
            RDF2D(IP) = RDF2D(IP) + ROD*DF2
            RDTFD(IP) = RDTFD(IP) + ROD*DTF
            DRDTD(IP) = DRDTD(IP) + DRRD*DT1
            DRDFD(IP) = DRDFD(IP) + DRRD*DF1
crc             if (ip.eq.5.or.ip.eq.6) then
c             write(6,9907) lm,rylm,dt1,dt2,df1,df2,dtf
c9907         format(1x,' lmt ',i3,6d12.6)
c             write(6,*) 'nikos',dry(ip),ddry(ip)
c             end if


   20     CONTINUE
   30   CONTINUE
   40 CONTINUE

c

      DO 50 IP = 1,NP
        SINT1 = SIN(THET(IP))
        SINT2 = SINT1**2
        TANT1 = TAN(THET(IP))
        IF (SINT1.EQ.0.d0) THEN
          VXCP(IP,1) = 0.d0
          VXCP(IP,2) = 0.d0
          EXCP(IP) = 0.d0
c          WRITE (6,FMT=*) 'interpolate'
c
c set values later
        ELSE
          RV2 = RV**2
          RV3 = RV**3


          RVSIN1 = RV*SINT1
          RVSIN2 = RV2*SINT2

          GRR = DRY(IP)
          GRT = RDT1(IP)/RV
          GRF = RDF1(IP)/RVSIN1
          RY2 = RY(IP)**2

          AGR = SQRT(GRR**2+GRT**2+GRF**2)

          DAGRR = (DRY(IP)*DDRY(IP)*RV3+
     +            RDT1(IP)* (DRDT(IP)*RV-RDT1(IP))+
     +            RDF1(IP)* (DRDF(IP)*RV-RDF1(IP))/SINT2)/AGR/RV3

          DAGRT = (DRY(IP)*DRDT(IP)*RV2+RDT1(IP)*RDT2(IP)+
     +            RDF1(IP)* (-RDF1(IP)/TANT1+RDTF(IP))/SINT2)/ (AGR*RV3)

          DAGRF = (DRY(IP)*DRDF(IP)*RV2+RDT1(IP)*RDTF(IP)+
     +            RDF1(IP)*RDF2(IP)/SINT2)/ (AGR*RV3*SINT1)


          DZDR = ((DRYU(IP)-DRYD(IP))*RY(IP)-
     +           (RYU(IP)-RYD(IP))*DRY(IP))/RY2

          DZDTR = ((RDT1U(IP)-RDT1D(IP))*RY(IP)-
     +            (RYU(IP)-RYD(IP))*RDT1(IP))/RY2/RV

          DZDFS = ((RDF1U(IP)-RDF1D(IP))*RY(IP)-
     +            (RYU(IP)-RYD(IP))*RDF1(IP))/RY2/RVSIN1

          G2R = DDRY(IP) + 2.D0*DRY(IP)/RV +
     +          (RDT2(IP)+RDT1(IP)/TANT1+RDF2(IP)/SINT2)/RV2

          GGGR = GRR*DAGRR + GRT*DAGRT + GRF*DAGRF

          GZGR = DZDR*GRR + DZDTR*GRT + DZDFS*GRF

c
          CHG = RY(IP)
          SPI = RYU(IP) - RYD(IP)
          CHG = MAX(1.0D-12,CHG)
          SMAG = SIGN(1.0D0,SPI)
          SPI = SMAG*MIN(CHG-1.0D-12,ABS(SPI))
          ZTA = SPI/CHG
c

          GRRU = DRYU(IP)
          GRTU = RDT1U(IP)/RV
          GRFU = RDF1U(IP)/RVSIN1

          AGRU = SQRT(GRRU**2+GRTU**2+GRFU**2)

          DAGRRU = (DRYU(IP)*DDRYU(IP)*RV3+
     +             RDT1U(IP)* (DRDTU(IP)*RV-RDT1U(IP))+
     +             RDF1U(IP)* (DRDFU(IP)*RV-RDF1U(IP))/SINT2)/AGRU/RV3

          DAGRTU = (DRYU(IP)*DRDTU(IP)*RV2+RDT1U(IP)*RDT2U(IP)+
     +             RDF1U(IP)* (-RDF1U(IP)/TANT1+RDTFU(IP))/SINT2)/
     +             (AGRU*RV3)

          DAGRFU = (DRYU(IP)*DRDFU(IP)*RV2+RDT1U(IP)*RDTFU(IP)+
     +             RDF1U(IP)*RDF2U(IP)/SINT2)/ (AGRU*RV3*SINT1)



          G2RU = DDRYU(IP) + 2.D0*DRYU(IP)/RV +
     +           (RDT2U(IP)+RDT1U(IP)/TANT1+RDF2U(IP)/SINT2)/RV2

          GGGRU = GRRU*DAGRRU + GRTU*DAGRTU + GRFU*DAGRFU

          GRGRU = GRR*GRRU + GRT*GRTU + GRF*GRFU


          GRRD = DRYD(IP)
          GRTD = RDT1D(IP)/RV
          GRFD = RDF1D(IP)/RVSIN1

          AGRD = SQRT(GRRD**2+GRTD**2+GRFD**2)

          DAGRRD = (DRYD(IP)*DDRYD(IP)*RV3+
     +             RDT1D(IP)* (DRDTD(IP)*RV-RDT1D(IP))+
     +             RDF1D(IP)* (DRDFD(IP)*RV-RDF1D(IP))/SINT2)/AGRD/RV3

          DAGRTD = (DRYD(IP)*DRDTD(IP)*RV2+RDT1D(IP)*RDT2D(IP)+
     +             RDF1D(IP)* (-RDF1D(IP)/TANT1+RDTFD(IP))/SINT2)/
     +             (AGRD*RV3)

          DAGRFD = (DRYD(IP)*DRDFD(IP)*RV2+RDT1D(IP)*RDTFD(IP)+
     +             RDF1D(IP)*RDF2D(IP)/SINT2)/ (AGRD*RV3*SINT1)



          G2RD = DDRYD(IP) + 2.D0*DRYD(IP)/RV +
     +           (RDT2D(IP)+RDT1D(IP)/TANT1+RDF2D(IP)/SINT2)/RV2

          GGGRD = GRRD*DAGRRD + GRTD*DAGRTD + GRFD*DAGRFD

          GRGRD = GRR*GRRD + GRT*GRTD + GRF*GRFD


          IDSPR = 0
          IF (RV.GT.RDSPR) IDSPR = 1

c
c


c for debug
          CALL GXCPT(IDSPR,CHG,ZTA,AGR,AGRU,AGRD,G2R,G2RU,G2RD,GGGR,
     +               GGGRU,GGGRD,GRGRU,GRGRD,GZGR,VXCP(IP,2),VXCP(IP,1),
     +               EXCP(IP),VXL1,VXL2,VCL1,VCL2,XEDL,CEDL,VXG1,VXG2,
     +               VCG1,VCG2,XEDG,CEDG)
c

c     if(ip.eq.202) then
c     write(6,9912) ir,ip,ry(ip),zta,vxcp(ip,2),vxcp(ip,1)
c9912 format(1x,' ir ip ry zta',2i5,5e15.6)
c     write(6,*) 'mkxcpe',sint1,sint2,tant1,thet(ip)
c
c     write(6,9911) agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,
c    &              gggrd,grgru,grgrd,gzgr
c9911 format(1x,' agr  ',6E15.6)

c     write(6,7777) vxcp(ip,1),vxcp(ip,2),
c    &              vxl1,vxl2,vcl1,vcl2
c     write(6,7778) vxg1,vxg2,vcg1,vcg2
c7777 format(1x,'vxcp(1,2) vxl(1,2) vcl(1,2) ',6D15.6)
c7778 format(1x,'vxg(1,2) vcg(1,2)  (asada)  ',4D15.6)

c     end if
        END IF

   50 CONTINUE
c
c This is expected to work only for the Lebedev points
      NN = 0
      NN1 = 0
      VTOT1 = 0.d0
      VTOT2 = 0.d0
      ETOT0 = 0.d0
      VTOTA1 = 0.d0
      VTOTA2 = 0.d0
      ETOTA0 = 0.d0
      DO 60 IP = 1,NP
        COSX = COS(THET(IP))
        IF (COSX.GT.0.99d0 .AND. COSX.NE.1.d0) THEN
          NN = NN + 1
          VTOT1 = VTOT1 + VXCP(IP,1)
          VTOT2 = VTOT2 + VXCP(IP,2)
          ETOT0 = ETOT0 + EXCP(IP)
c        write(6,*) 'more',ip,vxcp(ip,1),nn
        END IF
        IF (COSX.LT.-0.99d0 .AND. COSX.NE.-1.d0) THEN
          NN1 = NN1 + 1
          VTOTA1 = VTOTA1 + VXCP(IP,1)
          VTOTA2 = VTOTA2 + VXCP(IP,2)
          ETOTA0 = ETOTA0 + EXCP(IP)
c           write(6,*) 'less',ip,vxcp(ip,1),nn1
        END IF
   60 CONTINUE
      DO 70 IP = 1,NP
        COSX = COS(THET(IP))
        IF (COSX.EQ.1.d0) THEN
          VXCP(IP,1) = VTOT1/NN
          VXCP(IP,2) = VTOT2/NN
          EXCP(IP) = ETOT0/NN
c     write(6,*) 'averaging ',ip,vxcp(ip,1),vxcp(ip,2),excp(ip)
c     write(6,*) 'averaging1 ',vtot1,vtot2,etot0,nn
        END IF
        IF (COSX.EQ.-1.d0) THEN
          VXCP(IP,1) = VTOTA1/NN1
          VXCP(IP,2) = VTOTA2/NN1
          EXCP(IP) = ETOTA0/NN1
c     write(6,*) 'averaging ',ip,cosx,vxcp(ip,1),vxcp(ip,2),excp(ip)
c     write(6,*)'averaging2 ',vtota1,vtota2,etota0,nn1
        END IF
   70 CONTINUE
      RETURN
      END SUBROUTINE
