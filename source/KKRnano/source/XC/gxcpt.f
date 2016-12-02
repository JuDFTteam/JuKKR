      SUBROUTINE GXCPT(IDSPR,RO,ZTA,AGR,AGRU,AGRD,G2R,G2RU,G2RD,GGGR,
     +                 GGGRU,GGGRD,GRGRU,GRGRD,GZGR,XCPTU,XCPTD,XCED,
     +                 VXLU,VXLD,VCLU,VCLD,XEDL,CEDL,VXGU,VXGD,VCGU,
     +                 VCGD,XEDG,CEDG)
      implicit none 
c.....-----------------------------------------------------------------
c.....gxcp: exchange-correlation potential in ry. also total-energy.
c.....-----------------------------------------------------------------
c     common/cxcf/igl,igh,imj,ibh,ica,icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,
c    &  iex,xlf
c     common/ctrns7/hugeo,huges,hugef,dspr,rdspr,idspr
c.....-----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION AGR,AGRD,AGRU,CEDG,CEDL,G2R,G2RD,G2RU,GGGR,GGGRD,
     +                 GGGRU,GRGRD,GRGRU,GZGR,RO,VCGD,VCGU,VCLD,VCLU,
     +                 VXGD,VXGU,VXLD,VXLU,XCED,XCPTD,XCPTU,XEDG,XEDL,
     +                 ZTA
      INTEGER IDSPR
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,A1,A2,A3,AF,ALC,ALF,ALFC,AP,B,B1,B1F,B1P,B2,
     +                 B2F,B2P,B3,BCR,BETA,BF,BP,BRS,BX,BXD,BXU,BZ41,C,
     +                 C1,C113,C115,C13,C1415,C2,C23,C2915,C2Q23,C3,C32,
     +                 C43,C53,C56,C76,C83,CA,CCF,CCP,CE,CEF,CEP,CF,CGZ,
     +                 CP,CRDC,CRF,CRO,CRP,CRR1,CRR2,D,DACDR,DBDR,DBROD,
     +                 DBROU,DCDR,DD,DECDRF,DECDRP,DF,DFDZ,DLTA,DP,
     +                 DSDFD,DSDFU,DSPR,DSPRS,DVDR1,DVDR2,DVDRD,DVDRU,
     +                 EC,ECF,ECP,ECRS,ECZTA,EF3VI,EXPFAI,F1D,F1U,F2D,
     +                 F2U,F3D,F3U,FAI,FAI2,FD,FDD0,FK,FU,FZ,G,GF,GP,
     +                 GR2,GR2D,GR2U,GZ,GZ2,GZ3,HUGEF,HUGES,PI,Q,
     +                 Q1,Q2,Q3,R,RNC,RO113,RO13,RO2,RO43,RO76,RO83,ROD,
     +                 ROD13,ROD23,ROD3,ROD43,ROD53,ROU,ROU13,ROU23,
     +                 ROU3,ROU43,ROU53,RS,RS2,RS3,SD,SD2,SD3,SD4,SD6,
     +                 SIDFD,SIDFU,SK,SML,SSFC,SU,SU2,SU3,SU4,SU6,TC,TD,
     +                 TKSG,TU,UC,UD,UU,VC,VC13,VC45D,VC45U,VC6,VCCF,
     +                 VCF,VCL1,VCL2,VCP,VXP,VZ,WC,X,X0,X01,X02,X03,
     +                 XEDGD,XEDGU,XEDLD,XEDLU,XF,XL,XL0,XL01,XL02,XL03,
     +                 XL1,XL2,XL3,XLD,XLD1,XLD2,XLD3,XLF,XP,XS,ZT13M,
     +                 ZT13P,ZTA3,ZTA4
      INTEGER IBH,ICA,ICG,IEX,IGD,IGH,IGL,IMJ,IP9,IPG,IVG,IVN,IXLF
C     ..
C     .. External Subroutines ..
      EXTERNAL CORLSD,CPW91,EXCH91
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ACOS,ATAN,EXP,LOG,SQRT
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION FBET,FDEDR,FDFDZ,FFZ,FNCECL,FNCECS,FNCF,FNCVCL,
     +                 FNCVCS,FVNEC,FVQ
C     ..
C     .. Save statement ..
      SAVE GP,GF,B1P,B1F,B2P,B2F,CP,CF,DP,DF,AP,BP,AF,BF,A1,X01,B1,C1,
     +     A2,X02,B2,C2,A3,X03,B3,C3,FDD0,HUGES,HUGEF,DSPR,IGL,
     +     IGH,IMJ,IBH,ICA,ICG,IVN,IPG,IVG,IP9,IGD,IXLF,IEX,XLF
C     ..
C     .. Statement Function definitions ..
      FNCF(X) = (1.d0+X*X*X)*LOG(1.d0+1.d0/X) + X/2.d0 - X*X -
     +          0.333333333d0
      FNCECL(R,G,B1,B2) = G/ (1.d0+B1*SQRT(R)+B2*R)
      FNCVCL(CE,R,B1,B2) = CE* (1.d0+1.16666667d0*B1*SQRT(R)+
     +                     1.33333333d0*B2*R)/ (1.d0+B1*SQRT(R)+B2*R)
      FNCECS(R,A,B,C,D) = A*LOG(R) + B + C*R*LOG(R) + D*R
      FNCVCS(R,A,B,C,D) = A*LOG(R) + (B-A/3.d0) +
     +                    0.666666667d0*C*R*LOG(R) + (2.d0*D-C)*R/3.d0
      FFZ(ZTA) = 1.923661051d0* ((1.d0+ZTA)**1.3333333333d0+
     +           (1.d0-ZTA)**1.3333333333d0-2.d0)
      FDFDZ(ZTA) = 2.564881401d0* ((1.d0+ZTA)**.333333333333d0-
     +             (1.d0-ZTA)**.333333333333d0)
      FVQ(B,C) = SQRT(4.d0*C-B**2)
      FVNEC(A,X,XL,X0,XL0,B,Q) = A* (LOG(X*X/XL)+
     +                           2.d0*B/Q*ATAN(Q/ (2.d0*X+B))-
     +                           B*X0/XL0* (LOG((X-X0)**2/XL)+2.d0* (B+
     +                           2.d0*X0)/Q*ATAN(Q/ (2.d0*X+B))))
      FBET(FDD0,ECF,ECP,ALC) = FDD0* (ECF-ECP)/ALC - 1.d0
      FDEDR(RO,X,A,X0,XL,XL0,XLD,B,Q) = -X/ (6.d0*RO)*A*
     +  ((2.d0*XL-X*XLD)/ (X*XL)-B* (4.d0/ (XLD**2+Q**2)+
     +  X0/XL0* ((2.d0*XL- (X-X0)*XLD)/ ((X-X0)*XL)-4.d0* (B+
     +  2.d0*X0)/ (XLD**2+Q**2))))
C     ..
C     .. Data statements ..
      DATA GP,GF,B1P,B1F,B2P,B2F,CP,CF,DP,DF/-.2846d0,-.1686d0,1.0529d0,
     +     1.3981d0,0.3334d0,0.2611d0,0.0040d0,0.0014d0,-.0232d0,
     +     -.0096d0/
      DATA AP,BP,AF,BF/0.0622d0,-.096d0,0.0311d0,-0.0538d0/
      DATA A1,X01,B1,C1/.0621814d0,-.10498d0,3.72744d0,12.9352d0/
      DATA A2,X02,B2,C2/.0310907d0,-.32500d0,7.06042d0,18.0578d0/
      DATA A3,X03,B3,C3/-.03377373d0,-.0047584d0,1.13107d0,13.0045d0/
      DATA FDD0/1.70992093d0/
      DATA HUGES,HUGEF,DSPR/1.d+6,50.D0,1.d-4/
      DATA IGL,IGH,IMJ,IBH,ICA,ICG,IVN,IPG,IVG,IP9,IGD,IXLF,IEX,
     +     XLF/0,0,0,0,0,0,0,0,0,0,0,0,0,0.00D0/
C     ..
c.....-----------------------------------------------------------------
c.....Perdew-zunger parametrization of Ceperley-Alder. g,a,b,c,d in ry.
c.....for vwn. a1,x01,b1,c1 for para(zta=0), -2 for zta=1, -3 for alfac.
c.....fdd0: fz''(o).
c.....-----------------------------------------------------------------
C
C     PW91   ip9=1  igd=1
C
C     PW86      imj=1  igd=1
C

      IP9 = 1
      IGD = 1
C
C
      PI = ACOS(-1.d0)
      SML = 1.d-12
      IF (ZTA.GT.1.d0-SML) ZTA = 1.d0 - SML
c.....
c     vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),(up,dw)
c     vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),(up,dw)
c     xedl,xedg: exchange energy density (local,grad.exp.) in ry.
c     cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
c.....
      VXLU = 0.0D0
      VCLU = 0.0D0
      VXLD = 0.0D0
      VCLD = 0.0D0
      XEDL = 0.0D0
      CEDL = 0.0D0
      VXGU = 0.0D0
      VCGU = 0.0D0
      VXGD = 0.0D0
      VCGD = 0.0D0
      XEDG = 0.0D0
      CEDG = 0.0D0
c.....
      IF (RO.LT.SML) GO TO 20
c.....
      C13 = 1.D0/3.D0
      C23 = 2.D0/3.D0
      C32 = 3.D0/2.D0
      C43 = 4.D0/3.D0
      C53 = 5.D0/3.D0
      C76 = 7.D0/6.D0
      C113 = 11.D0/3.D0
c.....ca=2.**(-.33333333)
      CA = 0.793700526d0
c.....alf=-3*(3/4*pai)**(1/3).
      ALF = -1.861051473d0
c.....
      RO2 = RO*RO
      RO13 = RO**C13
      RO43 = RO**C43
      RO83 = RO43**2
      RO76 = RO**C76
      RO113 = RO**C113
c.....
      ROU = RO* (1.D0+ZTA)/2.d0
      ROU3 = ROU**3
      ROU13 = ROU**C13
      ROU23 = ROU**C23
      ROU43 = ROU**C43
c.....
      ROD = RO - ROU
      ROD3 = ROD**3
      ROD13 = ROD**C13
      ROD23 = ROD**C23
      ROD43 = ROD**C43
c.....
c     gr2=drr*drr
c     gr2u=drru**2
c     drrd=drr-drru
c     gr2d=drrd**2
c     ddrrd=ddrr-ddrru
c.....
      FZ = FFZ(ZTA)
      DFDZ = FDFDZ(ZTA)
c.....
c.....gz,gz2,gz3: for Wang-Perdew ssf.
      GZ = ((1.D0+ZTA)**C23+ (1.D0-ZTA)**C23)/2.D0
      GZ2 = GZ**2
      GZ3 = GZ**3
c.....
      ZTA3 = ZTA**3
      ZTA4 = ZTA**4
      ZT13P = (1.D0+ZTA)**C13
      ZT13M = (1.D0-ZTA)**C13
c.....
      RS = 0.620350491d0/RO13
      RS2 = RS*RS
      RS3 = RS*RS2
c.....
c.....xedl: exchange-energy-density in ry.
      XEDL = ALF* (ROU43+ROD43)
c.....
c.....exchange-potential, vxp,vxlu,vxld: v-exchange-(para,up,dw).
      VXLU = C43*ALF*ROU13
      VXLD = C43*ALF*ROD13

c.....
      IF (IEX.EQ.1) GO TO 10
c.....

c.....xlfa.
      IF (IXLF.NE.0) THEN
        XLF = 2.D0/3.D0
        VCLU = (XLF*C32-1.D0)*VXLU
        VCLD = (XLF*C32-1.D0)*VXLD
        CEDL = (XLF*C32-1.D0)*XEDL
        GO TO 10
      END IF
c.....
c.....Gunnarson-Lundqvist.(p.r.b13('76),4274,eqs(54)~(56).) or
cc....  g-l but with beta by hedin-lundq.(j.phys.c.4('71),2064))
      IF ((IGL.NE.0) .OR. (IGH.NE.0)) THEN
        XP = RS/11.4d0
        XF = RS/15.9d0
        CEP = -0.0666d0*FNCF(XP)
        CEF = -0.0406d0*FNCF(XF)
        CE = CEP + (CEF-CEP)*FZ
        CEDL = CE*RO
c.....
        IF (IGL.NE.0) BETA = 1.D0 + 0.0545d0*RS*LOG(1.D0+11.4d0/RS)
        IF (IGH.NE.0) BETA = 1.D0 + .03683d0*RS*LOG(1.D0+21.D0/RS)
        DLTA = 1.D0 - 0.036d0*RS + 1.36d0*RS/ (1.D0+10.D0*RS)
c.....
        VXP = C43*ALF*CA*RO13
        VCLU = VXP* (BETA+DLTA/3.D0*ZTA/ (1.D0+0.297d0*ZTA)) - VXLU
        VCLD = VXP* (BETA-DLTA/3.D0*ZTA/ (1.D0-0.297d0*ZTA)) - VXLD
c.....
        GO TO 10
      END IF
c.....
c.....Hedin-von Barth. (j.phys.c.5('72),1629) or moruzzi-janak-williams.
      IF ((IBH.NE.0) .OR. (IMJ.NE.0)) THEN
        IF (IBH.NE.0) THEN
          CRP = 30.D0
          CRF = 75.D0
          CCP = 0.0504d0
          CCF = 0.0254d0
        ELSE IF (IMJ.NE.0) THEN
          CRP = 21.D0
          CRF = 52.916684d0
          CCP = 0.045d0
          CCF = 0.0225d0
c           write(6,*) 'MJW'
        END IF
        XP = RS/CRP
        XF = RS/CRF
        CEP = -CCP*FNCF(XP)
        CEF = -CCF*FNCF(XF)
        CE = CEP + (CEF-CEP)*FZ
        CEDL = CE*RO
c       vclu,vcld: v-correlation-(up,dw). potential.(ry)
        RNC = C43*CA/ (1.D0-CA)* (CEF-CEP)
        VCP = -CCP*LOG(1.D0+CRP/RS)
        BRS = VCP - RNC
        VCLU = RNC*ZT13P + BRS
        VCLD = RNC*ZT13M + BRS
c.....
        GO TO 10
      END IF
c.....
c.....Ceperley-Alder.(paramtrzd by Perdew-zunger.(p.r.23('81),5048)).
      IF (ICA.NE.0) THEN
c.....
        IF (RS.GE.1.d0) THEN
          CEP = FNCECL(RS,GP,B1P,B2P)
          CEF = FNCECL(RS,GF,B1F,B2F)
          VCP = FNCVCL(CEP,RS,B1P,B2P)
          VCF = FNCVCL(CEF,RS,B1F,B2F)
        ELSE
          CEP = FNCECS(RS,AP,BP,CP,DP)
          CEF = FNCECS(RS,AF,BF,CF,DF)
          VCP = FNCVCS(RS,AP,BP,CP,DP)
          VCF = FNCVCS(RS,AF,BF,CF,DF)
        END IF
c.....
        CE = CEP + (CEF-CEP)*FZ
        CEDL = CE*RO
c.....
c.....
        VCL2 = (CEF-CEP)*DFDZ
        VCL1 = VCP + (VCF-VCP)*FZ - VCL2*ZTA
        VCLU = VCL1 + VCL2
        VCLD = VCL1 - VCL2
c.....
        GO TO 10
      END IF
c.....
c.....Ceperley-Alder.with Wang-Perdew spin-scaling-factor.
      IF (ICG.NE.0) THEN
c.....
        IF (RS.GE.1.d0) THEN
          CEP = FNCECL(RS,GP,B1P,B2P)
          VCP = FNCVCL(CEP,RS,B1P,B2P)
        ELSE
          CEP = FNCECS(RS,AP,BP,CP,DP)
          VCP = FNCVCS(RS,AP,BP,CP,DP)
        END IF
c.....
        CE = CEP*GZ3
        CEDL = CE*RO
c.....
        CGZ = CEP*GZ2* (1.D0/ZT13P-1.D0/ZT13M)
        VCL1 = VCP*GZ3 - CGZ*ZTA
        VCLU = VCP*GZ3 + CGZ
        VCLD = VCP*GZ3 - CGZ
c.....
        GO TO 10
      END IF
c.....
c.....Vosko-Wilk-Nusair. Phys.Rev..22,3812,'80.
      IF (IVN.NE.0) THEN
c.....
c.....xl:x-large. xld:d(xl)/dx. xl0:x-large for x=x0.
        XS = SQRT(RS)
        XL1 = XS**2 + B1*XS + C1
        XL2 = XS**2 + B2*XS + C2
        XL3 = XS**2 + B3*XS + C3
        XLD1 = 2.D0*XS + B1
        XLD2 = 2.D0*XS + B2
        XLD3 = 2.D0*XS + B3
        XL01 = X01**2 + B1*X01 + C1
        XL02 = X02**2 + B2*X02 + C2
        XL03 = X03**2 + B3*X03 + C3
        Q1 = FVQ(B1,C1)
        Q2 = FVQ(B2,C2)
        Q3 = FVQ(B3,C3)
        ECP = FVNEC(A1,XS,XL1,X01,XL01,B1,Q1)
        ECF = FVNEC(A2,XS,XL2,X02,XL02,B2,Q2)
        ALC = FVNEC(A3,XS,XL3,X03,XL03,B3,Q3)
        BETA = FBET(FDD0,ECF,ECP,ALC)
        BZ41 = 1.D0 + BETA*ZTA4
c.....
        CE = ECP + ALC*FZ/FDD0*BZ41
        CEDL = CE*RO
c.....
c.....alc: alfac.
c.....decdrp,decdrf: d(ec)/dro-para(zta=0), -(zta=1).
c.....dacdr: d(alc)/dro.
c.....dbdr: d(beta)/dro.
        DECDRP = FDEDR(RO,XS,A1,X01,XL1,XL01,XLD1,B1,Q1)
        DECDRF = FDEDR(RO,XS,A2,X02,XL2,XL02,XLD2,B2,Q2)
        DACDR = FDEDR(RO,XS,A3,X03,XL3,XL03,XLD3,B3,Q3)
c.....
        DBDR = FDD0* ((DECDRF-DECDRP)*ALC- (ECF-ECP)*DACDR)/ALC**2
        VCL1 = CE + RO* (DECDRP+ (DACDR*FZ*BZ41+ALC*FZ*DBDR*ZTA4)/FDD0)
        VCL2 = 2.d0*ALC/ (FDD0*RO)* (DFDZ*BZ41+4.D0*FZ*BETA*ZTA3)
        VCLU = VCL1 + VCL2*ROD
        VCLD = VCL1 + VCL2* (-ROU)
c.....
        GO TO 10
      END IF
c.....
      IF (IP9.EQ.1) THEN
        GO TO 10
      END IF
c.....
   10 CONTINUE

c.....gradient expansion.
c.....
      IF (IGD.LE.0) GO TO 20
c       write(6,*)  '  GGA '
c.....
      GR2 = AGR**2
      GR2U = AGRU**2
      GR2D = AGRD**2

      C56 = 5.D0/6.D0
      C115 = 1.D0/15.D0
      C1415 = 14.D0/15.D0
      C2915 = 29.D0/15.D0
      C2Q23 = 2.D0**C23
      C83 = 8.D0/3.D0
c.....
c.....  dsprs: divergence-suppress-factor.
c       if((log(dspr)+2.*log(agr)-c83*log(ro)).gt.8.0) go to 200
      DSPRS = 1.D0
      IF (IDSPR.EQ.1) DSPRS = EXP(-DSPR*GR2/RO**C83)
c.....
c     agr,agru,agrd: abs(grad(rho)), for all, up, and down.
cc    gr2,gr2u,gr2d: grad(rho_all)**2, grad(rho_up)**2, grad(rho_d)**2.
c     g2r,g2ru,g2rd: laplacian rho_all, _up and _down.
c     gggru,-d: grad(rho)*grad(abs(grad(rho))) for all,up and down.
c     grgru,-d: grad(rho_all)*grad(rhor_up) and for down.

c       g2r=ddrr+2.*drr/rv
c.....
      ROU53 = ROU**C53
c.....
c.....  edrru: d(abs(d(rou)/dr))/dr, edrrd for down.
c       edrru=ddrru
c       if(drru.lt.0.) edrru=-ddrru
c.....
c       agr,agbru,-d: abs(grad(rho)),for rou, rod.
c       gggru,-d: grad(rho)*grad(abs(grad(rho))) for up and down.
c.....  su:at ro=2*rou. 1/(2(3*pai**2)**(1/3))*|grad(rou)|/rou**(4/3).
      SU = 0.128278244d0*AGRU/ROU43
      IF (SU.GT.HUGES) GO TO 20
c       g2ru=ddrru+2.*drru/rv
      TU = .016455307d0*G2RU/ROU53
      UU = 0.002110857d0*GGGRU/ROU3

      IF (IP9.NE.1) THEN

        SU2 = SU*SU
        SU3 = SU*SU2
        SU4 = SU2*SU2
        SU6 = SU2*SU4
c.....
        F1U = 1.d0 + 1.296d0*SU2 + 14.d0*SU4 + .2d0*SU6
        F2U = 2.592d0 + 56.d0*SU2 + 1.2d0*SU4
        F3U = 112.d0*SU + 4.8d0*SU3
c.....
c.....  fu: fgga(su) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
c.....  sidfu: su**(-1)*d(fu)/d(su)).
c.....  dsdfu: d(sidfu)/d(su).
c.....  xedgu; exchange energy density xe at ro=2*rou.(16) of p.w.
cc....      xedgu=ax*rou**(4/3)*(fu-1). ax=2**(4/3)*1.47711(ry).

        FU = F1U**C115
        SIDFU = C115*F1U** (-C1415)*F2U
        DSDFU = C115*F1U** (-C2915)* (-C1415*SU*F2U**2+F1U*F3U)
c.....
        XEDGU = -3.722102942d0* (FU-1.d0)*ROU43
c.....
        VXGU = DSPRS*ALF*ROU13* (C43* (FU-1.D0)-TU*SIDFU-
     +         (UU-C43*SU3)*DSDFU)

      ELSE

        DBROU = ROU*2.D0

        CALL EXCH91(DBROU,SU,UU,TU,XEDLU,XEDGU,VXLU,VXGU)

        XEDL = XEDLU/2.d0

      END IF

c.....
c.....bxu,bxd,bx: grad-coeff. for exchange.

      BXU = XEDGU/GR2U*ROU43
c.....
      ROD53 = ROD**C53
c       edrrd=ddrrd
c       if(drrd.lt.0.) edrrd=-ddrrd

      SD = 0.128278244d0*AGRD/ROD43
      IF (SD.GT.HUGES) GO TO 20

c       g2rd=ddrrd+2.*drrd/rv

      TD = .016455307d0*G2RD/ROD53
      UD = 0.002110857d0*GGGRD/ROD3

      IF (IP9.NE.1) THEN

        SD2 = SD*SD
        SD3 = SD*SD2
        SD4 = SD2*SD2
        SD6 = SD2*SD4
c.....
        F1D = 1.d0 + 1.296d0*SD2 + 14.d0*SD4 + .2d0*SD6
        F2D = 2.592d0 + 56.d0*SD2 + 1.2d0*SD4
        F3D = 112.d0*SD + 4.8d0*SD3
c.....
c.....  fd: fgga(sd) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
c.....  sidfd: sd**(-1)*d(fd)/d(sd)).
c.....  dsdfd: d(sidfd)/d(sd).
c.....  xedgd; exchange energy density xe at ro=2*rod.(16) of p.w.
cc....      xedgd=ax*rod**(4/3)*(fd-1). ax=2**(4/3)*1.47711(ry).

        FD = F1D**C115
        SIDFD = C115*F1D** (-C1415)*F2D
        DSDFD = C115*F1D** (-C2915)* (-C1415*SD*F2D**2+F1D*F3D)
c.....
        XEDGD = -3.722102942d0* (FD-1.d0)*ROD43
c.....
        VXGD = DSPRS*ALF*ROD13* (C43* (FD-1.D0)-TD*SIDFD-
     +         (UD-C43*SD3)*DSDFD)

      ELSE

        DBROD = ROD*2.D0

        CALL EXCH91(DBROD,SD,UD,TD,XEDLD,XEDGD,VXLD,VXGD)

        XEDL = XEDL + XEDLD/2.d0

      END IF

      BXD = XEDGD/GR2D*ROD43
c.....

      XEDG = DSPRS* (XEDGU+XEDGD)/2.d0

      BX = (BXU+BXD)/2.d0

      IF (IEX.EQ.1) GO TO 20

c.....
c.... cro: c(n) of (6),Phys.Rev..b33,8822('86). in ry.
c.... dcdr: d(cro)/d(ro).
c.....0.001625816=1.745*f(=0.11)*cro(rs=0).

      IF (IP9.NE.1) THEN

        CRR1 = .005136d0 + .046532d0*RS + 1.4778d-5*RS2
        CRR2 = 1.D0 + 8.723d0*RS + .472d0*RS2 + .07389d0*RS3
        CRO = .003334d0 + CRR1/CRR2
        DCDR = ((.046532d0+2.9556d-5*RS)*CRR2-
     +         CRR1* (8.723d0+.944d0*RS+.22167d0*RS2))/CRR2/CRR2*
     +         (-RS/RO/3.D0)
c.....
        FAI = 0.001625816d0/CRO*AGR/RO76
        IF (FAI.GT.HUGEF) GO TO 20
        FAI2 = FAI*FAI
        EXPFAI = EXP(-FAI)
c.....
c.....
        IF (IPG.EQ.0) THEN

          DD = 0.707106781d0*SQRT((1.D0+ZTA)**C53+ (1.D0-ZTA)**C53)
c.....    ssfc: spin-scaling-factor for gradient correlation energy.
          SSFC = 1.D0/DD
          CRDC = C56/ (RO113*DD**2)*C2Q23
          VC45U = -CRDC* (ROU23-ROD23)* ((1.D0-FAI)*ROD*GR2-
     +            (2.D0-FAI)*RO*GRGRD)
          VC45D = -CRDC* (ROD23-ROU23)* ((1.D0-FAI)*ROU*GR2-
     +            (2.D0-FAI)*RO*GRGRU)

        ELSE IF (IPG.EQ.1) THEN

          SSFC = GZ
          CRDC = C2Q23/ (3.D0*GZ*RO83)
          VC45U = CRDC* (1.D0/ROU13-1.D0/ROD13)*
     +            ((1.D0-FAI)*ROD*GR2- (2.D0-FAI)*RO*GRGRD)
          VC45D = CRDC* (1.D0/ROD13-1.D0/ROU13)*
     +            ((1.D0-FAI)*ROU*GR2- (2.D0-FAI)*RO*GRGRU)

        ELSE IF (IVG.EQ.1) THEN

          WRITE (6,FMT=
     +      '(/'' non-spher modification not completed for vg'')')
          STOP 30

          IF (IVN.EQ.0) THEN
            WRITE (6,FMT=9000) IVN,IVG
            STOP 16
          END IF
c.....
          DFDZ = FDFDZ(ZTA)
          VZ = (1.D0+ALC/ECP*FZ/FDD0*BZ41)**C13
c.....
          SSFC = VZ
c.....
c.....    dvdru,dvdrd: d(vz)/drou,-d.
          EF3VI = 1.d0/ (ECP*FDD0*3.d0*VZ**2)
          DVDR1 = (DACDR*BZ41-ALC/ECP*DECDRP*BZ41+ALC*DBDR*ZTA4)*FZ*
     +            EF3VI
          DVDR2 = 2.d0* (DFDZ*BZ41+4.d0*FZ*BETA*ZTA3)*ALC/RO2*EF3VI
          DVDRU = DVDR1 + DVDR2*ROD
          DVDRD = DVDR1 - DVDR2*ROU
c.....
          VC45U = ((1.D0-FAI)*GR2*DVDRU-
     +            (2.D0-FAI)*GRGRD* (DVDRU-DVDRD))/ (VZ*RO)
          VC45D = ((1.D0-FAI)*GR2*DVDRD-
     +            (2.D0-FAI)*GRGRU* (DVDRD-DVDRU))/ (VZ*RO)

        END IF

c.....  cedg: correlation-energy-density from grad.expansion.
c.....  bcr: grad-coeff. for correlation.
        BCR = SSFC*EXPFAI*CRO
        CEDG = DSPRS*BCR*GR2/RO43
c.....
c.....  vccf:v-correlation-coeff.
        VCCF = -SSFC*EXPFAI*CRO/RO13
        VC13 = (2.D0-FAI)*G2R/RO - (C43-C113*FAI+C76*FAI2)*GR2/RO2 +
     +         FAI* (FAI-3.D0)*GGGR/AGR/RO
c    &    fai*(fai-3.)*ddrr/ro
        VC6 = -GR2/RO* (FAI2-FAI-1.D0)/CRO*DCDR
c.....
        VCGU = DSPRS*VCCF* (VC13+VC6+VC45U)
c.....
        VCGD = DSPRS*VCCF* (VC13+VC6+VC45D)

      ELSE

c       PW91

        CALL CORLSD(RS,ZTA,EC,VCLU,VCLD,ECRS,ECZTA,ALFC)

        VCLU = VCLU*2.D0
        VCLD = VCLD*2.D0
        CEDL = EC*2.D0*RO

        FK = 1.91915829d0/RS
        SK = SQRT(4.D0*FK/PI)
        TKSG = 2.D0*SK*GZ
        TC = AGR/ (RO*TKSG)
c           gagr: d(ABS(d(ro)/dr))/dr.
c           gagr=ddrr
c           if(drr.lt.0.) gagr=-ddrr
        UC = GGGR/ (RO2*TKSG**3)
c         uc=drr*gagr/(ro2*tksg**3)
        VC = G2R/ (RO*TKSG**2)
        WC = GZGR/ (RO*TKSG**2)
c         wc=drr*dzr/(ro*tksg**2)

        CALL CPW91(FK,SK,GZ,EC,ECRS,ECZTA,RS,ZTA,TC,UC,VC,WC,CEDG,VCGU,
     +             VCGD)

        VCGU = VCGU*2.D0
        VCGD = VCGD*2.D0
        CEDG = CEDG*RO*2.D0*DSPRS

        BCR = CEDG/GR2*RO43

      END IF
c.....
   20 CONTINUE

      XCPTU = VXLU + VCLU + VXGU + VCGU
      XCPTD = VXLD + VCLD + VXGD + VCGD
check
c     ro is small
c
      XCED = 0.0D0
      IF (RO.GT.SML) XCED = (XEDL+CEDL+XEDG+CEDG)/RO

c     write(6,'(/'' vxlu,vxld,vclu,vcld,xedl,cedl ro='',7f11.5)') vxlu,
c    &    vxld,vclu,vcld,xedl,cedl,ro
c       write(6,'(/'' vxgu,vxgd,vcgu,vcgd,xedg,cedg='',6f12.7)') vxgu,
c    &    vxgd,vcgu,vcgd,xedg,cedg

      RETURN
 9000 FORMAT (/,' ivn should be 1 for ivg=1. ivn,ivg=',2i5,/)
      END
