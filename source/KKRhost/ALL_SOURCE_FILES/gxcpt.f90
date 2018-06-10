SUBROUTINE gxcpt(idspr,ro,zta,agr,agru,agrd,g2r,g2ru,g2rd,gggr,  &
        gggru,gggrd,grgru,grgrd,gzgr,xcptu,xcptd,xced,  &
        vxlu,vxld,vclu,vcld,xedl,cedl,vxgu,vxgd,vcgu,  &
        vcgd,xedg,cedg)
!.....-----------------------------------------------------------------
!.....gxcp: exchange-correlation potential in ry. also total-energy.
!.....-----------------------------------------------------------------
implicit none
!.. Scalar Arguments ..
DOUBLE PRECISION AGR,AGRD,AGRU,CEDG,CEDL,G2R,G2RD,G2RU,GGGR,GGGRD, &
                 GGGRU,GRGRD,GRGRU,GZGR,RO,VCGD,VCGU,VCLD,VCLU, &
                 VXGD,VXGU,VXLD,VXLU,XCED,XCPTD,XCPTU,XEDG,XEDL, &
                 ZTA
INTEGER IDSPR
!..
!.. Local Scalars ..
DOUBLE PRECISION A,A1,A2,A3,AF,ALC,ALF,ALFC,AP,B,B1,B1F,B1P,B2, &
                 B2F,B2P,B3,BCR,BETA,BF,BP,BRS,BX,BXD,BXU,BZ41,C, &
                 C1,C113,C115,C13,C1415,C2,C23,C2915,C2Q23,C3,C32, &
                 C43,C53,C56,C76,C83,CA,CCF,CCP,CE,CEF,CEP,CF,CGZ, &
                 CP,CRDC,CRF,CRO,CRP,CRR1,CRR2,D,DACDR,DBDR,DBROD, &
                 DBROU,DCDR,DD,DECDRF,DECDRP,DF,DFDZ,DLTA,DP, &
                 DSDFD,DSDFU,DSPR,DSPRS,DVDR1,DVDR2,DVDRD,DVDRU, &
                 EC,ECF,ECP,ECRS,ECZTA,EF3VI,EXPFAI,F1D,F1U,F2D, &
                 F2U,F3D,F3U,FAI,FAI2,FD,FDD0,FK,FU,FZ,G,GF,GP, &
                 GR2,GR2D,GR2U,GZ,GZ2,GZ3,HUGEF,HUGES,PI,Q, &
                 Q1,Q2,Q3,R,RNC,RO113,RO13,RO2,RO43,RO76,RO83,ROD, &
                 ROD13,ROD23,ROD3,ROD43,ROD53,ROU,ROU13,ROU23, &
                 ROU3,ROU43,ROU53,RS,RS2,RS3,SD,SD2,SD3,SD4,SD6, &
                 SIDFD,SIDFU,SK,SML,SSFC,SU,SU2,SU3,SU4,SU6,TC,TD, &
                 TKSG,TU,UC,UD,UU,VC,VC13,VC45D,VC45U,VC6,VCCF, &
                 VCF,VCL1,VCL2,VCP,VXP,VZ,WC,X,X0,X01,X02,X03, &
                 XEDGD,XEDGU,XEDLD,XEDLU,XF,XL,XL0,XL01,XL02,XL03, &
                 XL1,XL2,XL3,XLD,XLD1,XLD2,XLD3,XLF,XP,XS,ZT13M, &
                 ZT13P,ZTA3,ZTA4
INTEGER IBH,ICA,ICG,IEX,IGD,IGH,IGL,IMJ,IP9,IPG,IVG,IVN,IXLF
!..
!.. External Subroutines ..
EXTERNAL CORLSD,CPW91,EXCH91
!..
!.. Intrinsic Functions ..
INTRINSIC ACOS,ATAN,EXP,LOG,SQRT
!..
!.. Statement Functions ..
DOUBLE PRECISION FBET,FDEDR,FDFDZ,FFZ,FNCECL,FNCECS,FNCF,FNCVCL, &
                 FNCVCS,FVNEC,FVQ
!..
!.. Save statement ..
SAVE GP,GF,B1P,B1F,B2P,B2F,CP,CF,DP,DF,AP,BP,AF,BF,A1,X01,B1,C1, &
     A2,X02,B2,C2,A3,X03,B3,C3,FDD0,HUGES,HUGEF,DSPR,IGL, &
     IGH,IMJ,IBH,ICA,ICG,IVN,IPG,IVG,IP9,IGD,IXLF,IEX,XLF
!..
!.. Statement Function definitions ..
FNCF(X) = (1.d0+X*X*X)*LOG(1.d0+1.d0/X) + X/2.d0 - X*X - &
          0.333333333d0
FNCECL(R,G,B1,B2) = G/ (1.d0+B1*SQRT(R)+B2*R)
FNCVCL(CE,R,B1,B2) = CE* (1.d0+1.16666667d0*B1*SQRT(R)+ &
                     1.33333333d0*B2*R)/ (1.d0+B1*SQRT(R)+B2*R)
FNCECS(R,A,B,C,D) = A*LOG(R) + B + C*R*LOG(R) + D*R
FNCVCS(R,A,B,C,D) = A*LOG(R) + (B-A/3.d0) + &
                    0.666666667d0*C*R*LOG(R) + (2.d0*D-C)*R/3.d0
FFZ(ZTA) = 1.923661051d0* ((1.d0+ZTA)**1.3333333333d0+ &
           (1.d0-ZTA)**1.3333333333d0-2.d0)
FDFDZ(ZTA) = 2.564881401d0* ((1.d0+ZTA)**.333333333333d0- &
             (1.d0-ZTA)**.333333333333d0)
FVQ(B,C) = SQRT(4.d0*C-B**2)
FVNEC(A,X,XL,X0,XL0,B,Q) = A* (LOG(X*X/XL)+ &
                           2.d0*B/Q*ATAN(Q/ (2.d0*X+B))- &
                           B*X0/XL0* (LOG((X-X0)**2/XL)+2.d0* (B+ &
                           2.d0*X0)/Q*ATAN(Q/ (2.d0*X+B))))
FBET(FDD0,ECF,ECP,ALC) = FDD0* (ECF-ECP)/ALC - 1.d0
FDEDR(RO,X,A,X0,XL,XL0,XLD,B,Q) = -X/ (6.d0*RO)*A* &
  ((2.d0*XL-X*XLD)/ (X*XL)-B* (4.d0/ (XLD**2+Q**2)+ &
  X0/XL0* ((2.d0*XL- (X-X0)*XLD)/ ((X-X0)*XL)-4.d0* (B+ &
  2.d0*X0)/ (XLD**2+Q**2))))
!..
!.. Data statements ..
DATA GP,GF,B1P,B1F,B2P,B2F,CP,CF,DP,DF/-.2846d0,-.1686d0,1.0529d0, &
     1.3981d0,0.3334d0,0.2611d0,0.0040d0,0.0014d0,-.0232d0, &
     -.0096d0/
DATA AP,BP,AF,BF/0.0622d0,-.096d0,0.0311d0,-0.0538d0/
DATA A1,X01,B1,C1/.0621814d0,-.10498d0,3.72744d0,12.9352d0/
DATA A2,X02,B2,C2/.0310907d0,-.32500d0,7.06042d0,18.0578d0/
DATA A3,X03,B3,C3/-.03377373d0,-.0047584d0,1.13107d0,13.0045d0/
DATA FDD0/1.70992093d0/
DATA HUGES,HUGEF,DSPR/1.d+6,50.D0,1.d-4/
DATA IGL,IGH,IMJ,IBH,ICA,ICG,IVN,IPG,IVG,IP9,IGD,IXLF,IEX, &
     XLF/0,0,0,0,0,0,0,0,0,0,0,0,0,0.00D0/
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


pi = ACOS(-1.d0)
sml = 1.d-12
IF (zta > 1.d0-sml) zta = 1.d0 - sml
!.....
!     vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),(up,dw)
!     vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),(up,dw)
!     xedl,xedg: exchange energy density (local,grad.exp.) in ry.
!     cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
!.....
vxlu = 0.0D0
vclu = 0.0D0
vxld = 0.0D0
vcld = 0.0D0
xedl = 0.0D0
cedl = 0.0D0
vxgu = 0.0D0
vcgu = 0.0D0
vxgd = 0.0D0
vcgd = 0.0D0
xedg = 0.0D0
cedg = 0.0D0
!.....
IF (ro < sml) GO TO 20
!.....
c13 = 1.d0/3.d0
c23 = 2.d0/3.d0
c32 = 3.d0/2.d0
c43 = 4.d0/3.d0
c53 = 5.d0/3.d0
c76 = 7.d0/6.d0
c113 = 11.d0/3.d0
!.....ca=2.**(-.33333333)
ca = 0.793700526D0
!.....alf=-3*(3/4*pai)**(1/3).
alf = -1.861051473D0
!.....
ro2 = ro*ro
ro13 = ro**c13
ro43 = ro**c43
ro83 = ro43**2
ro76 = ro**c76
ro113 = ro**c113
!.....
rou = ro* (1.d0+zta)/2.d0
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
gz = ((1.d0+zta)**c23+ (1.d0-zta)**c23)/2.d0
gz2 = gz**2
gz3 = gz**3
!.....
zta3 = zta**3
zta4 = zta**4
zt13p = (1.d0+zta)**c13
zt13m = (1.d0-zta)**c13
!.....
rs = 0.620350491D0/ro13
rs2 = rs*rs
rs3 = rs*rs2
!.....
!.....xedl: exchange-energy-density in ry.
xedl = alf* (rou43+rod43)
!.....
!.....exchange-potential, vxp,vxlu,vxld: v-exchange-(para,up,dw).
vxlu = c43*alf*rou13
vxld = c43*alf*rod13

!.....
IF (iex == 1) GO TO 10
!.....

!.....xlfa.
IF (ixlf /= 0) THEN
  xlf = 2.d0/3.d0
  vclu = (xlf*c32-1.d0)*vxlu
  vcld = (xlf*c32-1.d0)*vxld
  cedl = (xlf*c32-1.d0)*xedl
  GO TO 10
END IF
!.....
!.....Gunnarson-Lundqvist.(p.r.b13('76),4274,eqs(54)~(56).) or
!c....  g-l but with beta by hedin-lundq.(j.phys.c.4('71),2064))
IF ((igl /= 0) .OR. (igh /= 0)) THEN
  xp = rs/11.4D0
  xf = rs/15.9D0
  cep = -0.0666D0*fncf(xp)
  cef = -0.0406D0*fncf(xf)
  ce = cep + (cef-cep)*fz
  cedl = ce*ro
!.....
  IF (igl /= 0) beta = 1.d0 + 0.0545D0*rs*LOG(1.d0+11.4D0/rs)
  IF (igh /= 0) beta = 1.d0 + .03683D0*rs*LOG(1.d0+21.d0/rs)
  dlta = 1.d0 - 0.036D0*rs + 1.36D0*rs/ (1.d0+10.d0*rs)
!.....
  vxp = c43*alf*ca*ro13
  vclu = vxp* (beta+dlta/3.d0*zta/ (1.d0+0.297D0*zta)) - vxlu
  vcld = vxp* (beta-dlta/3.d0*zta/ (1.d0-0.297D0*zta)) - vxld
!.....
  GO TO 10
END IF
!.....
!.....Hedin-von Barth. (j.phys.c.5('72),1629) or moruzzi-janak-williams.
IF ((ibh /= 0) .OR. (imj /= 0)) THEN
  IF (ibh /= 0) THEN
    crp = 30.d0
    crf = 75.d0
    ccp = 0.0504D0
    ccf = 0.0254D0
  ELSE IF (imj /= 0) THEN
    crp = 21.d0
    crf = 52.916684D0
    ccp = 0.045D0
    ccf = 0.0225D0
!           write(6,*) 'MJW'
  END IF
  xp = rs/crp
  xf = rs/crf
  cep = -ccp*fncf(xp)
  cef = -ccf*fncf(xf)
  ce = cep + (cef-cep)*fz
  cedl = ce*ro
!       vclu,vcld: v-correlation-(up,dw). potential.(ry)
  rnc = c43*ca/ (1.d0-ca)* (cef-cep)
  vcp = -ccp*LOG(1.d0+crp/rs)
  brs = vcp - rnc
  vclu = rnc*zt13p + brs
  vcld = rnc*zt13m + brs
!.....
  GO TO 10
END IF
!.....
!.....Ceperley-Alder.(paramtrzd by Perdew-zunger.(p.r.23('81),5048)).
IF (ica /= 0) THEN
!.....
  IF (rs >= 1.d0) THEN
    cep = fncecl(rs,gp,b1p,b2p)
    cef = fncecl(rs,gf,b1f,b2f)
    vcp = fncvcl(cep,rs,b1p,b2p)
    vcf = fncvcl(cef,rs,b1f,b2f)
  ELSE
    cep = fncecs(rs,ap,bp,cp,dp)
    cef = fncecs(rs,af,bf,cf,df)
    vcp = fncvcs(rs,ap,bp,cp,dp)
    vcf = fncvcs(rs,af,bf,cf,df)
  END IF
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
  GO TO 10
END IF
!.....
!.....Ceperley-Alder.with Wang-Perdew spin-scaling-factor.
IF (icg /= 0) THEN
!.....
  IF (rs >= 1.d0) THEN
    cep = fncecl(rs,gp,b1p,b2p)
    vcp = fncvcl(cep,rs,b1p,b2p)
  ELSE
    cep = fncecs(rs,ap,bp,cp,dp)
    vcp = fncvcs(rs,ap,bp,cp,dp)
  END IF
!.....
  ce = cep*gz3
  cedl = ce*ro
!.....
  cgz = cep*gz2* (1.d0/zt13p-1.d0/zt13m)
  vcl1 = vcp*gz3 - cgz*zta
  vclu = vcp*gz3 + cgz
  vcld = vcp*gz3 - cgz
!.....
  GO TO 10
END IF
!.....
!.....Vosko-Wilk-Nusair. Phys.Rev..22,3812,'80.
IF (ivn /= 0) THEN
!.....
!.....xl:x-large. xld:d(xl)/dx. xl0:x-large for x=x0.
  xs = SQRT(rs)
  xl1 = xs**2 + b1*xs + c1
  xl2 = xs**2 + b2*xs + c2
  xl3 = xs**2 + b3*xs + c3
  xld1 = 2.d0*xs + b1
  xld2 = 2.d0*xs + b2
  xld3 = 2.d0*xs + b3
  xl01 = x01**2 + b1*x01 + c1
  xl02 = x02**2 + b2*x02 + c2
  xl03 = x03**2 + b3*x03 + c3
  q1 = fvq(b1,c1)
  q2 = fvq(b2,c2)
  q3 = fvq(b3,c3)
  ecp = fvnec(a1,xs,xl1,x01,xl01,b1,q1)
  ecf = fvnec(a2,xs,xl2,x02,xl02,b2,q2)
  alc = fvnec(a3,xs,xl3,x03,xl03,b3,q3)
  beta = fbet(fdd0,ecf,ecp,alc)
  bz41 = 1.d0 + beta*zta4
!.....
  ce = ecp + alc*fz/fdd0*bz41
  cedl = ce*ro
!.....
!.....alc: alfac.
!.....decdrp,decdrf: d(ec)/dro-para(zta=0), -(zta=1).
!.....dacdr: d(alc)/dro.
!.....dbdr: d(beta)/dro.
  decdrp = fdedr(ro,xs,a1,x01,xl1,xl01,xld1,b1,q1)
  decdrf = fdedr(ro,xs,a2,x02,xl2,xl02,xld2,b2,q2)
  dacdr = fdedr(ro,xs,a3,x03,xl3,xl03,xld3,b3,q3)
!.....
  dbdr = fdd0* ((decdrf-decdrp)*alc- (ecf-ecp)*dacdr)/alc**2
  vcl1 = ce + ro* (decdrp+ (dacdr*fz*bz41+alc*fz*dbdr*zta4)/fdd0)
  vcl2 = 2.d0*alc/ (fdd0*ro)* (dfdz*bz41+4.d0*fz*beta*zta3)
  vclu = vcl1 + vcl2*rod
  vcld = vcl1 + vcl2* (-rou)
!.....
  GO TO 10
END IF
!.....
IF (ip9 == 1) THEN
  GO TO 10
END IF
!.....
10 CONTINUE

!.....gradient expansion.
!.....
IF (igd <= 0) GO TO 20
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
IF (idspr == 1) dsprs = EXP(-dspr*gr2/ro**c83)
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
su = 0.128278244D0*agru/rou43
IF (su > huges) GO TO 20
!       g2ru=ddrru+2.*drru/rv
tu = .016455307D0*g2ru/rou53
uu = 0.002110857D0*gggru/rou3

IF (ip9 /= 1) THEN
  
  su2 = su*su
  su3 = su*su2
  su4 = su2*su2
  su6 = su2*su4
!.....
  f1u = 1.d0 + 1.296D0*su2 + 14.d0*su4 + .2D0*su6
  f2u = 2.592D0 + 56.d0*su2 + 1.2D0*su4
  f3u = 112.d0*su + 4.8D0*su3
!.....
!.....  fu: fgga(su) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
!.....  sidfu: su**(-1)*d(fu)/d(su)).
!.....  dsdfu: d(sidfu)/d(su).
!.....  xedgu; exchange energy density xe at ro=2*rou.(16) of p.w.
!c....      xedgu=ax*rou**(4/3)*(fu-1). ax=2**(4/3)*1.47711(ry).
  
  fu = f1u**c115
  sidfu = c115*f1u** (-c1415)*f2u
  dsdfu = c115*f1u** (-c2915)* (-c1415*su*f2u**2+f1u*f3u)
!.....
  xedgu = -3.722102942D0* (fu-1.d0)*rou43
!.....
  vxgu = dsprs*alf*rou13* (c43* (fu-1.d0)-tu*sidfu- (uu-c43*su3)*dsdfu)
  
ELSE
  
  dbrou = rou*2.d0
  
  CALL exch91(dbrou,su,uu,tu,xedlu,xedgu,vxlu,vxgu)
  
  xedl = xedlu/2.d0
  
END IF

!.....
!.....bxu,bxd,bx: grad-coeff. for exchange.

bxu = xedgu/gr2u*rou43
!.....
rod53 = rod**c53
!       edrrd=ddrrd
!       if(drrd.lt.0.) edrrd=-ddrrd

sd = 0.128278244D0*agrd/rod43
IF (sd > huges) GO TO 20

!       g2rd=ddrrd+2.*drrd/rv

td = .016455307D0*g2rd/rod53
ud = 0.002110857D0*gggrd/rod3

IF (ip9 /= 1) THEN
  
  sd2 = sd*sd
  sd3 = sd*sd2
  sd4 = sd2*sd2
  sd6 = sd2*sd4
!.....
  f1d = 1.d0 + 1.296D0*sd2 + 14.d0*sd4 + .2D0*sd6
  f2d = 2.592D0 + 56.d0*sd2 + 1.2D0*sd4
  f3d = 112.d0*sd + 4.8D0*sd3
!.....
!.....  fd: fgga(sd) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
!.....  sidfd: sd**(-1)*d(fd)/d(sd)).
!.....  dsdfd: d(sidfd)/d(sd).
!.....  xedgd; exchange energy density xe at ro=2*rod.(16) of p.w.
!c....      xedgd=ax*rod**(4/3)*(fd-1). ax=2**(4/3)*1.47711(ry).
  
  fd = f1d**c115
  sidfd = c115*f1d** (-c1415)*f2d
  dsdfd = c115*f1d** (-c2915)* (-c1415*sd*f2d**2+f1d*f3d)
!.....
  xedgd = -3.722102942D0* (fd-1.d0)*rod43
!.....
  vxgd = dsprs*alf*rod13* (c43* (fd-1.d0)-td*sidfd- (ud-c43*sd3)*dsdfd)
  
ELSE
  
  dbrod = rod*2.d0
  
  CALL exch91(dbrod,sd,ud,td,xedld,xedgd,vxld,vxgd)
  
  xedl = xedl + xedld/2.d0
  
END IF

bxd = xedgd/gr2d*rod43
!.....

xedg = dsprs* (xedgu+xedgd)/2.d0

bx = (bxu+bxd)/2.d0

IF (iex == 1) GO TO 20

!.....
!.... cro: c(n) of (6),Phys.Rev..b33,8822('86). in ry.
!.... dcdr: d(cro)/d(ro).
!.....0.001625816=1.745*f(=0.11)*cro(rs=0).

IF (ip9 /= 1) THEN
  
  crr1 = .005136D0 + .046532D0*rs + 1.4778D-5*rs2
  crr2 = 1.d0 + 8.723D0*rs + .472D0*rs2 + .07389D0*rs3
  cro = .003334D0 + crr1/crr2
  dcdr = ((.046532D0+2.9556D-5*rs)*crr2-  &
      crr1* (8.723D0+.944D0*rs+.22167D0*rs2))/crr2/crr2* (-rs/ro/3.d0)
!.....
  fai = 0.001625816D0/cro*agr/ro76
  IF (fai > hugef) GO TO 20
  fai2 = fai*fai
  expfai = EXP(-fai)
!.....
!.....
  IF (ipg == 0) THEN
    
    dd = 0.707106781D0*SQRT((1.d0+zta)**c53+ (1.d0-zta)**c53)
!.....    ssfc: spin-scaling-factor for gradient correlation energy.
    ssfc = 1.d0/dd
    crdc = c56/ (ro113*dd**2)*c2q23
    vc45u = -crdc* (rou23-rod23)* ((1.d0-fai)*rod*gr2- (2.d0-fai)*ro*grgrd)
    vc45d = -crdc* (rod23-rou23)* ((1.d0-fai)*rou*gr2- (2.d0-fai)*ro*grgru)
    
  ELSE IF (ipg == 1) THEN
    
    ssfc = gz
    crdc = c2q23/ (3.d0*gz*ro83)
    vc45u = crdc* (1.d0/rou13-1.d0/rod13)*  &
        ((1.d0-fai)*rod*gr2- (2.d0-fai)*ro*grgrd)
    vc45d = crdc* (1.d0/rod13-1.d0/rou13)*  &
        ((1.d0-fai)*rou*gr2- (2.d0-fai)*ro*grgru)
    
  ELSE IF (ivg == 1) THEN
    
    WRITE (6,FMT= '(/'' NON-SPHER MODIFICATION NOT COMPLETED FOR VG'')')
    STOP 30
    
    IF (ivn == 0) THEN
      WRITE (6,FMT=9000) ivn,ivg
      STOP 16
    END IF
!.....
    dfdz = fdfdz(zta)
    vz = (1.d0+alc/ecp*fz/fdd0*bz41)**c13
!.....
    ssfc = vz
!.....
!.....    dvdru,dvdrd: d(vz)/drou,-d.
    ef3vi = 1.d0/ (ecp*fdd0*3.d0*vz**2)
    dvdr1 = (dacdr*bz41-alc/ecp*decdrp*bz41+alc*dbdr*zta4)*fz* ef3vi
    dvdr2 = 2.d0* (dfdz*bz41+4.d0*fz*beta*zta3)*alc/ro2*ef3vi
    dvdru = dvdr1 + dvdr2*rod
    dvdrd = dvdr1 - dvdr2*rou
!.....
    vc45u = ((1.d0-fai)*gr2*dvdru- (2.d0-fai)*grgrd* (dvdru-dvdrd))/ (vz*ro)
    vc45d = ((1.d0-fai)*gr2*dvdrd- (2.d0-fai)*grgru* (dvdrd-dvdru))/ (vz*ro)
    
  END IF
  
!.....  cedg: correlation-energy-density from grad.expansion.
!.....  bcr: grad-coeff. for correlation.
  bcr = ssfc*expfai*cro
  cedg = dsprs*bcr*gr2/ro43
!.....
!.....  vccf:v-correlation-coeff.
  vccf = -ssfc*expfai*cro/ro13
  vc13 = (2.d0-fai)*g2r/ro - (c43-c113*fai+c76*fai2)*gr2/ro2 +  &
      fai* (fai-3.d0)*gggr/agr/ro
!    &    fai*(fai-3.)*ddrr/ro
  vc6 = -gr2/ro* (fai2-fai-1.d0)/cro*dcdr
!.....
  vcgu = dsprs*vccf* (vc13+vc6+vc45u)
!.....
  vcgd = dsprs*vccf* (vc13+vc6+vc45d)
  
ELSE
  
!       PW91
  
  CALL corlsd(rs,zta,ec,vclu,vcld,ecrs,eczta,alfc)
  
  vclu = vclu*2.d0
  vcld = vcld*2.d0
  cedl = ec*2.d0*ro
  
  fk = 1.91915829D0/rs
  sk = SQRT(4.d0*fk/pi)
  tksg = 2.d0*sk*gz
  tc = agr/ (ro*tksg)
!           gagr: d(ABS(d(ro)/dr))/dr.
!           gagr=ddrr
!           if(drr.lt.0.) gagr=-ddrr
  uc = gggr/ (ro2*tksg**3)
!         uc=drr*gagr/(ro2*tksg**3)
  vc = g2r/ (ro*tksg**2)
  wc = gzgr/ (ro*tksg**2)
!         wc=drr*dzr/(ro*tksg**2)
  
  CALL cpw91(fk,sk,gz,ec,ecrs,eczta,rs,zta,tc,uc,vc,wc,cedg,vcgu, vcgd)
  
  vcgu = vcgu*2.d0
  vcgd = vcgd*2.d0
  cedg = cedg*ro*2.d0*dsprs
  
  bcr = cedg/gr2*ro43
  
END IF
!.....
20 CONTINUE

xcptu = vxlu + vclu + vxgu + vcgu
xcptd = vxld + vcld + vxgd + vcgd
!heck
!     ro is small

xced = 0.0D0
IF (ro > sml) xced = (xedl+cedl+xedg+cedg)/ro

!     write(6,'(/'' vxlu,vxld,vclu,vcld,xedl,cedl ro='',7f11.5)') vxlu,
!    &    vxld,vclu,vcld,xedl,cedl,ro
!       write(6,'(/'' vxgu,vxgd,vcgu,vcgd,xedg,cedg='',6f12.7)') vxgu,
!    &    vxgd,vcgu,vcgd,xedg,cedg

RETURN
9000 FORMAT (/,' ivn should be 1 for ivg=1. ivn,ivg=',2I5,/)
END SUBROUTINE gxcpt
