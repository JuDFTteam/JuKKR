      SUBROUTINE MKXCPE2(IR,NP,RV,RHOLM,VXCP,EXCP,YLM,DYLMT1,DYLMF1,
     +                   DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL,IRMD,
     +                   LMPOTD,LMMAX,USE_SOL)
c     ------------------------------------------------------------------
c     Calculation of the exchange-correlation potential.
c     coded by M. Ogura, Apr. 2015, Munich
c     ------------------------------------------------------------------
      implicit none
      integer ijd
      parameter(ijd=434)
c
      double precision rv
      integer ir,irmd,lmmax,lmpotd,np
      double precision ddrrl(irmd,lmpotd),ddrrul(irmd,lmpotd),
     &       drrl(irmd,lmpotd)
     &      ,drrul(irmd,lmpotd),dylmf1(ijd,lmpotd),dylmf2(ijd,lmpotd)
     &      ,dylmt1(ijd,lmpotd),dylmtf(ijd,lmpotd),excp(ijd)
     &      ,rholm(lmpotd,2),vxcp(ijd,2),ylm(ijd,lmpotd)
c
      double precision c,pi,s
      integer ip,ispin,l1,lm,n
      double precision d(2),d1(3,2),d2(5,2),dl(2)
c     use_sol=0 -> PBE, use_sol=1 -> PBEsol
      logical use_sol
      double precision :: um,bet
c
      pi=4d0*atan(1d0)
c   
c Set 'UM' for subroutine 'excpbex' and 'BET' for subroutine 'excpbec' for PBE or PBEsol
      IF (use_sol) THEN
        um=0.123456790123456D0
        bet=0.046D0
      ELSE
        um=0.2195149727645171D0
        bet=0.06672455060314922D0
      END IF
c
c
c     --- surface integration
      do 20 ip=1,np
      do 50 ispin=1,2
      d(ispin)=0d0
      dl(ispin)=0d0
      do 51 n=1,3
 51   d1(n,ispin)=0d0
      do 50 n=1,5
 50   d2(n,ispin)=0d0
      do 53 lm=1,lmmax
      l1=sqrt(dble(lm)-5d-1)
      d(1)=d(1)+rholm(lm,1)*ylm(ip,lm)
      d(2)=d(2)+rholm(lm,2)*ylm(ip,lm)
      dl(1)=dl(1)+dble(l1*(l1+1))*rholm(lm,1)*ylm(ip,lm)
      dl(2)=dl(2)+dble(l1*(l1+1))*rholm(lm,2)*ylm(ip,lm)
      d1(1,2)=d1(1,2)+drrul(ir,lm)*ylm(ip,lm)
      d1(1,1)=d1(1,1)+(drrl(ir,lm)-drrul(ir,lm))*ylm(ip,lm)
      d1(2,1)=d1(2,1)+rholm(lm,1)*dylmt1(ip,lm)
      d1(2,2)=d1(2,2)+rholm(lm,2)*dylmt1(ip,lm)
      d1(3,1)=d1(3,1)+rholm(lm,1)*dylmf1(ip,lm)
      d1(3,2)=d1(3,2)+rholm(lm,2)*dylmf1(ip,lm)
      d2(1,2)=d2(1,2)+ddrrul(ir,lm)*ylm(ip,lm)
      d2(1,1)=d2(1,1)+(ddrrl(ir,lm)-ddrrul(ir,lm))*ylm(ip,lm)
      d2(2,2)=d2(2,2)+drrul(ir,lm)*dylmt1(ip,lm)
      d2(2,1)=d2(2,1)+(drrl(ir,lm)-drrul(ir,lm))*dylmt1(ip,lm)
      d2(3,2)=d2(3,2)+drrul(ir,lm)*dylmf1(ip,lm)
      d2(3,1)=d2(3,1)+(drrl(ir,lm)-drrul(ir,lm))*dylmf1(ip,lm)
      d2(4,1)=d2(4,1)+rholm(lm,1)*dylmtf(ip,lm)
      d2(4,2)=d2(4,2)+rholm(lm,2)*dylmtf(ip,lm)
      d2(5,1)=d2(5,1)+rholm(lm,1)*dylmf2(ip,lm)
 53   d2(5,2)=d2(5,2)+rholm(lm,2)*dylmf2(ip,lm)
      c=ylm(ip,3)*sqrt(4d0*pi/3d0)
      s=sqrt(1d0-c**2)
      do 52 ispin=1,2
      dl(ispin)=dl(ispin)/rv**2
      d1(2,ispin)=d1(2,ispin)/rv
      d2(2,ispin)=d2(2,ispin)/rv
      d2(4,ispin)=d2(4,ispin)/rv**2
      if(s.gt.1d-8)then
      d1(3,ispin)=d1(3,ispin)/rv/s
      d2(3,ispin)=d2(3,ispin)/rv/s
      d2(5,ispin)=d2(5,ispin)/rv**2/s
      call fpexcpbe(d,dl,d1,d2,rv,s,c,vxcp(ip,1),vxcp(ip,2),excp(ip),
     +              um,bet)
      else
      d1(3,ispin)=0d0
      d2(3,ispin)=0d0
      d2(5,ispin)=0d0
      vxcp(ip,1)=0d0 
      vxcp(ip,2)=0d0 
      excp(ip)=0d0 
      endif
 52   continue
 20   continue
c
      return
      end
c
      subroutine fpexcpbe(ro,rol,ro1,ro2,xr,s,c,v1,v2,exc,um,bet)
c----------------------------------------------------------------------
c     driver routine for PBE GGA subroutines.
c     based on excpbe.f in Munich code (version on 20 Dec 2009)
c     coded by M. Ogura, Jun. 2011, Munich
c----------------------------------------------------------------------
      implicit none
c
      double precision c,exc,s,v1,v2,xr
      double precision ro(2),rol(2),ro1(3,2),ro2(5,2)
      double precision um,bet
c
      double precision conf,conrs,d,drv1,drv2,drv2s,drv3,drv4,ec,ex,fk,
     &       g,pi,rs,sk
     &      ,ss,thrd,thrd2,tt,uu,vcdn,vcup,vv,vx,vxcdn,vxcup,ww,x,xd,xu
     &      ,y,z,zet
      integer jsp,llda
c
      pi=4d0*atan(1d0)
      thrd=1d0/3d0
      thrd2=2d0/3d0
      conf=(3d0*pi**2)**thrd
      conrs=(3d0/(4d0*pi))**thrd
      llda=0
      exc=0d0
      vxcup=0d0
      vxcdn=0d0
      if(ro(1) .gt. 1d-12 .and. ro(2) .gt. 1d-12 )then
c
c     ---begin the spin loop for exchange
      if(ro(1) .le. 1d-6 .or. ro(2) .le. 1d-6)llda=1
      do 10 jsp=1,2
      d=2d0*ro(jsp)
      fk=conf*d**thrd
      x=ro1(1,jsp)
      y=ro1(2,jsp)
      z=ro1(3,jsp)
      drv1=sqrt(x**2+y**2+z**2)*2d0
      if(abs(drv1) .lt. 1d-8)then
      drv2=0d0
      else
      drv2s=2d0*y*z*ro2(4,jsp)+(z**2-y**2)*ro2(5,jsp)-c/xr*y*(z**2+y**2)
      drv2=x**2*ro2(1,jsp)+2d0*x*y*ro2(2,jsp)+2d0*x*z*ro2(3,jsp)
     &     -x*y**2/xr-x*z**2/xr-y**2*rol(jsp)
      if(abs(drv2s) .ge. 1d-10)drv2=drv2+drv2s/s
      drv2=drv2/drv1*8d0
      endif
      drv3=ro2(1,jsp)+2d0/xr*x-rol(jsp)
      drv3=drv3*2d0
      ss=drv1/(d*2d0*fk)
      uu=drv2/(d**2*(2d0*fk)**3)
      vv=drv3/(d*(2d0*fk)**2)
      call excpbex(d,ss,uu,vv,ex,vx,llda,um)
      exc=exc+ex*(d/2d0)/(ro(1)+ro(2))
      if(jsp .eq. 1)vxcup=vx
      if(jsp .eq. 2)vxcdn=vx
 10   continue
c
c     ---correlation
      d=ro(1)+ro(2)
      zet=(ro(1)-ro(2))/d
      rs=conrs/d**thrd
      fk=1.91915829d0/rs
      sk=dsqrt(4d0*fk/pi)
      g=((1d0+zet)**thrd2+(1d0-zet)**thrd2)/2d0
      x=ro1(1,1)+ro1(1,2)
      y=ro1(2,1)+ro1(2,2)
      z=ro1(3,1)+ro1(3,2)
      drv1=sqrt(x**2+y**2+z**2)
      if(drv1 .lt. 1d-8)then
      drv2=0d0
      else
      drv2s=2d0*y*z*(ro2(4,1)+ro2(4,2))+(z**2-y**2)*(ro2(5,1)+ro2(5,2))
     &      -c/xr*y*(z**2+y**2)
      drv2=x**2*(ro2(1,1)+ro2(1,2))+2d0*x*y*(ro2(2,1)+ro2(2,2))
     &     +2d0*x*z*(ro2(3,1)+ro2(3,2))-x*y**2/xr-x*z**2/xr
     &     -y**2*(rol(1)+rol(2))
      if(abs(drv2s) .ge. 1d-10)drv2=drv2+drv2s/s
      drv2=drv2/drv1
      endif
      drv3=ro2(1,1)+ro2(1,2)+2d0/xr*x-rol(1)-rol(2)
      drv4=x*(ro1(1,1)-ro1(1,2)-zet*x)+y*(ro1(2,1)-ro1(2,2)-zet*y)
     &     +z*(ro1(3,1)-ro1(3,2)-zet*z)
      tt=drv1/(d*2d0*sk*g)
      uu=drv2/(d**2*(2d0*sk*g)**3)
      vv=drv3/(d*(2d0*sk*g)**2)
      ww=drv4/(d**2*(2d0*sk*g)**2)
      call excpbec(rs,zet,tt,uu,vv,ww,ec,vcup,vcdn,llda,bet)
      exc=exc+ec
      vxcup=vxcup+vcup
      vxcdn=vxcdn+vcdn
c
      endif
c     ---convert from h to ry
      exc=2d0*exc
      xu=2d0*vxcup
      xd=2d0*vxcdn
c
      v1=xu
      v2=xd
      end
c
C*==excpbex.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009
      SUBROUTINE EXCPBEX(RHO,S,U,V,EX,VX,LLDA,UM)
C----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  K Burke's modification of PW91 codes, May 14, 1996
C  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
C----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
C   (for U,V, see PW86(24))
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
C----------------------------------------------------------------------
C References:
C [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
C     Phys. Rev. Lett. 77, 3865 (1996).
C [b] J.P. Perdew and Y. Wang,
C     Phys. Rev. B33, 8800 (1986); B40, 3399 (1989) (E).
C----------------------------------------------------------------------
C Formulas:
C   	e_x[unif]=ax*rho^(4/3)  [LDA]
C ax = -0.75*(3/pi)^(1/3)
C	e_x[PBE]=e_x[unif]*FxPBE(s)
C	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
C uk, ul defined after [a](13)
C----------------------------------------------------------------------
C
C  All input and output is in atomic units!
C
C  Modifications by: E. Engel
C  Last revision:    May 9, 2001
Cengel
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8, PARAMETER :: thrd=1.d0/3.d0
      REAL*8, PARAMETER :: thrd4=4.d0/3.d0
      REAL*8, PARAMETER :: ax=-0.738558766382022405884230032680836D0
      REAL*8, PARAMETER :: uk=0.8040D0
      DOUBLE PRECISION            :: ul
C
C Dummy arguments
C
      DOUBLE PRECISION EX,RHO,S,U,V,VX,UM
      INTEGER LLDA
C
C Local variables
C
      DOUBLE PRECISION EXUNIF,FS,FSS,FXPBE,P0,S2
C
C*** End of declarations rewritten by SPAG
C
C----------------------------------------------------------------------
c    Define UL with via UM and UK
      UL=UM/UK
c
C----------------------------------------------------------------------
C construct LDA exchange energy density
      EXUNIF = AX*RHO**THRD
      IF ( LLDA.EQ.1 ) THEN
         EX = EXUNIF
         VX = EX*THRD4
         RETURN
      END IF
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C construct PBE enhancement factor
      S2 = S*S
      P0 = 1.D0 + UL*S2
      FXPBE = 1D0 + UK - UK/P0
      EX = EXUNIF*FXPBE
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
C  find first and second derivatives of Fx w.r.t s.
C  Fs=(1/s)*d FxPBE/ ds
C  Fss=d Fs/ds
      FS = 2.D0*UK*UL/(P0*P0)
      FSS = -4.D0*UL*S*FS/P0
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C calculate potential from [b](24)
      VX = EXUNIF*(THRD4*FXPBE-(U-THRD4*S2*S)*FSS-V*FS)
      END
C*==excpbec.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009
      SUBROUTINE EXCPBEC(RS,ZETA,T,UU,VV,WW,EC,VCUP,VCDN,LLDA,BET)
Cengel
C  This subroutine evaluates the correlation energy per particle and
C  spin-up and spin-dn correlation potentials within the Perdew-Burke-
C  Ernzerhof GGA. It is a slightly modified version of K. Burke's
C  official PBE subroutine.
C
C  Input:  RS   = WIGNER-SEITZ RADIUS = ( 3 / (4*PI*(DUP+DDN)) )**(1/3)
C          ZETA = RELATIVE SPIN POLARIZATION = (DUP-DDN)/(DUP+DDN)
C          T    = ABS(GRAD D) / ( (2*SK*G) * D )
C          UU   = (GRAD D)*GRAD(ABS(GRAD D)) / ( (2*SK*G)**3 * D**2 )
C          VV   = (LAPLACIAN D) / ( (2*SK*G)**2 * D )
C          WW   = (GRAD D)*(GRAD ZETA) / ( (2*SK*G)**2 * D )
C  where:  FK   = LOCAL FERMI MOMENTUM = (3*PI**2*(DUP+DDN))**(1/3)
C          SK   = LOCAL SCREENING MOMENTUM = (4*FK/PI)**(1/2)
C
C  Output: EC   = correlation energy per particle
C          VCUP = spin-up correlation potential
C          VCDN = spin-dn correlation potential
C
C  All input and output is in atomic units!
C
C References:
C [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
C     Phys. Rev. Lett. 77, 3865 (1996).
C [b] J. P. Perdew, K. Burke, and Y. Wang,
C     Phys. Rev. B54, 16533 (1996).
C [c] J. P. Perdew and Y. Wang,
C     Phys. Rev. B45, 13244 (1992).
C
C
C  Last revision:    May 9, 2001
C  Written by:       K. Burke, May 14, 1996.
C  Modifications by: E. Engel
Cengel
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8, PARAMETER :: thrd=1.d0/3.d0
      REAL*8, PARAMETER :: thrdm=-thrd
      REAL*8, PARAMETER :: thrd2=2.d0*thrd
      REAL*8, PARAMETER :: sixthm=thrdm/2.d0
      REAL*8, PARAMETER :: thrd4=4.d0*thrd
      REAL*8, PARAMETER :: gam=0.5198420997897463295344212145565D0
      REAL*8, PARAMETER :: fzz=8.d0/(9.d0*gam)
      REAL*8, PARAMETER :: gamma=0.03109069086965489503494086371273D0
      REAL*8, PARAMETER :: eta=1.d-12
C
C Dummy arguments
C
      DOUBLE PRECISION EC,RS,T,UU,VCDN,VCUP,VV,WW,ZETA,BET
      INTEGER LLDA
C
C Local variables
C
      DOUBLE PRECISION ALFM,ALFRSM,B,B2,BEC,BG,COMM,ECRS,ECZETA,EP,
     &       EPRS,EU,EURS,F,
     &       FAC,FACT0,FACT1,FACT2,FACT3,FACT5,FZ,G,G3,G4,GZ,H,HB,HBT,
     &       HRS,HRST,HT,HTT,HZ,HZT,PON,PREF,Q4,Q5,Q8,Q9,RSTHRD,RTRS,T2,
     &       T4,T6,Z4,DELT
      EXTERNAL EXCGCOR2
C
C*** End of declarations rewritten by SPAG
C
C thrd*=various multiples of 1/3
C numbers for use in LSD energy spin-interpolation formula, [c](9).
C      GAM= 2^(4/3)-2
C      FZZ=f''(0)= 8/(9*GAM)
C numbers for construction of PBE
C      gamma=(1-log(2))/pi^2
C      bet=coefficient in gradient expansion for correlation, [a](4).
C      eta=small number to stop d phi/ dzeta from blowing up at
C          |zeta|=1.
C----------------------------------------------------------------------
c    Define DELT via BET and GAMMA
      DELT=BET/GAMMA
c
C----------------------------------------------------------------------
C find LSD energy contributions, using [c](10) and Table I[c].
C EU=unpolarized LSD correlation energy
C EURS=dEU/drs
C EP=fully polarized LSD correlation energy
C EPRS=dEP/drs
C ALFM=-spin stiffness, [c](3).
C ALFRSM=-dalpha/drs
C F=spin-scaling factor from [c](9).
C construct ec, using [c](8)
      IF ( RS.LT.3.D5 ) THEN
         RTRS = SQRT(RS)
         CALL EXCGCOR2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     &                 0.49294D0,RTRS,EU,EURS)
         CALL EXCGCOR2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,
     &                 3.3662D0,0.62517D0,RTRS,EP,EPRS)
         CALL EXCGCOR2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,
     &                 0.88026D0,0.49671D0,RTRS,ALFM,ALFRSM)
         Z4 = ZETA**4
         F = ((1.D0+ZETA)**THRD4+(1.D0-ZETA)**THRD4-2.D0)/GAM
         EC = EU*(1.D0-F*Z4) + EP*F*Z4 - ALFM*F*(1.D0-Z4)/FZZ
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C LSD potential from [c](A1)
C ECRS = dEc/drs [c](A2)
C ECZETA=dEc/dzeta [c](A3)
C FZ = dF/dzeta [c](A4)
         ECRS = EURS*(1.D0-F*Z4) + EPRS*F*Z4 - ALFRSM*F*(1.D0-Z4)/FZZ
         FZ = THRD4*((1.D0+ZETA)**THRD-(1.D0-ZETA)**THRD)/GAM
         ECZETA = 4.D0*(ZETA**3)*F*(EP-EU+ALFM/FZZ)
     &            + FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)
         COMM = EC - RS*ECRS/3.D0 - ZETA*ECZETA
         VCUP = COMM + ECZETA
         VCDN = COMM - ECZETA
         IF ( LLDA.EQ.1 ) RETURN
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C PBE correlation energy
C G=phi(zeta), given after [a](3)
C DELT=bet/gamma
C B=A of [a](8)
         G = ((1.D0+ZETA)**THRD2+(1.D0-ZETA)**THRD2)/2.D0
         G3 = G**3
         PON = -EC/(G3*GAMMA)
         B = DELT/(EXP(PON)-1.D0)
         B2 = B*B
         T2 = T*T
         T4 = T2*T2
         Q4 = 1.D0 + B*T2
         Q5 = 1.D0 + B*T2 + B2*T4
         H = G3*(BET/DELT)*LOG(1.D0+DELT*Q4*T2/Q5)
         EC = EC + H
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
         G4 = G3*G
         T6 = T4*T2
         RSTHRD = RS/3.D0
         GZ = (((1.D0+ZETA)**2+ETA)**SIXTHM-((1.D0-ZETA)**2+ETA)
     &        **SIXTHM)/3.D0
         FAC = DELT/B + 1.D0
         BG = -3.D0*B2*EC*FAC/(BET*G4)
         BEC = B2*FAC/(BET*G3)
         Q8 = Q5*Q5 + DELT*Q4*Q5*T2
         Q9 = 1.D0 + 2.D0*B*T2
         HB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
         HRS = -RSTHRD*HB*BEC*ECRS
         FACT0 = 2.D0*DELT - 6.D0*B
         FACT1 = Q5*Q9 + Q4*Q9*Q9
         HBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
         HRST = RSTHRD*T2*HBT*BEC*ECRS
         HZ = 3.D0*GZ*H/G + HB*(BG*GZ+BEC*ECZETA)
         HT = 2.D0*BET*G3*Q9/Q8
         HZT = 3.D0*GZ*HT/G + HBT*(BG*GZ+BEC*ECZETA)
         FACT2 = Q4*Q5 + B*T2*(Q4*Q9+Q5)
         FACT3 = 2.D0*B*Q5*Q9 + DELT*FACT2
         HTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
         COMM = H + HRS + HRST + T2*HT/6.D0 + 7.D0*T2*T*HTT/6.D0
         PREF = HZ - GZ*T2*HT/G
         FACT5 = GZ*(2.D0*HT+T*HTT)/G
         COMM = COMM - PREF*ZETA - UU*HTT - VV*HT - WW*(HZT-FACT5)
         VCUP = VCUP + COMM + PREF
         VCDN = VCDN + COMM - PREF
      ELSE
         VCUP = 0.D0
         VCDN = 0.D0
      END IF
      END
C*==excgcor2.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009
      SUBROUTINE EXCGCOR2(A,A1,B1,B2,B3,B4,RTRS,GG,GGRS)
C----------------------------------------------------------------------
C######################################################################
C----------------------------------------------------------------------
C slimmed down version of GCOR used in PW91 routines, to interpolate
C LSD correlation energy, as given by (10) of
C J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
C K. Burke, May 11, 1996.
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      DOUBLE PRECISION A,A1,B1,B2,B3,B4,GG,GGRS,RTRS
C
C Local variables
C
      DOUBLE PRECISION Q0,Q1,Q2,Q3
C
C*** End of declarations rewritten by SPAG
C
      Q0 = -2.D0*A*(1.D0+A1*RTRS*RTRS)
      Q1 = 2.D0*A*RTRS*(B1+RTRS*(B2+RTRS*(B3+B4*RTRS)))
      Q2 = LOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RTRS+2.D0*B2+RTRS*(3.D0*B3+4.D0*B4*RTRS))
      GGRS = -2.D0*A*A1*Q2 - Q0*Q3/(Q1*(1.D0+Q1))
      END
