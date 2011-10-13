      SUBROUTINE CPW91(FK,SK,GZ,EC,ECRS,ECZTA,RS,ZTA,T,UU,VV,WW,H,DVCUP,
     +                 DVCDN)
c.....-----------------------------------------------------------------
c     gga91 correlation
c.....-----------------------------------------------------------------
c     input
c           rs: seitz radius
c         zta: relative spin polarization
c            t: abs(grad d)/(d*2.*ks*gz)
c           uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*gz)**3)
c           vv: (laplacian d)/(d * (2*ks*gz)**2)
c           ww: (grad d)*(gradzta)/(d * (2*ks*gz)**2
c      output
c                  h: nonlocal part of correlation energy per electron
c          dvcup,-dn: nonlocal parts of correlation potentials.

c      with ks=sqrt(4*kf/pai), gz=[(1+zta)**(2/3)+(1-zta)**(2/3)]/2, &
c           kf=cbrt(3*pai**2*d).
c.....-----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION DVCDN,DVCUP,EC,ECRS,ECZTA,FK,GZ,H,RS,SK,T,UU,VV,
     +                 WW,ZTA
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A4,ALF,ARGMX,B,B2,B2FAC,BEC,BET,BG,C1,C2,C3,C4,
     +                 C5,C6,CC,CC0,CCRS,COEFF,COMM,CX,DELT,FACT0,FACT1,
     +                 FACT2,FACT3,FACT4,FACT5,GM,GZ3,GZ4,H0,H0B,H0BT,
     +                 H0RS,H0RST,H0T,H0TT,H0Z,H0ZT,H1,H1RS,H1RST,H1T,
     +                 H1TT,H1Z,H1ZT,HRS,HRST,HT,HTT,HZ,HZT,PON,PREF,Q4,
     +                 Q5,Q6,Q7,Q8,Q9,R0,R1,R2,R3,R4,RS2,RS3,RSTHRD,T2,
     +                 T4,T6,THRD2,THRDM,XNU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,LOG
C     ..
C     .. Save statement ..
      SAVE XNU,CC0,CX,ALF,C1,C2,C3,C4,C5,C6,A4,THRDM,THRD2
C     ..
C     .. Data statements ..
c.....-----------------------------------------------------------------
      DATA XNU,CC0,CX,ALF/15.75592d0,0.004235d0,-0.001667212d0,0.09d0/
      DATA C1,C2,C3,C4/0.002568d0,0.023266d0,7.389d-6,8.723d0/
      DATA C5,C6,A4/0.472d0,7.389d-2,100.d0/
      DATA THRDM,THRD2/-0.333333333333d0,0.666666666667d0/
C     ..
c.....-----------------------------------------------------------------
      ARGMX = 174.0D0
      BET = XNU*CC0
      DELT = 2.d0*ALF/BET
      GZ3 = GZ**3
      GZ4 = GZ3*GZ
      PON = -DELT*EC/ (GZ3*BET)
      IF (PON.GT.ARGMX) THEN
        B = 0.D0
      ELSE
        B = DELT/ (EXP(PON)-1.d0)
      END IF
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.d0 + B*T2
      Q5 = 1.d0 + B*T2 + B2*T4
      Q6 = C1 + C2*RS + C3*RS2
      Q7 = 1.d0 + C4*RS + C5*RS2 + C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*GZ4
      COEFF = CC - CC0 - 3.d0*CX/7.d0
      R2 = XNU*COEFF*GZ3
      R3 = EXP(-R1*T2)
      H0 = GZ3* (BET/DELT)*LOG(1.d0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0 + H1
c  local correlation option:
c     h = 0.0d0

c  energy done. now the potential:

      CCRS = (C2+2.D0*C3*RS)/Q7 - Q6* (C4+2.D0*C5*RS+3.D0*C6*RS2)/Q7**2
      RSTHRD = RS/3.d0
      R4 = RSTHRD*CCRS/COEFF
      GM = ((1.d0+ZTA)**THRDM- (1.d0-ZTA)**THRDM)/3.d0
      IF (PON.GT.ARGMX) THEN
        B2FAC = 0.D0
      ELSE
        B2FAC = B2* (DELT/B+1.d0)
      END IF
      BG = -3.d0*EC*B2FAC/ (BET*GZ4)
      BEC = B2FAC/ (BET*GZ3)
      Q8 = Q5*Q5 + DELT*Q4*Q5*T2
      Q9 = 1.d0 + 2.d0*B*T2
      H0B = -BET*GZ3*B*T6* (2.d0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.d0*DELT - 6.d0*B
      FACT1 = Q5*Q9 + Q4*Q9*Q9
      H0BT = 2.d0*BET*GZ3*T4* ((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.d0*GM*H0/GZ + H0B* (BG*GM+BEC*ECZTA)
      H0T = 2.D0*BET*GZ3*Q9/Q8
      H0ZT = 3.d0*GM*H0T/GZ + H0BT* (BG*GM+BEC*ECZTA)
      FACT2 = Q4*Q5 + B*T2* (Q4*Q9+Q5)
      FACT3 = 2.d0*B*Q5*Q9 + DELT*FACT2
      H0TT = 4.d0*BET*GZ3*T* (2.d0*B/Q8- (Q9* (FACT3/Q8))/Q8)
      H1RS = R3*R2*T2* (-R4+R1*T2/3.d0)
      FACT4 = 2.d0 - R1*T2
      H1RST = R3*R2*T2* (2.d0*R4* (1.d0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GM*R3*R2*T2* (3.d0-4.d0*R1*T2)/GZ
      H1T = 2.d0*R3*R2* (1.d0-R1*T2)
      H1ZT = 2.d0*GM*R3*R2* (3.d0-11.d0*R1*T2+4.d0*R1*R1*T4)/GZ
      H1TT = 4.d0*R3*R2*R1*T* (-2.d0+R1*T2)
      HRS = H0RS + H1RS
      HRST = H0RST + H1RST
      HT = H0T + H1T
      HTT = H0TT + H1TT
      HZ = H0Z + H1Z
      HZT = H0ZT + H1ZT
      COMM = H + HRS + HRST + T2*HT/6.d0 + 7.d0*T2*T*HTT/6.d0
      PREF = HZ - GM*T2*HT/GZ
      FACT5 = GM* (2.d0*HT+T*HTT)/GZ
      COMM = COMM - PREF*ZTA - UU*HTT - VV*HT - WW* (HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF

c  local correlation option:
c     dvcup = 0.0d0
c     dvcdn = 0.0d0

      RETURN
      END
