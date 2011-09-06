      SUBROUTINE GRADR(NSPIN,IST1,MESH,DX,DRDI,DRDI2,RO,ZTA,DRR,DDRR,
     +                 DRRU,DDRRU,ROU,IRMD)
C.....-----------------------------------------------------------------
C     evaluates d(ro)/dr,d{d(ro)/dr}/dr.
C     drr=d(ro)/dr, ddrr=d(drr)/dr.
C     coded by T.Asada. Feb.1994.
C.....-----------------------------------------------------------------
C.....-----------------------------------------------------------------
c.....------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION DX
      INTEGER IRMD,IST1,MESH,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DDRR(IRMD),DDRRU(IRMD),DRDI(IRMD),DRDI2(IRMD),
     +                 DRR(IRMD),DRRU(IRMD),RO(IRMD),ROU(IRMD),ZTA(IRMD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION D,DRX,DRX0,DRX1,DRX2,DRX3,DRXU,DRXU0,DRXU1,DRXU2,
     +                 DRXU3,DRXX,DRXX0,DRXX1,DRXX2,DRXX3,DRXXU,DRXXU0,
     +                 DRXXU1,DRXXU2,DRXXU3,F0,F1,F2,F3,F4,F5,G1,G2,G3,
     +                 G4,G5,XLF
      INTEGER I,I1,I2,I3,I4,I5,I6,IBH,ICA,ICG,IEX,IGD,IGH,IGL,IHB,IMJ,
     +        IP9,IPG,IPW,IST,IVG,IVN,IWR,IXLF,J,NDVPT,NRED
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION F131,F132,F133,F141,F142,F143,F144,F151,F152,
     +                 F153,F154,F155,F161,F162,F163,F164,F165,F166,
     +                 F231,F232,F233,F241,F242,F243,F244,F251,F252,
     +                 F253,F254,F255,F261,F262,F263,F264,F265,F266
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
C     .. Save statement ..
      SAVE NDVPT,IGL,IGH,IMJ,IBH,ICA,ICG,IVN,IPW,IPG,IVG,IP9,IGD,IXLF,
     +     IEX,XLF,IWR
C     ..
C     .. Data statements ..
c.....-----------------------------------------------------------------
C     double function
C      double precision f131,f132,f133,f141,f142,f143,f144
C      double precision fl61,fl62,fl63,fl64,fl65,fl66
C      double precision fl51,fl52,fl53,fl54,fl55
C      double precision f231,f232,f233,f241,f242,f243,f244
C      double precision f251,f252,f253,f254,f255
C      double precision f261,f262,f263,f264,f265,f266

      DATA NDVPT/6/
      DATA IGL,IGH,IMJ,IBH,ICA,ICG,IVN,IPW,IPG,IVG,IP9,IGD,IXLF,IEX,
     +     XLF/0,0,0,0,0,1,0,0,1,0,0,1,0,0,0.00D0/
      DATA IWR/0/
C     ..
C     .. Statement Function definitions ..
c  statement functions:

c.....three point formula for the 1st deriv.

c.....four point formula for the 1st deriv.

c.....five point formula for the 1st deriv.

c.....six point formula for the 1st deriv.

c.....three point formula for the 2nd deriv.

c.....four point formula for the 2nd deriv.

c.....five point formula for the 2nd deriv.

c.....six point formula for the 2nd deriv.
      F131(F0,F1,F2,D) = (-3*F0+4*F1-F2)/ (2*D)
      F132(G1,F0,F1,D) = (-1*G1-0*F0+F1)/ (2*D)
      F133(G2,G1,F0,D) = (G2-4*G1+3*F0)/ (2*D)
      F141(F0,F1,F2,F3,D) = (-11*F0+18*F1-9*F2+2*F3)/ (6*D)
      F142(G1,F0,F1,F2,D) = (-2*G1-3*F0+6*F1-F2)/ (6*D)
      F143(G2,G1,F0,F1,D) = (G2-6*G1+3*F0+2*F1)/ (6*D)
      F144(G3,G2,G1,F0,D) = (-2*G3+9*G2-18*G1+11*F0)/ (6*D)
      F151(F0,F1,F2,F3,F4,D) = (-50*F0+96*F1-72*F2+32*F3-6*F4)/ (24*D)
      F152(G1,F0,F1,F2,F3,D) = (-6*G1-20*F0+36*F1-12*F2+2*F3)/ (24*D)
      F153(G2,G1,F0,F1,F2,D) = (2*G2-16*G1-0*F0+16*F1-2*F2)/ (24*D)
      F154(G3,G2,G1,F0,F1,D) = (-2*G3+12*G2-36*G1+20*F0+6*F1)/ (24*D)
      F155(G4,G3,G2,G1,F0,D) = (6*G4-32*G3+72*G2-96*G1+50*F0)/ (24*D)
      F161(F0,F1,F2,F3,F4,F5,D) = (-274*F0+600*F1-600*F2+400*F3-150*F4+
     +                            24*F5)/ (120*D)
      F162(G1,F0,F1,F2,F3,F4,D) = (-24*G1-130*F0+240*F1-120*F2+40*F3-
     +                            6*F4)/ (120*D)
      F163(G2,G1,F0,F1,F2,F3,D) = (6*G2-60*G1-40*F0+120*F1-30*F2+4*F3)/
     +                            (120*D)
      F164(G3,G2,G1,F0,F1,F2,D) = (-4*G3+30*G2-120*G1+40*F0+60*F1-6*F2)/
     +                            (120*D)
      F165(G4,G3,G2,G1,F0,F1,D) = (6*G4-40*G3+120*G2-240*G1+130*F0+
     +                            24*F1)/ (120*D)
      F166(G5,G4,G3,G2,G1,F0,D) = (-24*G5+150*G4-400*G3+600*G2-600*G1+
     +                            274*F0)/ (120*D)
      F231(F0,F1,F2,D) = (F0-2*F1+F2)/ (D*D)
      F232(G1,F0,F1,D) = (G1-2*F0+F1)/ (D*D)
      F233(G2,G1,F0,D) = (G2-2*G1+F0)/ (D*D)
      F241(F0,F1,F2,F3,D) = (6*F0-15*F1+12*F2-3*F3)/ (3*D*D)
      F242(G1,F0,F1,F2,D) = (3*G1-6*F0+3*F1+0*F2)/ (3*D*D)
      F243(G2,G1,F0,F1,D) = (0*G2+3*G1-6*F0+3*F1)/ (3*D*D)
      F244(G3,G2,G1,F0,D) = (-3*G3+2*G2+15*G1+6*F0)/ (3*D*D)
      F251(F0,F1,F2,F3,F4,D) = (35*F0-104*F1+114*F2-56*F3+11*F4)/
     +                         (12*D*D)
      F252(G1,F0,F1,F2,F3,D) = (11*G1-20*F0+6*F1+4*F2-F3)/ (12*D*D)
      F253(G2,G1,F0,F1,F2,D) = (-G2+16*G1-30*F0+16*F1-F2)/ (12*D*D)
      F254(G3,G2,G1,F0,F1,D) = (-G3+4*G2+6*G1-20*F0+11*F1)/ (12*D*D)
      F255(G4,G3,G2,G1,F0,D) = (11*G4-56*G3+114*G2-104*G1+35*F0)/
     +                         (12*D*D)
      F261(F0,F1,F2,F3,F4,F5,D) = (225*F0-770*F1+1070*F2-780*F3+305*F4-
     +                            50*F5)/ (60*D*D)
      F262(G1,F0,F1,F2,F3,F4,D) = (50*G1-75*F0-20*F1+70*F2-30*F3+5*F4)/
     +                            (60*D*D)
      F263(G2,G1,F0,F1,F2,F3,D) = (-5*G2+80*G1-150*F0+80*F1-5*F2+0*F3)/
     +                            (60*D*D)
      F264(G3,G2,G1,F0,F1,F2,D) = (0*G3-5*G2+80*G1-150*F0+80*F1-5*F2)/
     +                            (60*D*D)
      F265(G4,G3,G2,G1,F0,F1,D) = (5*G4-30*G3+70*G2-20*G1-75*F0+50*F1)/
     +                            (60*D*D)
      F266(G5,G4,G3,G2,G1,F0,D) = (-50*G5+305*G4-780*G3+1070*G2-770*G1+
     +                            225*F0)/ (60*D*D)
C     ..

c.....-----------------------------------------------------------------
      IF (IWR.EQ.1) WRITE (6,FMT=
     +'(/''  igl,igh,imj,ihb,ica,icg,ivn,ipw,ipg,'',            ''ivg,ip
     +9,igd,ixlf,iex,xlf='',14i2,f10.4)') IGL,IGH,IMJ,IHB,ICA,ICG,IVN,
     +    IPW,IPG,IVG,IP9,IGD,IXLF,IEX,XLF
      IWR = 0

      IST = IST1
c     write(6,*) 'ndvpt ist mesh dx drdi2' ,ndvpt,ist,mesh,dx,
c    &            drdi2(ist)

      IF (NDVPT.LT.3 .OR. NDVPT.GT.6) THEN
        WRITE (6,FMT=9000) NDVPT
        STOP 18
      END IF
c.....
c.....ro: total(core+val)(up+down) charge density.

      DO 10 I = IST,MESH
        ROU(I) = RO(I)* (ZTA(I)+1.D0)/2.D0
   10 CONTINUE
c.....
      IF (IGD.LE.0) THEN

        DO 20 I = IST,MESH
          DRR(I) = 0.D0
          DDRR(I) = 0.D0
          DRRU(I) = 0.D0
          DDRRU(I) = 0.D0
   20   CONTINUE
        GO TO 60

      END IF

      I1 = IST
      I2 = IST + 1
      I3 = IST + 2
      I4 = IST + 3
      I5 = IST + 4
      I6 = IST + 5

c.....drr:d(ro)/dr, ddrr=d(d(ro)/dr)/dr
cc.... drru,ddrru: for up   spin,
c.....

      IF (NSPIN.EQ.1) GO TO 40
c.....
      IF (NDVPT.EQ.3) THEN

        DRX1 = F131(RO(I1),RO(I2),RO(I3),DX)
        DRXU1 = F131(ROU(I1),ROU(I2),ROU(I3),DX)
        DRXX1 = F231(RO(I1),RO(I2),RO(I3),DX)
        DRXXU1 = F231(ROU(I1),ROU(I2),ROU(I3),DX)

      ELSE IF (NDVPT.EQ.4) THEN

        DRX1 = F141(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRXU1 = F141(ROU(I1),ROU(I2),ROU(I3),ROU(I4),DX)
        DRXX1 = F241(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRXXU1 = F241(ROU(I1),ROU(I2),ROU(I3),ROU(I4),DX)
        DRX2 = F142(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRXU2 = F142(ROU(I1),ROU(I2),ROU(I3),ROU(I4),DX)
        DRXX2 = F242(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRXXU2 = F242(ROU(I1),ROU(I2),ROU(I3),ROU(I4),DX)

      ELSE IF (NDVPT.EQ.5) THEN

        DRX1 = F151(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRXU1 = F151(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),DX)
        DRXX1 = F251(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRXXU1 = F251(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),DX)
        DRX2 = F152(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRXU2 = F152(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),DX)
        DRXX2 = F252(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRXXU2 = F252(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),DX)

      ELSE IF (NDVPT.EQ.6) THEN

        DRX1 = F161(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXU1 = F161(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),ROU(I6),DX)
        DRXX1 = F261(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXXU1 = F261(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),ROU(I6),
     +           DX)
        DRX2 = F162(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXU2 = F162(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),ROU(I6),DX)
        DRXX2 = F262(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXXU2 = F262(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),ROU(I6),
     +           DX)
        DRX3 = F163(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXU3 = F163(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),ROU(I6),DX)
        DRXX3 = F263(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXXU3 = F263(ROU(I1),ROU(I2),ROU(I3),ROU(I4),ROU(I5),ROU(I6),
     +           DX)

      END IF

      DRR(I1) = DRX1/DRDI(I1)
      DDRR(I1) = (DRXX1-DRX1*DRDI2(I1))/DRDI(I1)**2
      DRRU(I1) = DRXU1/DRDI(I1)
      DDRRU(I1) = (DRXXU1-DRXU1*DRDI2(I1))/DRDI(I1)**2

      IF (NDVPT.GT.3) THEN

        DRR(I2) = DRX2/DRDI(I2)
        DDRR(I2) = (DRXX2-DRX2*DRDI2(I2))/DRDI(I2)**2
        DRRU(I2) = DRXU2/DRDI(I2)
        DDRRU(I2) = (DRXXU2-DRXU2*DRDI2(I2))/DRDI(I2)**2

        IF (NDVPT.EQ.6) THEN
          DRR(I3) = DRX3/DRDI(I3)
          DDRR(I3) = (DRXX3-DRX3*DRDI2(I3))/DRDI(I3)**2
          DRRU(I3) = DRXU3/DRDI(I3)
          DDRRU(I3) = (DRXXU3-DRXU3*DRDI2(I3))/DRDI(I3)**2
        END IF

      END IF

      NRED = DBLE(NDVPT)/2 + .1D0

      DO 30 J = NRED + IST,MESH - NRED

        IF (NDVPT.EQ.3) THEN

          DRX = F132(RO(J-1),RO(J),RO(J+1),DX)
          DRXU = F132(ROU(J-1),ROU(J),ROU(J+1),DX)
          DRXX = F232(RO(J-1),RO(J),RO(J+1),DX)
          DRXXU = F232(ROU(J-1),ROU(J),ROU(J+1),DX)

        ELSE IF (NDVPT.EQ.4) THEN

          DRX = F142(RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXU = F142(ROU(J-1),ROU(J),ROU(J+1),ROU(J+2),DX)
          DRXX = F242(RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXXU = F242(ROU(J-1),ROU(J),ROU(J+1),ROU(J+2),DX)

        ELSE IF (NDVPT.EQ.5) THEN

          DRX = F153(RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXU = F153(ROU(J-2),ROU(J-1),ROU(J),ROU(J+1),ROU(J+2),DX)
          DRXX = F253(RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXXU = F253(ROU(J-2),ROU(J-1),ROU(J),ROU(J+1),ROU(J+2),DX)

        ELSE IF (NDVPT.EQ.6) THEN

          DRX = F164(RO(J-3),RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXU = F164(ROU(J-3),ROU(J-2),ROU(J-1),ROU(J),ROU(J+1),
     +           ROU(J+2),DX)
          DRXX = F264(RO(J-3),RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXXU = F264(ROU(J-3),ROU(J-2),ROU(J-1),ROU(J),ROU(J+1),
     +            ROU(J+2),DX)

        END IF

        DRR(J) = DRX/DRDI(J)
        DDRR(J) = (DRXX-DRX*DRDI2(J))/DRDI(J)**2
        DRRU(J) = DRXU/DRDI(J)
        DDRRU(J) = (DRXXU-DRXU*DRDI2(J))/DRDI(J)**2

   30 CONTINUE
c.....
      IF (NDVPT.EQ.3) THEN

        DRX0 = F133(RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXU0 = F133(ROU(MESH-2),ROU(MESH-1),ROU(MESH),DX)
        DRXX0 = F233(RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXXU0 = F233(ROU(MESH-2),ROU(MESH-1),ROU(MESH),DX)

      ELSE IF (NDVPT.EQ.4) THEN

        DRX1 = F143(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXU1 = F143(ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),ROU(MESH),DX)
        DRXX1 = F243(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXXU1 = F243(ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),ROU(MESH),DX)
        DRX0 = F144(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXU0 = F144(ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),ROU(MESH),DX)
        DRXX0 = F244(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXXU0 = F244(ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),ROU(MESH),DX)

      ELSE IF (NDVPT.EQ.5) THEN

        DRX1 = F154(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +         RO(MESH),DX)
        DRXU1 = F154(ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),
     +          ROU(MESH),DX)
        DRXX1 = F254(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +          RO(MESH),DX)
        DRXXU1 = F254(ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),
     +           ROU(MESH),DX)
        DRX0 = F155(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +         RO(MESH),DX)
        DRXU0 = F155(ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),
     +          ROU(MESH),DX)
        DRXX0 = F255(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +          RO(MESH),DX)
        DRXXU0 = F255(ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),ROU(MESH-1),
     +           ROU(MESH),DX)

      ELSE IF (NDVPT.EQ.6) THEN

        DRX2 = F164(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +         RO(MESH-1),RO(MESH),DX)
        DRXU2 = F164(ROU(MESH-5),ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),
     +          ROU(MESH-1),ROU(MESH),DX)
        DRXX2 = F264(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +          RO(MESH-1),RO(MESH),DX)
        DRXXU2 = F264(ROU(MESH-5),ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),
     +           ROU(MESH-1),ROU(MESH),DX)

        DRX1 = F165(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +         RO(MESH-1),RO(MESH),DX)
        DRXU1 = F165(ROU(MESH-5),ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),
     +          ROU(MESH-1),ROU(MESH),DX)
        DRXX1 = F265(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +          RO(MESH-1),RO(MESH),DX)
        DRXXU1 = F265(ROU(MESH-5),ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),
     +           ROU(MESH-1),ROU(MESH),DX)

        DRX0 = F166(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +         RO(MESH-1),RO(MESH),DX)
        DRXU0 = F166(ROU(MESH-5),ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),
     +          ROU(MESH-1),ROU(MESH),DX)
        DRXX0 = F266(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +          RO(MESH-1),RO(MESH),DX)
        DRXXU0 = F266(ROU(MESH-5),ROU(MESH-4),ROU(MESH-3),ROU(MESH-2),
     +           ROU(MESH-1),ROU(MESH),DX)


      END IF

      IF (NDVPT.GT.3) THEN

        IF (NDVPT.EQ.6) THEN
          DRR(MESH-2) = DRX2/DRDI(MESH-2)
          DRRU(MESH-2) = DRXU2/DRDI(MESH-2)
          DDRR(MESH-2) = (DRXX2-DRX2*DRDI2(MESH-2))/DRDI(MESH-2)**2
          DDRRU(MESH-2) = (DRXXU2-DRXU2*DRDI2(MESH-2))/DRDI(MESH-2)**2
        END IF

        DRR(MESH-1) = DRX1/DRDI(MESH-1)
        DRRU(MESH-1) = DRXU1/DRDI(MESH-1)
        DDRR(MESH-1) = (DRXX1-DRX1*DRDI2(MESH-1))/DRDI(MESH-1)**2
        DDRRU(MESH-1) = (DRXXU1-DRXU1*DRDI2(MESH-1))/DRDI(MESH-1)**2

      END IF

      DRR(MESH) = DRX0/DRDI(MESH)
      DRRU(MESH) = DRXU0/DRDI(MESH)
      DDRR(MESH) = (DRXX0-DRX0*DRDI2(MESH))/DRDI(MESH)**2
      DDRRU(MESH) = (DRXXU0-DRXU0*DRDI2(MESH))/DRDI(MESH)**2

      GO TO 60

   40 CONTINUE

c.....
      IF (NDVPT.EQ.3) THEN

        DRX1 = F131(RO(I1),RO(I2),RO(I3),DX)
        DRXX1 = F231(RO(I1),RO(I2),RO(I3),DX)

      ELSE IF (NDVPT.EQ.4) THEN

        DRX1 = F141(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRXX1 = F241(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRX2 = F142(RO(I1),RO(I2),RO(I3),RO(I4),DX)
        DRXX2 = F242(RO(I1),RO(I2),RO(I3),RO(I4),DX)

      ELSE IF (NDVPT.EQ.5) THEN

        DRX1 = F151(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRXX1 = F251(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRX2 = F152(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)
        DRXX2 = F252(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),DX)

      ELSE IF (NDVPT.EQ.6) THEN

        DRX1 = F161(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXX1 = F261(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRX2 = F162(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXX2 = F262(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRX3 = F163(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)
        DRXX3 = F263(RO(I1),RO(I2),RO(I3),RO(I4),RO(I5),RO(I6),DX)

      END IF

      DRR(I1) = DRX1/DRDI(I1)
      DDRR(I1) = (DRXX1-DRX1*DRDI2(I1))/DRDI(I1)**2

      IF (NDVPT.GT.3) THEN

        DRR(I2) = DRX2/DRDI(I2)
        DDRR(I2) = (DRXX2-DRX2*DRDI2(I2))/DRDI(I2)**2

        IF (NDVPT.EQ.6) THEN
          DRR(I3) = DRX3/DRDI(I3)
          DDRR(I3) = (DRXX3-DRX3*DRDI2(I3))/DRDI(I3)**2
        END IF

      END IF

      NRED = DBLE(NDVPT)/2 + .1D0

      IF (MESH-NRED.LE.IST) THEN
         WRITE (6,FMT='(/'' mesh-nred.lt.ist. mesh,nred,ist='',3i4)')
     +        MESH,NRED,IST
         STOP 13
      END IF

      DO 50 J = NRED + IST,MESH - NRED

        IF (NDVPT.EQ.3) THEN

          DRX = F132(RO(J-1),RO(J),RO(J+1),DX)
          DRXX = F232(RO(J-1),RO(J),RO(J+1),DX)

        ELSE IF (NDVPT.EQ.4) THEN

          DRX = F142(RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXX = F242(RO(J-1),RO(J),RO(J+1),RO(J+2),DX)

        ELSE IF (NDVPT.EQ.5) THEN

          DRX = F153(RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXX = F253(RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)

        ELSE IF (NDVPT.EQ.6) THEN

          DRX = F164(RO(J-3),RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)
          DRXX = F264(RO(J-3),RO(J-2),RO(J-1),RO(J),RO(J+1),RO(J+2),DX)

        END IF

        DRR(J) = DRX/DRDI(J)
        DDRR(J) = (DRXX-DRX*DRDI2(J))/DRDI(J)**2
c           write(6,9000) j,drr(j)
c9000       format(1x,' j drr(j)',i5,e15.5)
   50 CONTINUE
c.....
      IF (NDVPT.EQ.3) THEN

        DRX0 = F133(RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXX0 = F233(RO(MESH-2),RO(MESH-1),RO(MESH),DX)

      ELSE IF (NDVPT.EQ.4) THEN

        DRX1 = F143(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXX1 = F243(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRX0 = F144(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)
        DRXX0 = F244(RO(MESH-3),RO(MESH-2),RO(MESH-1),RO(MESH),DX)

      ELSE IF (NDVPT.EQ.5) THEN

        DRX1 = F154(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +         RO(MESH),DX)
        DRXX1 = F254(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +          RO(MESH),DX)
        DRX0 = F155(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +         RO(MESH),DX)
        DRXX0 = F255(RO(MESH-4),RO(MESH-3),RO(MESH-2),RO(MESH-1),
     +          RO(MESH),DX)

      ELSE IF (NDVPT.EQ.6) THEN

        DRX2 = F164(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +         RO(MESH-1),RO(MESH),DX)
        DRXX2 = F264(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +          RO(MESH-1),RO(MESH),DX)

        DRX1 = F165(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +         RO(MESH-1),RO(MESH),DX)
        DRXX1 = F265(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +          RO(MESH-1),RO(MESH),DX)

        DRX0 = F166(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +         RO(MESH-1),RO(MESH),DX)
        DRXX0 = F266(RO(MESH-5),RO(MESH-4),RO(MESH-3),RO(MESH-2),
     +          RO(MESH-1),RO(MESH),DX)


      END IF

      IF (NDVPT.GT.3) THEN

        IF (NDVPT.EQ.6) THEN
          DRR(MESH-2) = DRX2/DRDI(MESH-2)
          DDRR(MESH-2) = (DRXX2-DRX2*DRDI2(MESH-2))/DRDI(MESH-2)**2
        END IF

        DRR(MESH-1) = DRX1/DRDI(MESH-1)
        DDRR(MESH-1) = (DRXX1-DRX1*DRDI2(MESH-1))/DRDI(MESH-1)**2

      END IF

      DRR(MESH) = DRX0/DRDI(MESH)
      DDRR(MESH) = (DRXX0-DRX0*DRDI2(MESH))/DRDI(MESH)**2

   60 CONTINUE
C
C
C
C      write(6,8000) nspin,ist1,mesh,dx
C8000 format(1x,' nspin ist1 mesh dx',3i5,2d20.10)
C     write(6,8001) (ro(kk),drr(kk),ddrr(kk),
c    &  drdi(kk),drdi2(kk), kk=ist1,mesh,20)
c8001 format(1x,' ro drr ddrr drdi drdi2',5f12.5)
      RETURN
 9000 FORMAT (/,' ndvpt should be ge.4 .or. le.6. ndvpt=',i3)
      END
