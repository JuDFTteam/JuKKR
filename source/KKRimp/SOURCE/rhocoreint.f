      MODULE MOD_RHOCOREINT

      CONTAINS

      !----------------------------------------------------------------------
      !> Summary: Calculation of core states
      !> Category: core-electrons, KKRimp, initialization
      !>     
      !> Performs sum over L-channels.
      !> 
      !> Explanation of some arrays:
      !>     lmxc = lmaxcore = (0,1,2,...), .e.g, argon core : lmxc = 1
      !>                                        krypton core : lmxc = 2
      !>     kfg = configuration of core, e.g., argon core: 3300=3s,3p,0d
      !>                                      krypton core: 4430=4s,4p,3d
      !>                                        xenon core: 5540=5s,5p,4d
      !> @note Similar to routine 'corel' of KKRhost code @endnote
      !----------------------------------------------------------------------
      SUBROUTINE RHOCOREINT(NSRA,IPR,IP,RHOC,V,ECORE,
     +                      LCORE,NCORE,DRDI,Z,QC,
     +                      A,B,IS,NSPIN,NR,RMAX,IRMD)

C     .. Parameters ..
      USE MOD_SIMP3
      IMPLICIT NONE
      INTEGER NITMAX,IRNUMX
      PARAMETER (NITMAX=40,IRNUMX=10)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,QC,RMAX,Z
      INTEGER IP,IPR,IRMD,IS,NCORE,NR,NSPIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),ECORE(*),RHOC(*),V(*)
      INTEGER LCORE(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION E,E1,E2,EDIFF,EI,SLOPE,SUM,TOL,VALUE,WGT
      INTEGER IC,IN,INUC,IR,L,LMP1,LMXC,LP1,NC,NMAX,NN,NRE
      LOGICAL VLNC
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION F(IRMD),G(IRMD),RHO(IRMD)
      INTEGER KFG(4)
      CHARACTER*4 SPN(2),TEXT(5)
C     ..
C     .. External Subroutines ..
!        EXTERNAL RHOCOREINT_INTCOR !,RHOCOREINT_SIMP3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,REAL
C     ..
C     .. Save statement ..
      SAVE SPN,TEXT
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','s   ','p   ','d   ','f   ','g   '/
C     ..
      VLNC = .false.
      VALUE = 1.D-8
      SLOPE = -1.D-8
      E2 = 50.0D0
c
      DO 10 IC = 1,4
        KFG(IC) = 0
   10 CONTINUE
      DO 20 IC = 1,NCORE
        IF (LCORE(IC).EQ.0) KFG(1) = KFG(1) + 1
        IF (LCORE(IC).EQ.1) KFG(2) = KFG(2) + 1
        IF (LCORE(IC).EQ.2) KFG(3) = KFG(3) + 1
        IF (LCORE(IC).EQ.3) KFG(4) = KFG(4) + 1
   20 CONTINUE
      IF (KFG(2).NE.0) KFG(2) = KFG(2) + 1
      IF (KFG(3).NE.0) KFG(3) = KFG(3) + 2
      IF (KFG(4).NE.0) KFG(4) = KFG(4) + 3
      LMXC = 0
      IF (KFG(2).NE.0) LMXC = 1
      IF (KFG(3).NE.0) LMXC = 2
      IF (KFG(4).NE.0) LMXC = 3
c
      TOL = 1.0D-12* (Z*Z+1.D0)
      LMP1 = LMXC + 1
      NC = 0
      INUC = -IRNUMX
c
      DO 30 IR = 1,IRMD
        RHOC(IR) = ZERO
        RHO(IR) = ZERO
   30 CONTINUE
c
      DO 70 LP1 = 1,LMP1
        L = LP1 - 1
        E1 = (-5.D0- ((Z+1.D0)/DBLE(LP1))**2)*1.5D0 - 50.D0
        NMAX = KFG(LP1)
        IF (NMAX.NE.0) THEN
          DO 60 IN = LP1,NMAX
            NN = IN - LP1
            NC = NC + 1
            INUC = INUC + IRNUMX
            E = ECORE(NC)
            EI = ECORE(NC)
            IF (IPR.NE.0) WRITE (6,FMT=9000) IN,TEXT(LP1),NN,SPN(IS),
     +          IP,E
            CALL RHOCOREINT_INTCOR(E1,E2,RHO,G,F,V,VALUE,SLOPE,
     +                             L,NN,E,SUM,NRE,VLNC,A,B,Z,
     +                             RMAX,NR,TOL,IRMD,IPR,NITMAX,NSRA)
            EDIFF = E - EI
            ECORE(NC) = E
            WGT = REAL(L+L+1)/SUM*2.D0/REAL(NSPIN)
            IF (IPR.NE.0) WRITE (6,FMT=9010) EI,EDIFF,E
   40       CONTINUE
c
c---> sum up contributions to total core charge
c
            DO 50 IR = 2,NRE
              RHOC(IR) = RHOC(IR) + RHO(IR)*WGT
              RHO(IR) = ZERO
   50       CONTINUE
   60     CONTINUE
        END IF

   70 CONTINUE
      IF (NC*IRNUMX.GT.150 .OR. IRNUMX.GT.10) STOP 'corel'
c
c---> integrate core density to get core charge
c
      CALL SIMP3(RHOC,QC,1,NR,DRDI)

 9000 FORMAT (1x,90 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1,
     +       '  spin=',a4,i5,'th cell','    einput = ',1p,d16.8)
 9010 FORMAT (1x,'  einput =',1p,d16.8,'   eout - ein =',1p,d16.8,
     +       '   eoutput = ',1p,d16.8)
      END SUBROUTINE RHOCOREINT


      !----------------------------------------------------------------------
      !> Summary: Integrate core density
      !> Category: core-electrons, KKRimp
      !> 
      !> @note Similar to routine 'intcor' of KKRhost code @endnote
      !> @warning uses own (hardcoded) value of `cvlight` which is changed for `nsra==1` @endwarning
      !----------------------------------------------------------------------
      SUBROUTINE RHOCOREINT_INTCOR(F1,F2,RHO,G,F,V,VALUE,
     +                             SLOPE,L,NN,E,SUM,NRE,VLNC,A,B,Z,
     +                             RN,NR,TOL,IRM,IPR,NITMAX,NSRA)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,E,F1,F2,RN,SLOPE,SUM,TOL,VALUE,Z
      INTEGER IPR,IRM,L,NITMAX,NN,NR,NRE,NSRA
      LOGICAL VLNC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION F(*),G(*),RHO(*),V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ARG,CAPPAI,DOFE
      DOUBLE PRECISION CVLIGHT,DE,DG1,DG2,DPSI1,DPSI2,DRDIKC,E1,E2,EA,
     +                 GKC2,PI,PKC1,PKC2,PSI1,PSI2,Q,QKC1,QKC2,RATIO,
     +                 RATIO1,RE,RKC,RPB,SLOP,TSRME,VALU,VME,XXX,ZZ
      INTEGER IR,K,K2,KC,NITER,NNE,NREM1,NREM2
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX HL(6)
C     ..
C     .. External Subroutines ..
!       EXTERNAL HANKEL,INTIN,INTOUT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DCMPLX,DSQRT,EXP,LOG,MAX0,MIN0,REAL,SQRT
C     ..
      PI = 4.D0*ATAN(1.D0)
      ZZ = Z + Z
      CVLIGHT = 274.0720442D0
      IF (NSRA.EQ.1) CVLIGHT = 1.0D0
      EA = EXP(A)
      NITER = 0
      E1 = F1
      E2 = F2
      IF (IPR.EQ.2) WRITE (6,FMT=9000) L,NN,NR,F1,E,F2,TOL,VALUE,SLOPE
   10 CONTINUE
      NITER = NITER + 1
      DO 20 IR = 1,IRM
        G(IR) = 0.0D0
        F(IR) = 0.0D0
   20 CONTINUE
      IF (NITER.GT.NITMAX) THEN
        GO TO 80

      ELSE
        IF (E.LE.E1 .OR. E.GE.E2) E = .5D0* (E1+E2)
        NRE = NR
        IF (E.LE.-1.D-8) THEN
          TSRME = 2.D0*SQRT(-E)
          RE = (LOG(-TSRME*E/1.D-8)/TSRME-ZZ/E)*2.D0
          NRE = LOG(RE/B+1.D0)/A + 1.D0
          NRE = (NRE/2)*2 + 1
          NRE = MIN0(NRE,NR)
          NRE = MAX0(NRE,35)
        END IF
        XXX = 1.D0
        VALU = 1.D-1
        SLOP = -1.D-1
        IF (NRE.LT.NR .AND. NITER.EQ.1 .AND. IPR.NE.0) WRITE (6,
     +      FMT=9010)
        IF (NRE.GE.NR) THEN
          VALU = VALUE
          SLOP = SLOPE
          IF (.NOT.VLNC) THEN
c--->   single site  boundary condition
            VME = -E
            IF (NSRA.EQ.1) THEN
              CAPPAI = DCMPLX(0.d0,DSQRT(VME))
            ELSE
              CAPPAI = DCMPLX(0.d0,DSQRT((1.D0-
     +                 VME/CVLIGHT/CVLIGHT)*VME))
            END IF
            ARG = CAPPAI*RN
            CALL RHOCOREINT_HANKEL(HL,L+2,ARG)
            DOFE = REAL(L+1)/RN - CAPPAI*HL(L+2)/HL(L+1)
            VALU = 1.D-10
            SLOP = VALU*DOFE
          END IF

        END IF
        K2 = 30
        IF (NN.EQ.0) K2 = NRE/3
        NNE = 0

        CALL RHOCOREINT_INTIN(G,F,V,E,L,NNE,VALU,SLOP,
     +                        NRE,K2,KC,DG2,A,B,Z,NSRA)

        RKC = B*EXP(A*KC-A) - B
        DRDIKC = A* (RKC+B)
        GKC2 = G(KC)
        PSI2 = G(KC)
        DPSI2 = DG2/DRDIKC
        QKC2 = PSI2*PSI2 + DPSI2*DPSI2*RKC*RKC
        PKC2 = .5D0 - ATAN(RKC*DPSI2/PSI2)/PI

        CALL RHOCOREINT_INTOUT(G,F,V,E,L,NNE,KC,DG1,
     +                         A,B,Z,NSRA)

        PSI1 = G(KC)
        DPSI1 = DG1/DRDIKC
        QKC1 = PSI1*PSI1 + DPSI1*DPSI1*RKC*RKC
        PKC1 = .5D0 - ATAN(RKC*DPSI1/PSI1)/PI
        IF (NNE.EQ.9) NNE = 0
        IF (NNE.EQ.NN) THEN
          RATIO1 = GKC2/G(KC)
          RATIO = SQRT(QKC2/QKC1)
          IF (RATIO1.LT.0.D0) RATIO = -RATIO
          DO 30 K = 1,KC
            G(K) = G(K)*RATIO
            F(K) = F(K)*RATIO
   30     CONTINUE
          SUM = 0.D0
          IF (NSRA.EQ.1) THEN
            DO 40 K = 1,NRE
              F(K) = 0.0D0
   40       CONTINUE
          END IF
          RPB = B/EA
          Q = EA*EA
          NREM1 = NRE - 1
          DO 50 K = 2,NREM1,2
            RPB = RPB*Q
            SUM = SUM + RPB* (G(K)*G(K)+F(K)*F(K))
   50     CONTINUE
          RPB = B
          SUM = SUM + SUM
          NREM2 = NRE - 2
          DO 60 K = 3,NREM2,2
            RPB = RPB*Q
            SUM = SUM + RPB* (G(K)*G(K)+F(K)*F(K))
   60     CONTINUE
          SUM = SUM + SUM + RPB*Q* (G(NRE)*G(NRE)+F(NRE)*F(NRE))
          SUM = A*SUM/3.D0
          DE = PI*QKC2* (PKC2-PKC1)/SUM/RKC
          IF (NITER.GE.NITMAX-10 .OR. IPR.EQ.2) WRITE (6,
     +        FMT=9020) NITER,NNE,NRE,KC,E1,E,E2,DE
          IF (DE.GT.0.D0) E1 = E
          IF (DE.LT.0.D0) E2 = E
          E = E + DE
          IF (ABS(DE).GT.TOL .AND. NITER.LT.NITMAX) GO TO 10

        ELSE
          IF (NITER.GE.NITMAX-10 .OR. IPR.EQ.2) WRITE (6,
     +        FMT=9020) NITER,NNE,NRE,KC,E1,E,E2
          IF (NNE.GT.NN) E2 = E
          IF (NNE.LT.NN) E1 = E
          E = .5D0* (E1+E2)
          GO TO 10

        END IF

      END IF

      E = E - DE
      DO 70 K = 1,NRE
        RHO(K) = G(K)*G(K) + F(K)*F(K)
   70 CONTINUE
      IF (XXX.LE.0.D0) WRITE (6,FMT=9030)
      IF (NITER.GE.NITMAX-10 .OR. IPR.GE.1 .OR. XXX.LE.0.D0) WRITE (6,
     +    FMT=9040) L,NN,NITER,KC,NRE,VALU,SLOP,E,DE,SUM
      RETURN

   80 WRITE (6,FMT=9050) NITMAX
      STOP 'INTCOR'


 9000 FORMAT (' l=',i3,'  nn=',i2,'  nr=',i4,'  f1/e/f2=',3f10.3,/,
     +       ' tol=',1p,d12.3,'  value/slope=',2d12.3)
 9010 FORMAT (13x,'  no boundary condition had to be used')
 9020 FORMAT (2i3,2i4,1p,3d16.8,1x,2d9.2)
 9030 FORMAT (/,' **** int: 0-pressure bcs not real')
 9040 FORMAT (' state',i2,',',i1,':',i4,'x,',i5,'/',i3,',  bc=',1p,
     +       2d12.3,/,14x,'e=',d14.6,'   de=',d11.2,'   sum=',d12.4)
 9050 FORMAT (' *** int: stop after',i4,' iterations')
      END SUBROUTINE RHOCOREINT_INTCOR


      !----------------------------------------------------------------------
      !> Summary: 
      !> Category: core-electrons, KKRimp
      !> 
      !> @note Similar to routine 'intin' of KKRhost code @endnote
      !> @warning uses own (hardcoded) value of `cvlight` which is changed for `nsra==1` @endwarning
      !----------------------------------------------------------------------
      SUBROUTINE RHOCOREINT_INTIN(G,F,V,E,L,NNE,VALU,SLOP,
     +                            K1,K2,KC,DG,A,B,Z,NSRA)

C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,DG,E,SLOP,VALU,Z
      INTEGER K1,K2,KC,L,NNE,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION F(*),G(*),V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AF1,AF2,AF3,AG1,AG2,AG3,B1,B2,CVLIGHT,DET,DF1,
     +                 DF2,DF3,DG1,DG2,DG3,DR,EA,FF,FLLP1,GG,H83,PHI,Q,
     +                 R,R1,R2,R3,R83SQ,RPB,SDG3,SG,SGP1,U,VB,X,Y,ZZ
      INTEGER I,K,KP1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION D(2,3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSIGN,EXP,SQRT
C     ..
      ZZ = Z + Z
      CVLIGHT = 274.0720442D0
      IF (NSRA.EQ.1) CVLIGHT = 1.0D0
      FLLP1 = L* (L+1.D0)
      R83SQ = 64.D0/9.D0
      R1 = 1.D0/9.D0
      R2 = -5.D0*R1
      R3 = 19.D0*R1
      H83 = -8.D0/3.D0
      EA = EXP(A)
      RPB = B*EXP(A*K1-A)
      R = RPB - B
      DR = A*RPB
      PHI = (E+ZZ/R-V(K1))*DR/CVLIGHT
      U = DR*CVLIGHT + PHI
      IF (NSRA.EQ.1) U = DR
      X = -DR/R
      Y = -FLLP1*X*X/U + PHI
      G(K1) = VALU
      F(K1) = (SLOP*DR+X*VALU)/U
      Q = 1.D0/SQRT(EA)
      AG1 = SLOP*DR
      AF1 = X*F(K1) - Y*G(K1)
      K = K1
      DG3 = AG1
      IF (K2.NE.K1) THEN
        DO 10 I = 1,3
          KP1 = K
          K = K - 1
          RPB = RPB*Q
          DR = RPB*A
          R = RPB - B
          GG = G(KP1) - .5D0*AG1
          FF = F(KP1) - .5D0*AF1
          VB = (3.D0*V(KP1)+6.D0*V(K)-V(K-1))*.125D0
          PHI = (E+ZZ/R-VB)*DR/CVLIGHT
          U = DR*CVLIGHT + PHI
          IF (NSRA.EQ.1) U = DR
          X = -DR/R
          Y = -FLLP1*X*X/U + PHI
          AG2 = U*FF - X*GG
          AF2 = X*FF - Y*GG
          GG = G(KP1) - .5D0*AG2
          FF = F(KP1) - .5D0*AF2
          AG3 = U*FF - X*GG
          AF3 = X*FF - Y*GG
          RPB = RPB*Q
          DR = A*RPB
          R = RPB - B
          PHI = (E+ZZ/R-V(K))*DR/CVLIGHT
          U = DR*CVLIGHT + PHI
          IF (NSRA.EQ.1) U = DR
          X = -DR/R
          Y = -FLLP1*X*X/U + PHI
          GG = G(KP1) - AG3
          FF = F(KP1) - AF3
          G(K) = G(KP1) - (AG1+2.D0* (AG2+AG3)+U*FF-X*GG)/6.D0
          F(K) = F(KP1) - (AF1+2.D0* (AF2+AF3)+X*FF-Y*GG)/6.D0
          SG = DSIGN(1.0D0,G(K))
          SGP1 = DSIGN(1.0D0,G(KP1))
          IF (SG*SGP1.LT.0.D0) NNE = NNE + 1
          AG1 = U*F(K) - X*G(K)
          AF1 = X*F(K) - Y*G(K)
          IF (K.EQ.K2) THEN
            GO TO 30

          ELSE
            D(1,I) = AG1
            D(2,I) = AF1
          END IF

   10   CONTINUE
        Q = 1.D0/EA
        DG1 = D(1,1)
        DG2 = D(1,2)
        DG3 = D(1,3)
        DF1 = D(2,1)
        DF2 = D(2,2)
        DF3 = D(2,3)
   20   CONTINUE
        KP1 = K
        K = K - 1
        RPB = RPB*Q
        DR = A*RPB
        R = RPB - B
        PHI = (E+ZZ/R-V(K))*DR/CVLIGHT
        U = DR*CVLIGHT + PHI
        IF (NSRA.EQ.1) U = DR
        X = -DR/R
        Y = -FLLP1*X*X/U + PHI
        DET = R83SQ - X*X + U*Y
        B1 = G(KP1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
        B2 = F(KP1)*H83 + R1*DF1 + R2*DF2 + R3*DF3
        G(K) = (B1* (H83-X)+B2*U)/DET
        F(K) = (B2* (H83+X)-B1*Y)/DET
        SG = DSIGN(1.0D0,G(K))
        SGP1 = DSIGN(1.0D0,G(KP1))
        IF (SG*SGP1.LT.0.D0) NNE = NNE + 1
        DG1 = DG2
        DF1 = DF2
        DG2 = DG3
        DF2 = DF3
        DG3 = U*F(K) - X*G(K)
        DF3 = X*F(K) - Y*G(K)
        SDG3 = DSIGN(1.0D0,DG3)
        IF (K.GT.K2 .AND. SG*SDG3.LT.0.D0) GO TO 20
      END IF
   30 KC = K
      DG = DG3
      END SUBROUTINE RHOCOREINT_INTIN


      !-------------------------------------------------------------------------------
      !> Summary:
      !> Author:
      !> Category: KKRimp, core-electrons
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !> @note Similar to routine 'intout' of KKRhost code @endnote
      !> @warning uses own (hardcoded) value of `cvlight` which is changed for `nsra==1` @endwarning
      !-------------------------------------------------------------------------------
      SUBROUTINE RHOCOREINT_INTOUT(G,F,V,E,L,NNE,K2,DG,A,B,Z,NSRA)

C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,DG,E,Z
      INTEGER K2,L,NNE,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION F(*),G(*),V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AA,ALFA,B1,B2,BB,BETA,CVLIGHT,DET,DF1,DF2,DF3,
     +                 DG1,DG2,DG3,DR,EA,FLLP1,H83,P12,P21,PHI,PP,QQ,R,
     +                 R1,R2,R3,R83SQ,RPB,S,SG,SGM1,U,X,Y,ZZ
      INTEGER I,I1,K,KM1,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION D(2,3),PX(20),QX(20)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSIGN,EXP,SQRT
C     ..
      ZZ = Z + Z
      CVLIGHT = 274.0720442D0
      IF (NSRA.EQ.1) CVLIGHT = 1.0D0
      EA = EXP(A)
      FLLP1 = L* (L+1.D0)
      R83SQ = 64.D0/9.D0
      R1 = 1.D0/9.D0
      R2 = -5.D0*R1
      R3 = 19.D0*R1
      H83 = 8.D0/3.D0
      AA = -ZZ/CVLIGHT
      BB = FLLP1 - AA*AA
      P21 = (V(1)-E)/CVLIGHT
      P12 = CVLIGHT - P21
      PX(1) = 0.D0
      QX(1) = 0.D0
      IF (Z.LE.20.D0 .OR. NSRA.EQ.1) THEN
        S = 1.D0*L
        PX(2) = 0.D0
        PX(3) = 1.D0
        DO 10 K = 2,9
          PX(K+2) = ((V(1)-E)*PX(K)-ZZ*PX(K+1))/ (K+L+L)/ (K-1.D0)
   10   CONTINUE
        DO 20 K = 2,10
          QX(K) = PX(K+1)* (L+K-2.D0)/CVLIGHT
   20   CONTINUE

      ELSE
        S = SQRT(FLLP1+1.D0-AA*AA)
        PX(2) = 1.D0
        QX(2) = (1.D0-S)/AA
        DO 30 I = 3,10
          N = I - 2
          ALFA = P12*QX(I-1)
          BETA = P21*AA*PX(I-1)
          IF (L.NE.0) BETA = BETA - P12*AA*PX(I-1) - P12*P21*PX(I-2) +
     +                       (N+S)*P12*QX(I-1)
          DET = N* (N+S+S)*AA
          PX(I) = (ALFA* (N+S+1.D0)*AA-AA*BETA)/DET
          QX(I) = (BETA* (N+S-1.D0)-BB*ALFA)/DET
   30   CONTINUE
      END IF

      G(1) = 0.D0
      F(1) = 0.D0
      RPB = B
      DO 50 K = 2,4
        RPB = RPB*EA
        R = RPB - B
        DR = A*RPB
        PHI = (E+ZZ/R-V(K))*DR/CVLIGHT
        U = DR*CVLIGHT + PHI
        IF (NSRA.EQ.1) U = DR
        X = -DR/R
        Y = -FLLP1*X*X/U + PHI
        PP = PX(10)
        QQ = QX(10)
        DO 40 I1 = 3,10
          I = 12 - I1
          PP = PP*R + PX(I)
          QQ = QQ*R + QX(I)
   40   CONTINUE
        G(K) = (R**S)*PP
        F(K) = (R**S)*QQ
        SG = DSIGN(1.0D0,G(K))
        SGM1 = DSIGN(1.0D0,G(K-1))
        IF (SG*SGM1.LT.0.D0) NNE = NNE + 1
        D(1,K-1) = U*F(K) - X*G(K)
        D(2,K-1) = X*F(K) - Y*G(K)
   50 CONTINUE
      DG1 = D(1,1)
      DG2 = D(1,2)
      DG3 = D(1,3)
      DF1 = D(2,1)
      DF2 = D(2,2)
      DF3 = D(2,3)
      DO 60 K = 5,K2
        KM1 = K - 1
        RPB = RPB*EA
        R = RPB - B
        DR = A*RPB
        PHI = (E+ZZ/R-V(K))*DR/CVLIGHT
        U = DR*CVLIGHT + PHI
        IF (NSRA.EQ.1) U = DR
        X = -DR/R
        Y = -FLLP1*X*X/U + PHI
        DET = R83SQ - X*X + U*Y
        B1 = G(KM1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
        B2 = F(KM1)*H83 + R1*DF1 + R2*DF2 + R3*DF3
        G(K) = (B1* (H83-X)+B2*U)/DET
        F(K) = (B2* (H83+X)-B1*Y)/DET
        SG = DSIGN(1.0D0,G(K))
        SGM1 = DSIGN(1.0D0,G(KM1))
        IF (SG*SGM1.LT.0.D0) NNE = NNE + 1
        DG1 = DG2
        DG2 = DG3
        DG3 = U*F(K) - X*G(K)
        DF1 = DF2
        DF2 = DF3
        DF3 = X*F(K) - Y*G(K)
   60 CONTINUE
      DG = DG3
      END SUBROUTINE RHOCOREINT_INTOUT

      !-------------------------------------------------------------------------------
      !> Summary:
      !> Author:
      !> Category: KKRimp, core-electrons
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !> This subroutine uses the explicit formulas for the hankel
      !> functions. for higher l-values these formulas may lead to
      !> loss of significant figures. This subroutine should be used
      !> only for core states.
      !>
      !> @note Similar to routine 'hankel' of KKRhost code @endnote
      !-------------------------------------------------------------------------------
      SUBROUTINE RHOCOREINT_HANKEL(H,L,ARG)

C     .. Scalar Arguments ..
      DOUBLE COMPLEX ARG
      INTEGER L
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX H(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX A1,A2,A3,A4
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     ..
C     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
      H(1) = -EXP(ARG*CI)/ARG
      IF (L.NE.1) THEN
        A1 = (1.D0,0.D0) - ARG*CI
        H(2) = H(1)*A1/ARG
        IF (L.NE.2) THEN
          A1 = 3.D0*A1
          A2 = ARG*ARG
          H(3) = H(1)* (A1-A2)/A2
          IF (L.NE.3) THEN
            A1 = 5.D0*A1
            A3 = A2*ARG*CI
            A4 = A2*ARG
            A2 = 6.D0*A2
            H(4) = H(1)* (A1-A2+A3)/A4
            IF (L.NE.4) THEN
              A1 = 7.D0*A1
              A2 = 7.5D0*A2
              A3 = 10.D0*A3
              A4 = A4*ARG
              H(5) = H(1)* (A1-A2+A3+A4)/A4
              IF (L.NE.5) THEN
                H(6) = (9.0D0,0.0D0)*H(5)/ARG - H(4)
                IF (L.NE.6) THEN
                  WRITE (6,FMT=9000) L
                  STOP 'HANKEL'

                END IF

              END IF

            END IF

          END IF

        END IF

      END IF

      RETURN


 9000 FORMAT (2x,' hankel :  l=',i2,' is too large')
      END SUBROUTINE RHOCOREINT_HANKEL

      END MODULE MOD_RHOCOREINT



