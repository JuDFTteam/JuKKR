!-------------------------------------------------------------------------------
!> Summary: This subroutine is used for calcualtion of the core charge density  
!> Author: Who wrote this subroutine
!> Category: core-electrons
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------
!> @note Similar to the intcor.f of  the host kkr code. 
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!-------------------------------------------------------------------------------


      SUBROUTINE INTCOR(F1,F2,RHO,G,F,V,VALUE,SLOPE,L,NN,E,SUM,NRE,
     +                    VLNC,A,B,Z,RN,NR,TOL,IRM,IPR,NITMAX,NSRA)
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
      EXTERNAL HANKEL,INTIN,INTOUT
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
            CALL HANKEL(HL,L+2,ARG)
            DOFE = REAL(L+1)/RN - CAPPAI*HL(L+2)/HL(L+1)
            VALU = 1.D-10
            SLOP = VALU*DOFE
          END IF

        END IF
        K2 = 30
        IF (NN.EQ.0) K2 = NRE/3
        NNE = 0

        CALL INTIN(G,F,V,E,L,NNE,VALU,SLOP,NRE,K2,KC,DG2,A,B,Z,NSRA)

        RKC = B*EXP(A*KC-A) - B
        DRDIKC = A* (RKC+B)
        GKC2 = G(KC)
        PSI2 = G(KC)
        DPSI2 = DG2/DRDIKC
        QKC2 = PSI2*PSI2 + DPSI2*DPSI2*RKC*RKC
        PKC2 = .5D0 - ATAN(RKC*DPSI2/PSI2)/PI

        CALL INTOUT(G,F,V,E,L,NNE,KC,DG1,A,B,Z,NSRA)

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
      END
