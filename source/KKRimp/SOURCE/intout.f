!-------------------------------------------------------------------------------
!> Summary: 
!> Author: Who wrote this subroutine
!> Category: core electron, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!-------------------------------------------------------------------------------

      SUBROUTINE INTOUT(G,F,V,E,L,NNE,K2,DG,A,B,Z,NSRA)

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
      END
