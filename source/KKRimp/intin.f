!-------------------------------------------------------------------------------
!> Summary: this subroutine is used for calculation of the core electron charge density
!>           
!> Author: Who wrote this subroutine
!> Category: core electron, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------



      SUBROUTINE INTIN(G,F,V,E,L,NNE,VALU,SLOP,K1,K2,KC,DG,A,B,Z,NSRA)

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
      END
