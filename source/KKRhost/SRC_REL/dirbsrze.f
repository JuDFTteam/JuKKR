      SUBROUTINE DIRBSRZE(IEST,XEST,YEST,YZ,DY,NV,NUSE)
C
C   ********************************************************************
C   *                                                                  *
C   *   diagonal rational function extrapolation to support the        *
C   *   Burlisch-Stoer method                                          *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C
C PARAMETER definitions
C
      INTEGER NCFMAX,ISEQMAX,NUSEMAX
      PARAMETER (NCFMAX=8,ISEQMAX=30,NUSEMAX=7)
C
C Dummy arguments
C
      INTEGER IEST,NUSE,NV
      REAL*8 XEST
      COMPLEX*16 DY(NV),YEST(NV),YZ(NV)
C
C Local variables
C
      COMPLEX*16 B,B1,C,D(NCFMAX,NUSEMAX),DDY,V,YY
      REAL*8 FX(NUSEMAX),X(ISEQMAX)
      INTEGER J,K,M1
      SAVE B,B1,C,D,DDY,FX,J,K,M1,V,X,YY
C
      X(IEST) = XEST
      IF ( IEST.EQ.1 ) THEN
         DO J = 1,NV
            YZ(J) = YEST(J)
            D(J,1) = YEST(J)
            DY(J) = YEST(J)
         END DO
      ELSE
         M1 = MIN(IEST,NUSE)
         DO K = 1,M1 - 1
            FX(K+1) = X(IEST-K)/XEST
         END DO
         DO J = 1,NV
            YY = YEST(J)
            V = D(J,1)
            C = YY
            D(J,1) = YY
            DO K = 2,M1
               B1 = FX(K)*V
               B = B1 - C
               IF ( B.NE.0. ) THEN
                  B = (C-V)/B
                  DDY = C*B
                  C = B1*B
               ELSE
                  DDY = V
               END IF
               V = D(J,K)
               D(J,K) = DDY
               YY = YY + DDY
            END DO
            DY(J) = DDY
            YZ(J) = YY
         END DO
      END IF
      END
