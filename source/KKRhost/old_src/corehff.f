      SUBROUTINE COREHFF(KAP1,KAP2,MJ,S,NSOL,BHF,GCK,FCK,RC,DRDIC,RNUC,
     &                   NZERO,NRC)
C
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE             *
C   *                  CURRENT  CORE STATE S                           *
C   *                                                                  *
C   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C
C PARAMETER definitions
C
      REAL*8 E0,A0,CAUTOG
C
C CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
C                                      ELECTRON CHARGE     IN ESU
C                                      BOHR-RADIUS         IN CM
C
C
      PARAMETER (E0=1.6021892D-19*2.997930D+09,A0=0.52917706D-08,
     &           CAUTOG=E0/(A0*A0))
C
C Dummy arguments
C
      INTEGER KAP1,KAP2,NRC,NSOL,NZERO,S
      REAL*8 MJ,RNUC
      REAL*8 BHF(2,2),DRDIC(NRC),FCK(2,2,NRC),GCK(2,2,NRC),RC(NRC)
C
C Local variables
C
      REAL*8 AME(2,2),XX(5),YI(NRC),YY(5),ZI(NRC)
      DOUBLE PRECISION DBLE,DSQRT
      INTEGER I,K1,K2,N
      REAL*8 YLAG
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
C        THE FACTOR  I  HAS BEEN OMITTED
C
      AME(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
      IF ( NSOL.EQ.2 ) THEN
         AME(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
         AME(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
         AME(2,1) = AME(1,2)
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
            DO N = 1,NZERO
               YI(N) = DRDIC(N)
     &                 *(GCK(K1,S,N)*FCK(K2,S,N)+FCK(K1,S,N)*GCK(K2,S,N)
     &                 )
            END DO
            CALL RINT4PTS(YI,NZERO,ZI)
            IF ( RNUC.NE.0.0D0 ) THEN
               DO I = 1,5
                  XX(I) = RC(NZERO-5+I)
                  YY(I) = ZI(NZERO-5+I)
               END DO
               ZI(NZERO) = YLAG(RNUC,XX,YY,0,4,5)
            END IF
            XX(1) = 1.0D0
C                      !RC( 1)
            XX(2) = 6.0D0
C                      !RC( 6)
            XX(3) = 11.0D0
C                      !RC(11)
            YY(1) = ZI(NZERO) - ZI(1)
            YY(2) = ZI(NZERO) - ZI(6)
            YY(3) = ZI(NZERO) - ZI(11)
            BHF(K1,K2) = CAUTOG*AME(K1,K2)*YLAG(0.0D0,XX,YY,0,2,3)
         END DO
      END DO
C
      END
