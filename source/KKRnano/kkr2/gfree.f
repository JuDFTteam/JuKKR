c ************************************************************************
      SUBROUTINE GFREE(RDIFF,E0,GMLL,CLEB,ICLEB,LOFLM,IEND,
C                      new input parameters after inc.p removal
     &                 lmax, ncleb)
c ************************************************************************
      IMPLICIT NONE

      INTEGER LMAX
      INTEGER ncleb

C      INTEGER LMGF0D
C      PARAMETER (LMGF0D= (LMAX+1)**2)
C      INTEGER LMAX2P,LMX2SQ
C      PARAMETER (LMAX2P=LMAX*2+1,LMX2SQ=LMAX2P**2)
C                                       =(LMAX*2+1)**2
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E0
      INTEGER IEND
C     ..
C     .. Array Arguments ..
C     DOUBLE COMPLEX GMLL(LMGF0D,LMGF0D)
      DOUBLE COMPLEX GMLL((LMAX +1 )**2, (LMAX + 1)**2)
      DOUBLE PRECISION CLEB(NCLEB),RDIFF(*)
      INTEGER ICLEB(NCLEB,3),LOFLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FPI,PI,RABS,RFPI,X,Y,Z
      INTEGER IFAC,J,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
C     DOUBLE COMPLEX BL(LMAX2P),HL(LMAX2P),HYL(LMX2SQ),NL(LMAX2P)
C     DOUBLE PRECISION YL(LMX2SQ)
C     INTEGER LF(LMX2SQ)

      DOUBLE COMPLEX BL(LMAX*2+1),HL(LMAX*2+1)
      DOUBLE COMPLEX HYL((LMAX*2+1)**2),NL(LMAX*2+1)
      DOUBLE PRECISION YL((LMAX*2+1)**2)
      INTEGER LF((LMAX*2+1)**2)

C     ..
C     .. External Subroutines ..
      EXTERNAL BESHAN,YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT

      INTEGER LMGF0D
      INTEGER LMX2SQ

      LMGF0D = (LMAX+1)**2
      LMX2SQ = (LMAX*2+1)**2

C     ..
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      RFPI = SQRT(FPI)
C-----------------------------------------------------------------------
C---- CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
C-----------------------------------------------------------------------
      DO 10 LM1 = 1,LMX2SQ
        LF(LM1) = LOFLM(LM1) + 1
   10 CONTINUE
      X = RDIFF(1)
      Y = RDIFF(2)
      Z = RDIFF(3)
      CALL YMY(X,Y,Z,RABS,YL,LMAX*2)
      CALL BESHAN(HL,BL,NL,SQRT(E0)*RABS,LMAX*2)
      DO 20 LM1 = 1,LMX2SQ
        HYL(LM1) = -FPI*CI*SQRT(E0)*YL(LM1)*HL(LF(LM1))
   20 CONTINUE
      DO 40 LM1 = 1,LMGF0D
        GMLL(LM1,LM1) = HYL(1)/RFPI
        DO 30 LM2 = 1,LM1 - 1
          GMLL(LM1,LM2) = CZERO
   30   CONTINUE
   40 CONTINUE
      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        GMLL(LM1,LM2) = GMLL(LM1,LM2) + CLEB(J)*HYL(LM3)
   50 CONTINUE
      DO 70 LM1 = 1,LMGF0D
        DO 60 LM2 = 1,LM1 - 1
          IFAC = (-1)** (LOFLM(LM1)+LOFLM(LM2))
          GMLL(LM2,LM1) = IFAC*GMLL(LM1,LM2)
   60   CONTINUE
   70 CONTINUE
      RETURN

      END
c ************************************************************************
      SUBROUTINE DGFREE(RDIFF,E0,DGMLL,CLEB,ICLEB,LOFLM,IEND,
C                       new input parameters after inc.p removal
     &                  lmax, ncleb)
c ************************************************************************
      IMPLICIT NONE

      INTEGER LMAX
      INTEGER ncleb

C     INTEGER LMGF0D
C     PARAMETER (LMGF0D= (LMAX+1)**2)
C     INTEGER LMAX2P,LMX2SQ
C     PARAMETER (LMAX2P=LMAX*2+1,LMX2SQ=LMAX2P**2)

      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E0
      INTEGER IEND
C     ..
C     .. Array Arguments ..
C     DOUBLE COMPLEX DGMLL(LMGF0D,LMGF0D)
      DOUBLE COMPLEX DGMLL((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE PRECISION CLEB(NCLEB),RDIFF(*)
      INTEGER ICLEB(NCLEB,3),LOFLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FPI,PI,RABS,RFPI,X,Y,Z
      INTEGER IFAC,J,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
C     DOUBLE COMPLEX BL(LMAX2P),HL(LMAX2P),NL(LMAX2P)
C     DOUBLE COMPLEX DHL(LMAX2P),DHYL(LMX2SQ)
C     DOUBLE PRECISION YL(LMX2SQ)
C     INTEGER LF(LMX2SQ)

      DOUBLE COMPLEX BL(LMAX*2+1),HL(LMAX*2+1),NL(LMAX*2+1)
      DOUBLE COMPLEX DHL(LMAX*2+1),DHYL((LMAX*2+1)**2)
      DOUBLE PRECISION YL((LMAX*2+1)**2)
      INTEGER LF((LMAX*2+1)**2)


C     ..
C     .. External Subroutines ..
      EXTERNAL BESHAN,YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT

      INTEGER LMGF0D
      INTEGER LMAX2P
      INTEGER LMX2SQ
      INTEGER LP1

      LMGF0D= (LMAX+1)**2
      LMAX2P=  LMAX*2+1
      LMX2SQ= (LMAX*2+1)**2

C     ..
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      RFPI = SQRT(FPI)
C-----------------------------------------------------------------------
C---- derivative of free electron green function matrix elements
C
C     the analytical formula for the derivative of spherical Hankel
C     functions is used:
C
C     d                     l+1
C     --  h (x) = h   (x) - --- h (x)   
C     dx   l       l-1       x   l
C
C     which for x = sqrt(E0)*r leads to
C
C      d                       r           rl
C     --- ( sqrt(E0) h (x) ) = - h   (x) - -- h (x) )
C     dE0             l        2  l-1      2x  l
C
C-----------------------------------------------------------------------
      DO 10 LM1 = 1,LMX2SQ
        LF(LM1) = LOFLM(LM1) + 1
   10 CONTINUE
      X = RDIFF(1)
      Y = RDIFF(2)
      Z = RDIFF(3)
      CALL YMY(X,Y,Z,RABS,YL,LMAX*2)
      CALL BESHAN(HL,BL,NL,SQRT(E0)*RABS,LMAX*2)
      DHL(1) = 0.5D0*CI*RABS*HL(1)
      DO LP1 = 2,LMAX2P
      DHL(LP1) = 0.5D0*(RABS*HL(LP1-1)-(LP1-1)*HL(LP1)/SQRT(E0))
      END DO 
      DO 20 LM1 = 1,LMX2SQ
        DHYL(LM1) = -FPI*CI*YL(LM1)*DHL(LF(LM1))
   20 CONTINUE
      DO 40 LM1 = 1,LMGF0D
        DGMLL(LM1,LM1) = DHYL(1)/RFPI
        DO 30 LM2 = 1,LM1 - 1
          DGMLL(LM1,LM2) = CZERO
   30   CONTINUE
   40 CONTINUE
      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        DGMLL(LM1,LM2) = DGMLL(LM1,LM2) + CLEB(J)*DHYL(LM3)
   50 CONTINUE
      DO 70 LM1 = 1,LMGF0D
        DO 60 LM2 = 1,LM1 - 1
          IFAC = (-1)** (LOFLM(LM1)+LOFLM(LM2))
          DGMLL(LM2,LM1) = IFAC*DGMLL(LM1,LM2)
   60   CONTINUE
   70 CONTINUE
      RETURN

      END
