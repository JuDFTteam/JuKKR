c ************************************************************************
      SUBROUTINE GFREE(RDIFF,E0,GMLL,CLEB,ICLEB,LOFLM,IEND)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMAX
c      PARAMETER (LMAX=4)
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER LMAX2P,LMX2SQ
      PARAMETER (LMAX2P=LMAX*2+1,LMX2SQ=LMAX2P**2)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMX2SQ*LMAXSQ)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E0
      INTEGER IEND
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GMLL(LMAXSQ,LMAXSQ)
      DOUBLE PRECISION CLEB(NCLEB),RDIFF(*)
      INTEGER ICLEB(NCLEB,4),LOFLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FPI,PI,RABS,RFPI,X,Y,Z
      INTEGER IFAC,J,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BL(14),HL(14),HYL(LMX2SQ),NL(14)
      DOUBLE PRECISION YL(LMX2SQ)
      INTEGER LF(LMX2SQ)
C     ..
C     .. External Subroutines ..
      EXTERNAL BESHAN,YSHYSH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DSQRT,ZSQRT
C     ..
      PI = 4.D0*DATAN(1.D0)
      FPI = 4.D0*PI
      RFPI = DSQRT(FPI)
C-----------------------------------------------------------------------
C---- CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
C-----------------------------------------------------------------------
      DO 10 LM1 = 1,LMX2SQ
        LF(LM1) = LOFLM(LM1) + 1
   10 CONTINUE
      X = RDIFF(1)
      Y = RDIFF(2)
      Z = RDIFF(3)
      CALL YSHYSH(X,Y,Z,RABS,YL)
      CALL BESHAN(HL,BL,NL,ZSQRT(E0)*RABS,0,LMAX*2)
      DO 20 LM1 = 1,LMX2SQ
        HYL(LM1) = -FPI*CI*ZSQRT(E0)*YL(LM1)*HL(LF(LM1))
   20 CONTINUE
      DO 40 LM1 = 1,LMAXSQ
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
      DO 70 LM1 = 1,LMAXSQ
        DO 60 LM2 = 1,LM1 - 1
          IFAC = (-1)** (LOFLM(LM1)+LOFLM(LM2))
          GMLL(LM2,LM1) = IFAC*GMLL(LM1,LM2)
   60   CONTINUE
   70 CONTINUE
      RETURN

      END
