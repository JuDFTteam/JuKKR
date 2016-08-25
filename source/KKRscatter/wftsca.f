C ************************************************************************
      SUBROUTINE WFTSCA(DRDI,EFAC,LMAX,PZ,QZ,IRMIN,IRWS,IPAN,IRCUT,FZ,
     +                  SZ,NSRA,PZLM,QZLM,PZEKDR,QZEKDR,EK,LOFLM)
C ************************************************************************
c                 R. Zeller      Oct. 1993
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EK
      INTEGER IPAN,IRMIN,IRWS,LMAX,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX EFAC(LMMAXD),FZ(IRMD,0:LMAXD),PZ(IRMD,0:LMAXD),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2),SZ(IRMD,0:LMAXD)
      DOUBLE PRECISION DRDI(*)
      INTEGER IRCUT(0:IPAND),LOFLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EFAC1,V1
      INTEGER IR,IRC1,IRS1,J,L,L1,LM,LM1,LMMAX,M
c
c
      LOGICAL TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..

c
c---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!
c
      EFAC(1) = CONE
      V1 = CONE
      DO 20 L = 1,LMAX
        V1 = V1*EK/REAL(2*L-1)
        DO 10 M = -L,L
          LM = L* (L+1) + M + 1
          EFAC(LM) = V1
   10   CONTINUE
   20 CONTINUE

c
      LMMAX = (LMAX+1)* (LMAX+1)
      IRC1 = IRCUT(IPAN)
      IF (IPAN.EQ.1) THEN
        IRS1 = IRWS

      ELSE
        IRS1 = IRC1
      END IF
c
c---> get wfts of same magnitude by scaling with efac
c
      DO 70 LM1 = 1,LMMAX
        L1 = LOFLM(LM1)
        EFAC1 = EFAC(LM1)
        DO 30 IR = IRMIN,IRS1
          PZLM(LM1,IR,1) = PZ(IR,L1)/EFAC1
          QZLM(LM1,IR,1) = QZ(IR,L1)*EFAC1
   30   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 40 IR = IRMIN,IRS1
            PZLM(LM1,IR,NSRA) = FZ(IR,L1)/EFAC1
            QZLM(LM1,IR,NSRA) = SZ(IR,L1)*EFAC1
   40     CONTINUE
        END IF
c
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIN,IRS1
            PZEKDR(LM1,IR,J) = PZLM(LM1,IR,J)*EK*DRDI(IR)
            QZEKDR(LM1,IR,J) = QZLM(LM1,IR,J)*EK*DRDI(IR)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE


      END
