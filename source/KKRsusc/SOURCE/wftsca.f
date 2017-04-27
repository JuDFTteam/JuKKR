      MODULE mod_WFTSCA
      CONTAINS
      SUBROUTINE WFTSCA(DRDI,EFAC,PZ,QZ,FZ,SZ,NSRA,PZLM,QZLM,PZEKDR,
     +                  QZEKDR,EK,LOFLM,IRMIND,IRMD,LMAXATOM,LMMAXATOM)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c                 R. Zeller      Oct. 1993
c-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EK
      INTEGER IRMD,IRMIND,LMAXATOM,LMMAXATOM,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX EFAC(LMMAXATOM),FZ(IRMD,0:LMAXATOM),
     +               PZ(IRMD,0:LMAXATOM),
     +               PZEKDR(LMMAXATOM,IRMIND:IRMD,2),
     +               PZLM(LMMAXATOM,IRMIND:IRMD,2),QZ(IRMD,0:LMAXATOM),
     +               QZEKDR(LMMAXATOM,IRMIND:IRMD,2),
     +               QZLM(LMMAXATOM,IRMIND:IRMD,2),SZ(IRMD,0:LMAXATOM)
      DOUBLE PRECISION DRDI(:)
      INTEGER LOFLM(:)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EFAC1,V1
      INTEGER IR,J,L,L1,LM,LM1,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..

c
c---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!
c
      EFAC(1) = CONE
      V1 = CONE
      DO 20 L = 1,LMAXATOM
        V1 = V1*EK/DBLE(2*L-1)
        DO 10 M = -L,L
          LM = L* (L+1) + M + 1
          EFAC(LM) = V1
   10   CONTINUE
   20 CONTINUE
c
c
c---> get wfts of same magnitude by scaling with efac
c
      DO 70 LM1 = 1,LMMAXATOM
        L1 = LOFLM(LM1)
        EFAC1 = EFAC(LM1)
        DO 30 IR = IRMIND,IRMD
          PZLM(LM1,IR,1) = PZ(IR,L1)/EFAC1
          QZLM(LM1,IR,1) = QZ(IR,L1)*EFAC1
   30   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 40 IR = IRMIND,IRMD
            PZLM(LM1,IR,NSRA) = FZ(IR,L1)/EFAC1
            QZLM(LM1,IR,NSRA) = SZ(IR,L1)*EFAC1
   40     CONTINUE
        END IF
c
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIND,IRMD
            PZEKDR(LM1,IR,J) = PZLM(LM1,IR,J)*EK*DRDI(IR)
            QZEKDR(LM1,IR,J) = QZLM(LM1,IR,J)*EK*DRDI(IR)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE


      END SUBROUTINE
      END MODULE mod_WFTSCA
