c ************************************************************************
      SUBROUTINE CONVOL(IMT1,IRC1,ICELL,IMAXSH,ILM,IFUNM,LMPOT,GSH,
     +                  THETAS,Z,RFPI,R,VONS,LMSP,
C                       new parameters after inc.p removal
     &                  irid, nfund, irmd, ngshd)
c ************************************************************************
      IMPLICIT NONE
C     .. Parameters ..
C      include 'inc.p'
C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)

C     Number of radial mesh points for shape functions
      INTEGER irid
C     maximal number of shape functions
      INTEGER nfund
      INTEGER irmd
      INTEGER ngshd
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION RFPI,Z
      INTEGER ICELL,IMAXSH,IMT1,IRC1,LMPOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GSH(*),R(*),THETAS(IRID,NFUND,*),VONS(IRMD,*)
      INTEGER IFUNM(*),ILM(NGSHD,3),LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ZZOR
      INTEGER I,IFUN,IR,IRH,LM,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION VSTORE(IRID,LMPOT)
C     ..
      DO 20 LM = 1,LMPOT
        DO 10 IR = 1,IRC1 - IMT1
          VSTORE(IR,LM) = 0.0D0
   10   CONTINUE
   20 CONTINUE

      DO 30 IR = IMT1 + 1,IRC1
        ZZOR = 2.0D0*Z/R(IR)*RFPI
        VONS(IR,1) = VONS(IR,1) - ZZOR
   30 CONTINUE

      DO 50 I = 1,IMAXSH
        LM1 = ILM(I,1)
        LM2 = ILM(I,2)
        LM3 = ILM(I,3)
        IF (LMSP(LM3).GT.0) THEN
          IFUN = IFUNM(LM3)
          DO 40 IR = IMT1 + 1,IRC1
            IRH = IR - IMT1
            VSTORE(IRH,LM1) = VSTORE(IRH,LM1) +
     +                        GSH(I)*VONS(IR,LM2)*THETAS(IRH,IFUN,ICELL)
   40     CONTINUE
        END IF
   50 CONTINUE

      DO 60 IR = IMT1 + 1,IRC1
        IRH = IR - IMT1
        ZZOR = 2.0D0*Z/R(IR)*RFPI
        VONS(IR,1) = VSTORE(IRH,1) + ZZOR
   60 CONTINUE

      DO 80 LM = 2,LMPOT
        DO 70 IR = IMT1 + 1,IRC1
          IRH = IR - IMT1
          VONS(IR,LM) = VSTORE(IRH,LM)
   70   CONTINUE
   80 CONTINUE

      RETURN

      END
