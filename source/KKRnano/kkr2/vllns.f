      SUBROUTINE VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,
C                      new input parameters after inc.p removal
     &                 lmax, irmd, irnsd, ncleb)
      IMPLICIT NONE
c-----------------------------------------------------------------------
C     Calculates V_LL' from V_L.
c     to determine the non - spherical wavefunctions the potential
c         has to be lm1 and lm2 dependent . the potential is stored
c         only as lm dependent , therefore a transformation in the
c         following way has to be done :
c
c        vnsll(r,lm1,lm2)   =   {  c(lm1,lm2,lm3) *vins(r,lm3)  }
c                                  (summed over lm3 at the right site )
c        where c(lm1,lm2,lm3) are the gaunt coeffients .
c
c             (see notes by b.drittler)
c
c     attention : the gaunt coeffients are stored in an index array
c                  only for lm1.gt.lm2
c                 (see subroutine gaunt)
c
c                               b.drittler   july 1988
c-----------------------------------------------------------------------
c                          modified by R. Zeller Sep. 2000
c-----------------------------------------------------------------------

      INTEGER lmax
      INTEGER irmd
      INTEGER irnsd
      INTEGER ncleb

C     INTEGER IRMIND
C     PARAMETER (IRMIND=IRMD-IRNSD)
C     INTEGER LMMAXD
C     PARAMETER (LMMAXD= (LMAXD+1)**2)
C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION CLEB(NCLEB,2),VINS(IRMIND:IRMD,LMPOTD),
C    +                 VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)

      DOUBLE PRECISION CLEB(NCLEB,2)
      DOUBLE PRECISION VINS(IRMD-IRNSD:IRMD,(2*LMAX+1)**2)
      DOUBLE PRECISION VNSPLL((LMAX+1)**2, (LMAX+1)**2,IRMD-IRNSD:IRMD)

      INTEGER ICLEB(NCLEB,3)
C     ..
C     .. Local Scalars ..
      INTEGER IR,J,LM1,LM2,LM3
C     ..
      INTEGER LMMAXD
      INTEGER IRMIND
      IRMIND=IRMD-IRNSD
      LMMAXD= (LMAX+1)**2

      DO 30 LM1 = 1,LMMAXD
        DO 20 LM2 = 1,LM1
          DO 10 IR = IRMIND,IRMD
            VNSPLL(LM1,LM2,IR) = 0.0D0
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
c
      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        DO 40 IR = IRMIND,IRMD
          VNSPLL(LM1,LM2,IR) = VNSPLL(LM1,LM2,IR) +
     +                         CLEB(J,1)*VINS(IR,LM3)
   40   CONTINUE
   50 CONTINUE

c
c---> use symmetry of the gaunt coef.
c
      DO 80 LM1 = 1,LMMAXD
        DO 70 LM2 = 1,LM1 - 1
          DO 60 IR = IRMIND,IRMD
            VNSPLL(LM2,LM1,IR) = VNSPLL(LM1,LM2,IR)
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE

      DO LM1 = 1,LMMAXD
        DO IR = IRMIND,IRMD
          VNSPLL(LM1,LM1,IR) = VNSPLL(LM1,LM1,IR) + VINS(IR,1)
        END DO
      END DO

      END
