      MODULE mod_VLLNS
      CONTAINS
      SUBROUTINE VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,
     + NCLEB,LMAXATOM,LMMAXATOM,LMPOTATOM,IRMIND,IRMD )
      IMPLICIT NONE
c-----------------------------------------------------------------------
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
C     .. Parameters ..
!       INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXATOM = 2 * (LMAXATOM+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
      INTEGER NCLEB,LMAXATOM,LMMAXATOM,LMPOTATOM,IRMIND,IRMD 
!       INTEGER IRMIND
!       PARAMETER (IRMIND=IRMD-IRNSD)
!       INTEGER LMMAXATOM
!       PARAMETER (LMMAXATOM= (krel+1) * (LMAXATOM+1)**2)
!       INTEGER LMPOTATOM
!       PARAMETER (LMPOTATOM= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(NCLEB,2),VINS(IRMIND:IRMD,LMPOTATOM),
     +                 VNSPLL(LMMAXATOM,LMMAXATOM,IRMIND:IRMD)
      INTEGER ICLEB(NCLEB,4)
C     ..
C     .. Local Scalars ..
      INTEGER IR,J,LM1,LM2,LM3
C     ..
!       DO 30 LM1 = 1,LMMAXATOM
!         DO 20 LM2 = 1,LM1
!           DO 10 IR = IRMIND,IRMD
!             VNSPLL(LM1,LM2,IR) = 0.0D0
!    10     CONTINUE
!    20   CONTINUE
!    30 CONTINUE
      VNSPLL(:,:,:)=0.0D0
c
      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        IF (LM1 <= LMMAXATOM .AND. LM2 <= LMMAXATOM .AND.
     +      LM3 >  1) THEN
          DO 40 IR = IRMIND,IRMD
            VNSPLL(LM1,LM2,IR) = VNSPLL(LM1,LM2,IR) +
     +                           CLEB(J,1)*VINS(IR,LM3)
   40     CONTINUE
        END IF
   50 CONTINUE
c
c---> use symmetry of the gaunt coef.
c
      DO 80 LM1 = 1,LMMAXATOM
        DO 70 LM2 = 1,LM1 - 1
            DO 60 IR = IRMIND,IRMD
              VNSPLL(LM2,LM1,IR) = VNSPLL(LM1,LM2,IR)
   60       CONTINUE
   70   CONTINUE
   80 CONTINUE

c    comment out by D. Bauer 10/12/09
c     DO LM1 = 1,LMMAXATOM
c       DO IR = IRMIND,IRMD
c         VNSPLL(LM1,LM1,IR) = VNSPLL(LM1,LM1,IR) + VINS(IR,1)
c       END DO
c     END DO

      END SUBROUTINE
      END MODULE mod_VLLNS
