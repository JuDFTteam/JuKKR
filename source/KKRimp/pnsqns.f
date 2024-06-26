      MODULE MOD_PNSQNS
      CONTAINS
      SUBROUTINE PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,QNS,NSRA,
     +                  VINS,IPAN,IRCUT,CLEB,ICLEB,IEND,LOFLM,LKONV,
     +                  IDOLDAU,LOPT,LMLO,LMHI,WLDAU,WLDAUAV,CUTOFF,
     +                  IRMD,IRMIND,NCLEB,
     +                  LMAXD,MMAXD,LMMAXD,LMPOTD)
      USE MOD_VLLNS
      USE MOD_WFTSCA
      USE MOD_REGNS
      USE MOD_IRWNS
      IMPLICIT NONE
C     .. Parameters ..
!       INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMAXD,MMAXD,LMMAXD,LMPOTD
      INTEGER NCLEB
      INTEGER IRMD
      INTEGER IRMIND
!       PARAMETER (IRMIND=IRMD-IRNSD)
!       INTEGER MMAXD
! !       PARAMETER (MMAXD=2*LMAXD+1)
!       INTEGER LMMAXD
!       PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EK
      DOUBLE PRECISION WLDAUAV
      INTEGER ICST,IDOLDAU,IEND,IPAN,LKONV,LMLO,LMHI,LOPT,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +               FZ(IRMD,0:LMAXD),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               PZ(IRMD,0:LMAXD),QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
      DOUBLE PRECISION CLEB(NCLEB,2),DRDI(IRMD),VINS(IRMIND:IRMD,LMPOTD)
      DOUBLE PRECISION CUTOFF(IRMD),WLDAU(MMAXD,MMAXD)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAN),LOFLM(:)
C     ..
C     .. Local Scalars ..
      INTEGER I,LM1,LM2,LMMKONV,M1,M2,IR
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),DR(LMMAXD,LMMAXD),
     +               EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2),TMATLL(LMMAXD,LMMAXD)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. External Subroutines ..
!       EXTERNAL IRWNS
C     ..
      CALL VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,
     +           NCLEB,LMAXD,LMMAXD,LMPOTD,IRMIND,IRMD )

      IF (LKONV.NE.LMAXD) THEN
        LMMKONV = (LKONV+1)* (LKONV+1)
        DO LM1 = 1,LMMAXD
          DO LM2 = LMMKONV + 1,LMMAXD
            DO I = IRMIND,IRMD
              VNSPLL(LM2,LM1,I) = 0.0D0
              VNSPLL(LM1,LM2,I) = 0.0D0
            END DO
          END DO
        END DO
      ELSE
        LMMKONV = LMMAXD
      END IF
C ======================================================================
C LDA+U
C Add WLDAU to non-spherical porential VINS in case of LDA+U
C Use the average wldau (=wldauav) and the deviation
C of wldau from this. Use the deviation in the Born series
C for the non-spherical wavefunction, while the average is
C used for the spherical wavefunction.
C
      IF ( IDOLDAU.EQ.1.AND.LOPT.GE.0 ) THEN
         DO IR = IRMIND,IRMD
C
C -> First add wldau to all elements ...
C
            DO LM2 = LMLO,LMHI
               M2 = LM2 - LMLO + 1
               DO LM1 = LMLO,LMHI
                  M1 = LM1 - LMLO + 1
                  VNSPLL(LM1,LM2,IR) =  VNSPLL(LM1,LM2,IR)
     &                            + WLDAU(M1,M2) * CUTOFF(IR)
               ENDDO
C
C ... and then subtract average from diag. elements
C
               VNSPLL(LM2,LM2,IR) =  VNSPLL(LM2,LM2,IR)
     &                            - WLDAUAV * CUTOFF(IR)
            ENDDO
         ENDDO
      END IF
C
C LDA+U
C ======================================================================
c
c---> get wfts of same magnitude by scaling with efac
c
      CALL WFTSCA(DRDI,EFAC,PZ,QZ,FZ,SZ,NSRA,PZLM,QZLM,PZEKDR,QZEKDR,
     +            EK,LOFLM,IRMIND,IRMD,LMAXD,LMMAXD)
c
c---> determine the irregular non sph. wavefunction
c
      CALL IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IPAN,IRCUT,NSRA,PZLM,QZLM,
     +             PZEKDR,QZEKDR,QNS(1,1,IRMIND,1),CMAT,
     +             QNS(1,1,IRMIND,2),DMAT,IRMIND,IRMD,IPAN,LMMAXD)
c
c---> determine the regular non sph. wavefunction
c
      CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,QZLM,
     +             PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,
     +             PNS(1,1,IRMIND,2),DMAT,NSRA,IRMIND,IRMD,IPAN,LMMAXD)
c

      RETURN

      END SUBROUTINE PNSQNS
      END MODULE MOD_PNSQNS
