      MODULE mod_pnstmat
      CONTAINS
      SUBROUTINE PNSTMAT(DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,TMATLL,VINS,IPAN,
     +                   IRCUT,NSRA,CLEB,ICLEB,IEND,LOFLM,NCLEB,TMAT,
     +                   LKONV,IDOLDAU,LOPT,LMLO,LMHI,WLDAU,WLDAUAV,
     +                   CUTOFF,NRMAX, NRMIN_NS,
     +                   LMAXATOM,MMAXATOM,LMMAXATOM,LMPOTATOM)
      USE mod_vllns
      USE mod_wftsca
      USE mod_regns

      IMPLICIT NONE
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
C *  LDA+U implementation     Mar. 2002-Dec.2004                      *
C *                           ph.mavropoulos, h. ebert, v. popescu    *
C *                                                                   *
C *********************************************************************
C
!       INTEGER NRMIN_NS
!       PARAMETER (NRMIN_NS=NRMAX-IRNSD)
!       INTEGER MMAXATOM
!       PARAMETER (MMAXATOM=2*LMAXATOM+1)
!       INTEGER LMMAXATOM
!       PARAMETER (LMMAXATOM= (KREL+1) * (LMAXATOM+1)**2)
!       INTEGER LMPOTATOM
!       PARAMETER (LMPOTATOM= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER NRMAX, NRMIN_NS,LMAXATOM,MMAXATOM,LMMAXATOM,
     +        LMPOTATOM,NCLEB
      DOUBLE COMPLEX EK
      INTEGER ICST,IDOLDAU,IEND,IPAN,LKONV,LOPT,NSRA,LMLO,LMHI
      DOUBLE PRECISION WLDAUAV
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX FZ(NRMAX,0:LMAXATOM),
     +               PNS(LMMAXATOM,LMMAXATOM,NRMIN_NS:NRMAX,2),
     +               PZ(NRMAX,0:LMAXATOM),QZ(NRMAX,0:LMAXATOM),
     +               SZ(NRMAX,0:LMAXATOM),
     +               TMAT(0:LMAXATOM),TMATLL(LMMAXATOM,LMMAXATOM)
      DOUBLE PRECISION CLEB(NCLEB,2),DRDI(NRMAX),
     +                 VINS(NRMIN_NS:NRMAX,LMPOTATOM)
      DOUBLE PRECISION WLDAU(MMAXATOM,MMAXATOM),CUTOFF(NRMAX)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAN),LOFLM(:)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,LM1,LM2,LMMKONV,M1,M2
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX AR(LMMAXATOM,LMMAXATOM),CMAT(LMMAXATOM,
     +               LMMAXATOM,NRMIN_NS:NRMAX),
     +               DMAT(LMMAXATOM,LMMAXATOM,NRMIN_NS:NRMAX),
     +               EFAC(LMMAXATOM),
     +               PZEKDR(LMMAXATOM,NRMIN_NS:NRMAX,2),
     +               PZLM(LMMAXATOM,NRMIN_NS:NRMAX,2),
     +               QZEKDR(LMMAXATOM,NRMIN_NS:NRMAX,2),
     +               QZLM(LMMAXATOM,NRMIN_NS:NRMAX,2)
      DOUBLE PRECISION VNSPLL(LMMAXATOM,LMMAXATOM,NRMIN_NS:NRMAX)
C     ..
C     .. External Subroutines ..
!       EXTERNAL REGNS,VLLNS,WFTSCA
C     ..
!       DO LM1=1,LMPOTATOM
!            write(*,*) 'lm ',lm1
!            write(*,*) 'lm ',vins(:,lm1)
!       END DO


      CALL VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,
     +           NCLEB,LMAXATOM,LMMAXATOM,LMPOTATOM,NRMIN_NS,NRMAX )
      IF (LKONV.NE.LMAXATOM) THEN
        LMMKONV = (LKONV+1)* (LKONV+1)
        DO LM1 = 1,LMMAXATOM
          DO LM2 = LMMKONV + 1,LMMAXATOM
            DO I = NRMIN_NS,NRMAX
              VNSPLL(LM2,LM1,I) = 0.0D0
              VNSPLL(LM1,LM2,I) = 0.0D0
            END DO
          END DO
        END DO
      ELSE
        LMMKONV = LMMAXATOM
      END IF
C ======================================================================
C LDA+U
C Add WLDAU to non-spherical porential VINS in case of LDA+U
C Use the average wldau (=wldauav) and calculate the deviation
C of wldau from this. Use the deviation in the Born series
C for the non-spherical wavefunction, while the average is
C used for the spherical wavefunction.
C
      IF ( IDOLDAU.EQ.1.AND.LOPT.GE.0 ) THEN
         DO IR = NRMIN_NS,NRMAX
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
     +              EK,LOFLM,NRMIN_NS,NRMAX,LMAXATOM,LMMAXATOM)
c
c---> determine the regular non sph. wavefunction
c
      CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,QZLM,
     +             PZEKDR,QZEKDR,EK,PNS(1,1,NRMIN_NS,1),CMAT,
     +             PNS(1,1,NRMIN_NS,2),DMAT,NSRA,NRMIN_NS,NRMAX,IPAN,
     +             LMMAXATOM)
c

      DO LM1 = 1,LMMKONV
        TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1))
      END DO

      RETURN
      END SUBROUTINE
      END MODULE mod_pnstmat
