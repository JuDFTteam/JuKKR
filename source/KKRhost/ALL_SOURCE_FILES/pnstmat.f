      SUBROUTINE PNSTMAT(DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,TMATLL,VINS,IRMIN,
     +                 IPAN,IRCUT,NSRA,CLEB,ICLEB,IEND,LOFLM,TMAT,LKONV,  ! Added IRMIN 1.7.2014
     +                   IDOLDAU,LOPT,LMLO,LMHI,WLDAU,WLDAUAV,CUTOFF,
     &                   ALPHA0)  ! LLY
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *  LDA+U implementation     Mar. 2002-Dec.2004                      *
C *                           ph.mavropoulos, h. ebert, v. popescu    *
C *                                                                   *
C *********************************************************************
C
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER MMAXD
      PARAMETER (MMAXD=2*LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EK
      INTEGER ICST,IDOLDAU,IEND,IPAN,LKONV,LOPT,NSRA,LMLO,LMHI,IRMIN
      DOUBLE PRECISION WLDAUAV
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX FZ(IRMD,0:LMAXD),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               PZ(IRMD,0:LMAXD),QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),
     +               TMAT(0:LMAXD),TMATLL(LMMAXD,LMMAXD),
     &               ALPHA0(LMMAXD,LMMAXD)  ! LLY
      DOUBLE PRECISION CLEB(NCLEB,2),DRDI(IRMD),VINS(IRMIND:IRMD,LMPOTD)
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD),CUTOFF(IRMD)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND),LOFLM(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,LM1,LM2,LMMKONV,M1,M2,IRMAX
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),EFAC(LMMAXD),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL REGNS,VLLNS,WFTSCA,ZGEMM
C     ..

      IRMAX = IRCUT(IPAN)

      call VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,
     + IRMD,NCLEB,LMPOTD,IRMIND,LMMAXD)
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
C Use the average wldau (=wldauav) and calculate the deviation
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
      PZLM(:,IRMIND:IRMD,:) = CZERO
      QZLM(:,IRMIND:IRMD,:) = CZERO
      PZEKDR(:,IRMIND:IRMD,:) = CZERO
      QZEKDR(:,IRMIND:IRMD,:) = CZERO
      CMAT(:,:,IRMIND:IRMD) = CZERO
      DMAT(:,:,IRMIND:IRMD) = CZERO
c
c---> get wfts of same magnitude by scaling with efac
c
      CALL WFTSCA(DRDI,EFAC,PZ,QZ,FZ,SZ,NSRA,PZLM,QZLM,PZEKDR,QZEKDR,
     +              EK,LOFLM,IRMIND,IRMD,IRMIN,IRMAX,LMAXD,LMMAXD)      ! Added IRMIN,IRMAX 1.7.2014
c
c---> determine the regular non sph. wavefunction
c
      CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,QZLM,
     +             PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,
     +             PNS(1,1,IRMIND,2),DMAT,NSRA,IRMIND,IRMD,IRMIN,IRMAX, ! Added IRMIN,IRMAX 1.7.2014
     &             IPAND,LMMAXD)
c

      DO LM1 = 1,LMMKONV
        TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1))
      END DO


      DO LM2 = 1,LMMKONV
         DO LM1 = 1,LMMKONV
            AR(LM1,LM2) = ALPHA0(LM1,LM1) * AR(LM1,LM2) ! LLY non-spher. contribution to alpha matrix
         ENDDO                                          ! LLY Drittler PhD eq. 3.106
      ENDDO
      ALPHA0(1:LMMAXD,1:LMMAXD) = AR(1:LMMAXD,1:LMMAXD) ! on output.


      RETURN
      END
