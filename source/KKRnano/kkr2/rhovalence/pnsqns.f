      SUBROUTINE PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,QNS,NSRA,
     +                  VINS,IPAN,IRCUT,CLEB,ICLEB,IEND,LOFLM,LKONV,
     &                  ISPIN,LDAU,NLDAU,LLDAU,
     &                  WMLDAU,WMLDAUAV,LDAUCUT,
C                       new after inc.p removal
     &                  lmaxd, nspind, irmd, irnsd, ipand, ncleb)

      IMPLICIT NONE

      INTEGER lmaxd
      INTEGER nspind
      INTEGER irmd
      INTEGER irnsd
      INTEGER ipand
      INTEGER ncleb

C     .. Parameters ..
C     INTEGER             IRMIND
C     PARAMETER          (IRMIND=IRMD-IRNSD)
C     INTEGER             LMMAXD
C     PARAMETER          (LMMAXD= (LMAXD+1)**2)
C     INTEGER             LMPOTD
C     PARAMETER          (LMPOTD= (LPOTD+1)**2)
C                               = (2*LMAXD+1)**2
C     INTEGER             LMAXD1
C     PARAMETER          (LMAXD1= LMAXD+1)
C     INTEGER             MMAXD
C     PARAMETER          (MMAXD = 2*LMAXD+1 )
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX     EK
      INTEGER            ICST,IEND,IPAN,LKONV,NSRA,NLDAU,ISPIN
      LOGICAL            LDAU
C     ..
C     .. Array Arguments ..


C     DOUBLE COMPLEX     PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
C     DOUBLE COMPLEX     PZ(IRMD,0:LMAXD)
C     DOUBLE COMPLEX     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
C     DOUBLE COMPLEX     QZ(IRMD,0:LMAXD)
C     DOUBLE COMPLEX     SZ(IRMD,0:LMAXD)

C     DOUBLE COMPLEX     AR(LMMAXD,LMMAXD)
C     DOUBLE COMPLEX     CR(LMMAXD,LMMAXD)
      DOUBLE COMPLEX     AR((LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX     CR((LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX     FZ(IRMD,0:LMAXD)

C     DOUBLE COMPLEX     PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX     PNS((LMAXD+1)**2,(LMAXD+1)**2,
     &                       (IRMD-IRNSD):IRMD,2)

      DOUBLE COMPLEX     PZ(IRMD,0:LMAXD)
C     DOUBLE COMPLEX     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX     QNS((LMAXD+1)**2,(LMAXD+1)**2,
     &                       (IRMD-IRNSD):IRMD,2)

      DOUBLE COMPLEX     QZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX     SZ(IRMD,0:LMAXD)

      DOUBLE PRECISION   CLEB(NCLEB,2)
      DOUBLE PRECISION   DRDI(IRMD)
C     DOUBLE PRECISION   VINS(IRMIND:IRMD,LMPOTD)
      DOUBLE PRECISION   VINS(IRMD-IRNSD:IRMD,(2*LMAXD+1)**2)
C     DOUBLE PRECISION   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      DOUBLE PRECISION   WMLDAU(2*LMAXD+1,2*LMAXD+1,LMAXD+1,NSPIND)
      DOUBLE PRECISION   LDAUCUT(IRMD)

      INTEGER            ICLEB(NCLEB,3),IRCUT(0:IPAND),LOFLM(*)
      INTEGER            LLDAU(LMAXD+1)
C     ..
C     .. Local Scalars ..
      INTEGER            I,IR,IRC1,LM1,LM2,LMMKONV,M1,M2,
     +                   LMLO,LMHI,ILDAU
C     ..
C     .. Local Arrays ..
C     DOUBLE COMPLEX     CMAT(LMMAXD,LMMAXD,IRMIND:IRMD)
      DOUBLE COMPLEX     CMAT((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD)
C     DOUBLE COMPLEX     DMAT(LMMAXD,LMMAXD,IRMIND:IRMD)
      DOUBLE COMPLEX     DMAT((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD)
C     DOUBLE COMPLEX     DR(LMMAXD,LMMAXD)
      DOUBLE COMPLEX     DR((LMAXD+1)**2, (LMAXD+1)**2)
C     DOUBLE COMPLEX     EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX     EFAC((LMAXD+1)**2)
      DOUBLE COMPLEX     PZEKDR((LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
C     DOUBLE COMPLEX     PZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX     PZLM((LMAXD+1)**2,(IRMD-IRNSD):IRMD,2)
C     DOUBLE COMPLEX     QZEKDR(LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX     QZEKDR((LMAXD+1)**2,(IRMD-IRNSD):IRMD,2)
C     DOUBLE COMPLEX     QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX     QZLM((LMAXD+1)**2, (IRMD-IRNSD):IRMD,2)
C     DOUBLE COMPLEX     TMATLL(LMMAXD,LMMAXD)
      DOUBLE COMPLEX     TMATLL((LMAXD+1)**2, (LMAXD+1)**2)
C     DOUBLE PRECISION   VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      DOUBLE PRECISION   VNSPLL((LMAXD+1)**2, (LMAXD+1)**2,
     &                          IRMD-IRNSD:IRMD)

      DOUBLE PRECISION   WMLDAUAV(LMAXD+1)
C     ..
C     .. External Subroutines ..
      EXTERNAL IRWNS,REGNS,VLLNS,WFTSCA

      INTEGER             LMMAXD
      INTEGER             IRMIND

      IRMIND= IRMD-IRNSD
      LMMAXD= (LMAXD+1)**2

      IRC1 = IRCUT(IPAN)
c

      CALL VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,
     &           lmaxd, irmd, irnsd, ncleb)


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

C-----------------------------------------------------------------------
C LDA+U
C Add WLDAU to non-spherical porential VINS in case of LDA+U
C Use the average wldau (=wldauav) and the deviation
C of wldau from this. Use the deviation in the Born series
C for the non-spherical wavefunction, while the average is
C used for the spherical wavefunction.
C
      IF (LDAU) THEN
      DO ILDAU=1,NLDAU
      IF (LLDAU(ILDAU).GE.0) THEN
C
        LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
        LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)
C
        DO IR = IRMIND,IRMD
C
C -> First add wldau to all elements ...
C
          DO LM2 = LMLO,LMHI
            M2 = LM2 - LMLO + 1
            DO LM1 = LMLO,LMHI
              M1 = LM1 - LMLO + 1
              VNSPLL(LM1,LM2,IR) =VNSPLL(LM1,LM2,IR)
     &                           +WMLDAU(M1,M2,ILDAU,ISPIN)*LDAUCUT(IR)
            ENDDO
C
C ... and then subtract average from diag. elements
C
            VNSPLL(LM2,LM2,IR) =  VNSPLL(LM2,LM2,IR)
     &                         - WMLDAUAV(ILDAU) * LDAUCUT(IR)
          ENDDO
        ENDDO
      ENDIF
      ENDDO
      ENDIF
C
C LDA+U
C-----------------------------------------------------------------------


c
c---> get wfts of same magnitude by scaling with efac
c
      CALL WFTSCA(DRDI,EFAC,PZ,QZ,FZ,SZ,NSRA,PZLM,QZLM,PZEKDR,QZEKDR,
     +              EK,LOFLM,IRMIND,IRMD,LMAXD,LMMAXD)
c
c---> determine the irregular non sph. wavefunction
c
      CALL IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IPAN,IRCUT,NSRA,PZLM,QZLM,
     +             PZEKDR,QZEKDR,QNS(1,1,IRMIND,1),CMAT,
     +             QNS(1,1,IRMIND,2),DMAT,IRMIND,IRMD,IPAND,LMMAXD)
c
c---> determine the regular non sph. wavefunction
c
      CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,QZLM,
     +             PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,
     +             PNS(1,1,IRMIND,2),DMAT,NSRA,IRMIND,IRMD,IPAND,LMMAXD)
c

      RETURN

      END
