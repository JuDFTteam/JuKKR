      SUBROUTINE PNSTMAT(DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,TMATLL,VINS,
     +                   IPAN,IRCUT,NSRA,CLEB,ICLEB,IEND,LOFLM,TMAT,
     +                   DET,LKONV,ISPIN,
     >                   LDAU,NLDAU,LLDAU,
     >                   WMLDAU,WMLDAUAV,LDAUCUT)
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
C
      INTEGER             IRMIND
      PARAMETER          (IRMIND=IRMD-IRNSD)
      INTEGER             LMMAXD
      PARAMETER          (LMMAXD= (LMAXD+1)**2)
      INTEGER             LMPOTD
      PARAMETER          (LMPOTD= (LPOTD+1)**2)
      INTEGER             LMAXD1
      PARAMETER          (LMAXD1 = LMAXD + 1)
      INTEGER             MMAXD
      PARAMETER          (MMAXD=2*LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX     EK,DET
      INTEGER            ICST,IEND,IPAN,LKONV,NSRA,
     +                   NLDAU
      LOGICAL            LDAU
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX     FZ(IRMD,0:LMAXD),
     +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +                   PZ(IRMD,0:LMAXD),QZ(IRMD,0:LMAXD),
     +                   SZ(IRMD,0:LMAXD),
     +                   TMAT(0:LMAXD),TMATLL(LMMAXD,LMMAXD)
      DOUBLE PRECISION   CLEB(NCLEB,2),DRDI(IRMD),
     +                   VINS(IRMIND:IRMD,LMPOTD),
     +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1),
     +                   LDAUCUT(IRMD),
     +                   WMLDAUAV(LMAXD1)
      INTEGER            ICLEB(NCLEB,3),IRCUT(0:IPAND),LOFLM(*),
     +                   LLDAU(LMAXD1)
C     ..
C     .. Local Scalars ..
      INTEGER            I,IRC1,LM1,LM2,LMMKONV,ISPIN,IR,M1,M2,
     +                   LMLO,LMHI,ILDAU
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX     AR(LMMAXD,LMMAXD),
     +                   CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),
     +                   DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),EFAC(LMMAXD),
     +                   PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +                   PZLM(LMMAXD,IRMIND:IRMD,2),
     +                   QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +                   QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION   VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER            IPVT(LMMAXD),INFO
C     ..
C     .. External Subroutines ..
      EXTERNAL REGNS,VLLNS,WFTSCA
C     ..
      IRC1 = IRCUT(IPAN)
c
C      DO LM1 = 1,LMMAXD
C        DO LM2 = 1,LMMAXD
C          DO I = IRMIND,IRMD
C            VNSPLL(LM1,LM2,I) = 0.0D0
C          ENDDO   
C        ENDDO
C      ENDDO

      CALL VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND)

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
C Use the average wldau (=wldauav) and calculate the deviation
C of wldau from this. Use the deviation in the Born series
C for the non-spherical wavefunction, while the average is
C used for the spherical wavefunction.
C
      IF (LDAU) THEN
      DO ILDAU=1,NLDAU
      IF (LLDAU(ILDAU).GE.0) THEN

        LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
        LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)

        DO IR = IRMIND,IRMD
C
C -> First add wldau to all elements ...
C
          DO LM2 = LMLO,LMHI
            M2 = LM2 - LMLO + 1
            DO LM1 = LMLO,LMHI
              M1 = LM1 - LMLO + 1
              VNSPLL(LM1,LM2,IR) = VNSPLL(LM1,LM2,IR)
     &                           + WMLDAU(M1,M2,ISPIN,ILDAU)*LDAUCUT(IR)
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
c---> determine the regular non sph. wavefunction
c
      CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,QZLM,
     +             PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,
     +             PNS(1,1,IRMIND,2),DMAT,NSRA,IRMIND,IRMD,IPAND,LMMAXD)

c

      DO LM1 = 1,LMMKONV
        TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1))
      END DO
        DET = (1.D0,0.D0)
        CALL ZGETRF(LMMAXD,LMMAXD,AR,LMMAXD,IPVT,INFO)
        DO LM1 = 1,LMMAXD
          IF (IPVT(LM1).NE.LM1) DET = -DET
          DET = AR(LM1,LM1)*DET
        ENDDO

      RETURN
      END
