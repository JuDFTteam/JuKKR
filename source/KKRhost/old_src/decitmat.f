      SUBROUTINE DECITMAT(ERYD,ZAT,IPAN,RR,DROR,VISP,IRCUT,RIRC,
     &                    KREL,NSRA,INS,TMATLL,LOFLM,
     &                    IDOLDAU,LOPT,WLDAUAV,
     &                    SOLVER,SOCTL,CTL,ZREL,VTREL,BTREL,DRDI,R2DRDI,
     &                    IPAND,IRMD,LMAXD,LMAXDP1,LM2D,LMMAXD)
C **********************************************************************
C *                                                                    *
C * A modified form of the CALCTMAT routine to deal with the host      *
C * t-matrices in case of decimation                                   *
C *                                                                    *
C * Non-spherical potential not implemented yet, neither LDA+U         *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C ..
C ..  Parameters ..
      DOUBLE PRECISION CVLIGHT
      PARAMETER ( CVLIGHT=274.0720442D0 )
      DOUBLE COMPLEX CI
      PARAMETER ( CI = (0D0,1D0) )
C ..
C ..  Scalar arguments ..
      INTEGER IDOLDAU,IPAN,KREL,LOPT,NSRA,INS,ZREL
      INTEGER IPAND,IRMD,LM2D,LMAXD,LMAXDP1,LMMAXD
      DOUBLE PRECISION ZAT,RIRC,WLDAUAV
      DOUBLE COMPLEX ERYD
      CHARACTER*10 SOLVER
C ..
C ..  Array arguments ..
      INTEGER IRCUT(0:IPAND),LOFLM(LM2D)
      DOUBLE PRECISION RR(IRMD),DROR(IRMD),VISP(IRMD)
      DOUBLE COMPLEX TMATLL(LMMAXD,LMMAXD)
      DOUBLE PRECISION SOCTL(KREL*LMAXD+1)
      DOUBLE PRECISION CTL(KREL*LMAXD+1)
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION DRDI(IRMD),R2DRDI(IRMD*KREL+(1-KREL))
C ..
C ..  Local scalars ..
      INTEGER LL,LM1
      DOUBLE PRECISION RIRC1
      DOUBLE COMPLEX EK,CARG,QF,HLW,BLW
C ..
C ..  Local arrays ..
      DOUBLE PRECISION CUTOFF(IRMD)
      DOUBLE PRECISION RS(:,:),S(:)
      DOUBLE COMPLEX BESSJW(:),BESSYW(:),HANKWS(:),DLOGDP(:)
      DOUBLE COMPLEX TMAT(:),MASS(:),HAMF(:,:),FZ(:,:),PZ(:,:)
      ALLOCATABLE RS,S
      ALLOCATABLE BESSJW,BESSYW,HANKWS,DLOGDP
      ALLOCATABLE TMAT,MASS,HAMF,FZ,PZ
C ..
C ..  External subroutines ..
      EXTERNAL BESHAN,CINIT,REGSOL,WFMESH
C ..
      CALL CINIT(LMMAXD*LMMAXD,TMATLL)
C ================================================================= KREL
      IF ( KREL.EQ.0 ) THEN
         ALLOCATE (BESSJW(0:LMAXDP1),BESSYW(0:LMAXDP1),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Allocate BESSJW/BESSYW'
         ALLOCATE (HANKWS(0:LMAXDP1),DLOGDP(0:LMAXD),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Allocate HANKWS/DLOGFP'
         ALLOCATE (TMAT(0:LMAXD),MASS(IRMD),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Allocate TMAT/MASS'
         ALLOCATE (HAMF(IRMD,0:LMAXD),FZ(IRMD,0:LMAXD),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Allocate HAMF/FZ'
         ALLOCATE (PZ(IRMD,0:LMAXD),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Allocate PZ'
         ALLOCATE (RS(IRMD,0:LMAXD),S(0:LMAXD),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Allocate RS/S'
         RIRC1 = 1D0/RIRC
         CALL WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,RR,S,RS,IRCUT(IPAN),
     &               IRMD,LMAXD)
C     
         CARG = RIRC * EK
         CALL BESHAN(HANKWS,BESSJW,BESSYW,CARG,LMAXDP1)
         DO LL = 0,LMAXDP1
            HANKWS(LL) = BESSYW(LL) - CI*BESSJW(LL)
         END DO
C
         CALL REGSOL(CVLIGHT,ERYD,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,
     &               DROR,RR,S,VISP,ZAT,IPAN,IRCUT,IDOLDAU,LOPT,WLDAUAV,
     &               CUTOFF,IRMD,IPAND,LMAXD)
C
C ----------------------------------------------------------------------
C --> determine KREL=0 t - matrix
C
         DO LL = 0,LMAXD
            QF = DBLE(LL)*RIRC1
            HLW = HANKWS(LL) * DLOGDP(LL)
            BLW = BESSJW(LL) * DLOGDP(LL)
C     
            HLW = QF*HANKWS(LL) - EK*HANKWS(LL+1) - HLW
            BLW = BLW - QF*BESSJW(LL) + EK*BESSJW(LL+1)
            HLW = HLW * EK
            TMAT(LL) = BLW/HLW
         END DO
C
C --> spherical/non-spherical
C
         IF ( INS.EQ.0 ) THEN
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            DO LM1 = 1,LMMAXD
               TMATLL(LM1,LM1) = TMAT(LOFLM(LM1))
            END DO
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         ELSE
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            STOP ' not implemented'
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         END IF
         DEALLOCATE (BESSJW,BESSYW,HANKWS,DLOGDP,STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Deallocate'
         DEALLOCATE (TMAT,MASS,HAMF,FZ,PZ,STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Deallocate'
         DEALLOCATE (RS,S,STAT=LM1)
         IF ( LM1.NE.0 ) STOP '    Deallocate'
C ----------------------------------------------------------------------
      ELSE                      ! KREL
         CALL DRVRELTMAT(ERYD,TMATLL,VTREL,BTREL,RR,
     &                   DRDI,R2DRDI,ZREL,IRCUT(IPAN),SOLVER,SOCTL,
     &                   CTL,LMMAXD,LMAXD,IRMD)
      END IF
C ================================================================= KREL
      END
