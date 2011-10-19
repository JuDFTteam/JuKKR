      SUBROUTINE RHOVAL(LDORHOEF,ICST,IELAST,NSRA,
     &                  ISPIN,NSPIN,
     &                  EZ,WEZ,DRDI,R,IRMIN,
     &                  VINS,VISP,ZAT,IPAN,IRCUT,
     &                  THETAS,IFUNM,LMSP,RHO2NS,R2NEF,DEN,
     &                  ESPV,CLEB,LOFLM,ICLEB,IEND,JEND,
     >                  GMATN,
     >                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,
     <                  DMATLDAU)
C
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER             LMXSPD
      PARAMETER          (LMXSPD= (2*LPOTD+1)**2)
      INTEGER             LMMAXD
      PARAMETER          (LMMAXD= (LMAXD+1)**2)
      INTEGER             LMAXD1
      PARAMETER          (LMAXD1= LMAXD+1)
      INTEGER             MMAXD
      PARAMETER          (MMAXD = 2*LMAXD+1 )
      INTEGER             LMPOTD
      PARAMETER          (LMPOTD= (LPOTD+1)**2)
      INTEGER             IRMIND
      PARAMETER          (IRMIND=IRMD-IRNSD)
      INTEGER             LM2D
      PARAMETER          (LM2D= (2*LMAXD+1)**2)
      DOUBLE PRECISION    CVLIGHT
      PARAMETER          (CVLIGHT=274.0720442D0)
      DOUBLE COMPLEX      CONE 
      PARAMETER          ( CONE=(1D0,0D0) )
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ZAT
      INTEGER            ICST,IELAST,IEND,IPAN,ISPIN,NSPIN,NSRA,
     +                   IRMIN,NLDAU
      LOGICAL            LDORHOEF,LDAU
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX     DEN(0:LMAXD1,IEMXD),EZ(IEMXD),
     +                   WEZ(IEMXD),
     +                   PHILDAU(IRMD,LMAXD1),
     +                   DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      DOUBLE PRECISION   CLEB(NCLEB,2),DRDI(IRMD),
     +                   ESPV(0:LMAXD1,1),
     +                   R(IRMD),RHO2NS(IRMD,LMPOTD,2),
     +                   R2NEF(IRMD,LMPOTD,2),   ! at fermi energy
     +                   THETAS(IRID,NFUND),VINS(IRMIND:IRMD,LMPOTD),
     +                   VISP(IRMD),
     +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      INTEGER            ICLEB(NCLEB,3),IFUNM(LMXSPD),IRCUT(0:IPAND),
     +                   JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +                   LMSP(LMXSPD),LOFLM(LM2D),
     +                   LLDAU(LMAXD1)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX     DF,ERYD,EK
      INTEGER            IDIM,IE,IR,L,LM1,LM2,
     +                   LMLO,LMHI,MMAX,IM,ILDAU
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX     ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),
     +                   DR(LMMAXD,LMMAXD),
     +                   CR(LMMAXD,LMMAXD),
     +                   EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
     +                   GMATLL(LMMAXD,LMMAXD),
     +                   GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND),
     +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +                   PZ(IRMD,0:LMAXD),
     +                   QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +                   QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),
     +                   TMAT(0:LMAXD)
      DOUBLE PRECISION   RS(IRMD,0:LMAXD),S(0:LMAXD),
     +                   LDAUCUT(IRMD),
     +                   WMLDAUAV(LMAXD1)
      DOUBLE COMPLEX     DENDUM(0:LMAXD1)
      LOGICAL            TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DSCAL,CRADWF,PNSQNS,RHONS,WFMESH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG,SQRT
C     ..
C
C-----------------------------------------------------------------------
CLDAU
C
      IF (LDAU) THEN
        DO ILDAU=1,NLDAU
          WMLDAUAV(ILDAU) = 0.D0
          LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
          LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)
          MMAX = LMHI - LMLO + 1
          DO IM = 1,MMAX
            WMLDAUAV(ILDAU)=WMLDAUAV(ILDAU)+WMLDAU(IM,IM,ISPIN,ILDAU)
          ENDDO
          WMLDAUAV(ILDAU) = WMLDAUAV(ILDAU)/DBLE(MMAX)
        ENDDO
C
C -> Note: Application if WLDAU makes the potential discontinuous.
C    A cutoff can be used if necessary to make the potential continuous
C    for example (array bounds should be adjusted):
C
        IF(TEST('CUTOFF  ')) THEN
          DO IR = 1,IRMD
            LDAUCUT(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
     &                   ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
            LDAUCUT(IR) = 1D0/LDAUCUT(IR)
          ENDDO
        ELSE
          DO IR = 1,IRMD
            LDAUCUT(IR) = 1.D0
          ENDDO
        ENDIF
      ENDIF
C
CLDAU
C-----------------------------------------------------------------------



         DO LM1 = 1,LMPOTD
            DO IR = 1,IRMD
               RHO2NS(IR,LM1,ISPIN) = 0.0D0
               R2NEF(IR,LM1,ISPIN) = 0.0D0
            END DO
         END DO
C
         DO L = 0,LMAXD1
            ESPV(L,1) = 0.0D0
         END DO
C
      DO IE = 1,IELAST
C
         DO LM2 = 1,LMMAXD
            DO LM1 = 1,LMMAXD
               GMATLL(LM1,LM2) = GMATN(LM1,LM2,IE,ISPIN)
            END DO
         END DO
C
         ERYD = EZ(IE)
         DF = WEZ(IE)/DBLE(NSPIN)
C
C=======================================================================
            CALL WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN),
     &                  IRMD,LMAXD)
            CALL CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S,
     +                  PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,R,ZAT,
     >                  LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT,
     &                  lmaxd, irmd, ipand)
C-----------------------------------------------------------------------
C non-spherical
C
            CALL PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ,
     &                  PNS,QNS,NSRA,VINS,IPAN,IRCUT,
     &                  CLEB,ICLEB,IEND,LOFLM,LMAXD,ISPIN,
     &                  LDAU,NLDAU,LLDAU,
     &                  WMLDAU,WMLDAUAV,LDAUCUT)
C

            DO L = 0,LMAXD
               EKL(L) = EK*DBLE(2*L+1)
            END DO
C-----------------------------------------------------------------------
            CALL RHONS(DEN(0,IE),DF,DRDI,GMATLL,EK,
     +           RHO2NS(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP,
     +           NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,
     +           JEND,IEND,EKL)

C-----------------------------------------------------------------------


            DO L = 0,LMAXD1
               ESPV(L,1) = ESPV(L,1) + DIMAG(ERYD*DEN(L,IE)*DF)
            END DO
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     get the charge at the Fermi energy (IELAST)
C     call with the energy weight CONE --> not overwrite DF
C          with the dummy DENDUM       --> not overwrite DEN
C
            IF ( (IE.EQ.IELAST) .AND. (LDORHOEF) ) THEN
               CALL RHONS(DENDUM,CONE,DRDI,GMATLL,EK,
     +              R2NEF(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP,
     +              NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,
     +              JEND,IEND,EKL)
            END IF

C
        IF (LDAU.AND.NLDAU.GE.1) THEN
          CALL LDAUDMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL,
     >                  IPAN,IRCUT,DRDI,EK,
     >                  IRMIN,LLDAU,PHILDAU,NLDAU,
     <                  DMATLDAU,ISPIN,
     &                  lmaxd, nspind, irmd, irnsd, ipand)
C
        ENDIF


      END DO

C
      IF (ISPIN.EQ.2) THEN
         IDIM = IRMD*LMPOTD
         CALL DSCAL(IDIM,2.D0,RHO2NS(1,1,1),1)
         CALL DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1),1,RHO2NS(1,1,2),1)
         CALL DAXPY(IDIM,1.0D0,RHO2NS(1,1,2),1,RHO2NS(1,1,1),1)
C
C --> do the same at the Fermi energy
C
         CALL DSCAL(IDIM,2.D0,R2NEF(1,1,1),1)
         CALL DAXPY(IDIM,-0.5D0,R2NEF(1,1,1),1,R2NEF(1,1,2),1)
         CALL DAXPY(IDIM,1.0D0,R2NEF(1,1,2),1,R2NEF(1,1,1),1)
      END IF
C
      END
