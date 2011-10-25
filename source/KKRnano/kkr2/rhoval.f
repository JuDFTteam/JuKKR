      SUBROUTINE RHOVAL(LDORHOEF,ICST,IELAST,NSRA,
     &                  ISPIN,NSPIN,
     &                  EZ,WEZ,DRDI,R,IRMIN,
     &                  VINS,VISP,ZAT,IPAN,IRCUT,
     &                  THETAS,IFUNM,LMSP,RHO2NS,R2NEF,DEN,
     &                  ESPV,CLEB,LOFLM,ICLEB,IEND,JEND,
     >                  GMATN,
     >                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,
     <                  DMATLDAU,
C                       new parameters after inc.p removal
     &                  iemxd,
     &                  lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
C
      IMPLICIT NONE

      INTEGER iemxd
      INTEGER lmaxd
      INTEGER irmd
      INTEGER ncleb
      INTEGER irnsd
      INTEGER irid
      INTEGER ipand
      INTEGER nfund
C
C     .. Parameters ..
C     INTEGER             LMXSPD
C     PARAMETER          (LMXSPD= (2*LPOTD+1)**2) ! = (4*LMAXD+1)**2
C     INTEGER             LMMAXD
C     PARAMETER          (LMMAXD= (LMAXD+1)**2)
C     INTEGER             LMAXD1
C     PARAMETER          (LMAXD1= LMAXD+1)
C     INTEGER             MMAXD
C     PARAMETER          (MMAXD = 2*LMAXD+1 )
C     INTEGER             LMPOTD
C     PARAMETER          (LMPOTD= (LPOTD+1)**2) ! = (2*LMAXD+1)**2
C     INTEGER             IRMIND
C     PARAMETER          (IRMIND=IRMD-IRNSD)
C     INTEGER             LM2D
C     PARAMETER          (LM2D= (2*LMAXD+1)**2)

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
C     DOUBLE COMPLEX     DEN(0:LMAXD1,IEMXD),EZ(IEMXD),
C    +                   WEZ(IEMXD),
C    +                   PHILDAU(IRMD,LMAXD1),
C    +                   DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
C     DOUBLE PRECISION   CLEB(NCLEB,2),DRDI(IRMD),
C    +                   ESPV(0:LMAXD1,1),
C    +                   R(IRMD),RHO2NS(IRMD,LMPOTD,2),
C    +                   R2NEF(IRMD,LMPOTD,2),   ! at fermi energy
C    +                   THETAS(IRID,NFUND),VINS(IRMIND:IRMD,LMPOTD),
C    +                   VISP(IRMD),
C    +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
C     INTEGER            ICLEB(NCLEB,3),IFUNM(LMXSPD),IRCUT(0:IPAND),
C    +                   JEND(LMPOTD,0:LMAXD,0:LMAXD),
C    +                   LMSP(LMXSPD),LOFLM(LM2D),
C    +                   LLDAU(LMAXD1)

      DOUBLE COMPLEX     DEN(0:LMAXD+1,IEMXD)
      DOUBLE COMPLEX     EZ(IEMXD)
      DOUBLE COMPLEX     WEZ(IEMXD)
      DOUBLE COMPLEX     PHILDAU(IRMD,LMAXD+1)
      DOUBLE COMPLEX     DMATLDAU(2*LMAXD+1,2*LMAXD+1,NSPIN,LMAXD+1)

      DOUBLE PRECISION   CLEB(NCLEB,2)
      DOUBLE PRECISION   DRDI(IRMD)
      DOUBLE PRECISION   ESPV(0:LMAXD+1,1)
      DOUBLE PRECISION   R(IRMD)
      DOUBLE PRECISION   RHO2NS(IRMD,(2*LMAXD+1)**2,2)
      DOUBLE PRECISION   R2NEF(IRMD,(2*LMAXD+1)**2,2)
      DOUBLE PRECISION   THETAS(IRID,NFUND)
      DOUBLE PRECISION   VINS(IRMD-IRNSD:IRMD,(2*LMAXD+1)**2)
      DOUBLE PRECISION   VISP(IRMD)
      DOUBLE PRECISION   WMLDAU(2*LMAXD+1,2*LMAXD+1,NSPIN,LMAXD+1)

      INTEGER            ICLEB(NCLEB,3)
      INTEGER            IFUNM((4*LMAXD+1)**2)
      INTEGER            IRCUT(0:IPAND)
      INTEGER            JEND((2*LMAXD+1)**2,0:LMAXD,0:LMAXD)
      INTEGER            LMSP((4*LMAXD+1)**2)
      INTEGER            LOFLM((2*LMAXD+1)**2)
      INTEGER            LLDAU(LMAXD+1)

C     .. Local Scalars ..
      DOUBLE COMPLEX     DF,ERYD,EK
      INTEGER            IDIM,IE,IR,L,LM1,LM2,
     +                   LMLO,LMHI,MMAX,IM,ILDAU
C     ..
C     .. Local Arrays ..
C     DOUBLE COMPLEX     ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),
C    +                   DR(LMMAXD,LMMAXD),
C    +                   CR(LMMAXD,LMMAXD),
C    +                   EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
C    +                   GMATLL(LMMAXD,LMMAXD),
C    +                   GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND),
C    +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
C    +                   PZ(IRMD,0:LMAXD),
C    +                   QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
C    +                   QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),
C    +                   TMAT(0:LMAXD)
C     DOUBLE PRECISION   RS(IRMD,0:LMAXD),S(0:LMAXD),
C    +                   LDAUCUT(IRMD),
C    +                   WMLDAUAV(LMAXD1)
C     DOUBLE COMPLEX     DENDUM(0:LMAXD1)

      DOUBLE COMPLEX    ALPHA(0:LMAXD)
      DOUBLE COMPLEX    AR((LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX    DR((LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX    CR((LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX    EKL(0:LMAXD)
      DOUBLE COMPLEX    FZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX    GMATLL((LMAXD+1)**2,(LMAXD+1)**2)
      DOUBLE COMPLEX    GMATN((LMAXD+1)**2,(LMAXD+1)**2,IEMXD,NSPIN)
      DOUBLE COMPLEX    PNS((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
      DOUBLE COMPLEX    PZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX    QNS((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
      DOUBLE COMPLEX    QZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX    SZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX    TMAT(0:LMAXD)

      DOUBLE PRECISION   RS(IRMD,0:LMAXD)
      DOUBLE PRECISION   S(0:LMAXD)
      DOUBLE PRECISION   LDAUCUT(IRMD)
      DOUBLE PRECISION   WMLDAUAV(LMAXD+1)

      DOUBLE COMPLEX     DENDUM(0:LMAXD+1)

      LOGICAL            TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DSCAL,CRADWF,PNSQNS,RHONS,WFMESH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG,SQRT
C     ..
C
      INTEGER             LMMAXD
      INTEGER             LMAXD1
      INTEGER             NSPIND
      INTEGER             LMPOTD

      LMPOTD = (2*LMAXD+1)**2
      LMMAXD= (LMAXD+1)**2
      LMAXD1= LMAXD+1
      NSPIND = NSPIN

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

         DO L = 0,LMAXD1
            ESPV(L,1) = 0.0D0
         END DO

      DO IE = 1,IELAST

         DO LM2 = 1,LMMAXD
            DO LM1 = 1,LMMAXD
               GMATLL(LM1,LM2) = GMATN(LM1,LM2,IE,ISPIN)
            END DO
         END DO

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
     &                  WMLDAU,WMLDAUAV,LDAUCUT,
     &                  lmaxd, nspind, irmd, irnsd, ipand, ncleb)


            DO L = 0,LMAXD
               EKL(L) = EK*DBLE(2*L+1)
            END DO
C-----------------------------------------------------------------------
            CALL RHONS(DEN(0,IE),DF,DRDI,GMATLL,EK,
     +           RHO2NS(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP,
     +           NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,
     +           JEND,IEND,EKL,
     &           lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

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
     +              JEND,IEND,EKL,
     &              lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
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
