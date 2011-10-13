      SUBROUTINE CALCTMAT(LDAU,NLDAU,ICST,
     +                   NSRA,EZ,
     +                   DRDI,R,VINS,VISP,ZAT,IPAN,
     +                   IRCUT,CLEB,LOFLM,ICLEB,IEND,
     <                   TMATN,TR_ALPH,LMAX,ISPIN,
     >                   LLDAU,WMLDAU)
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER             LMMAXD
      PARAMETER          (LMMAXD= (LMAXD+1)**2)
      INTEGER             LMAXD1
      PARAMETER          (LMAXD1 = LMAXD + 1)
      INTEGER             MMAXD
      PARAMETER          (MMAXD=2*LMAXD+1)
      INTEGER             LMPOTD
      PARAMETER          (LMPOTD= (LPOTD+1)**2)
      INTEGER             IRMIND
      PARAMETER          (IRMIND=IRMD-IRNSD)
      INTEGER             LM2D
      PARAMETER          (LM2D= (2*LMAXD+1)**2)
      DOUBLE PRECISION    CVLIGHT
      PARAMETER          (CVLIGHT=274.0720442D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ZAT
      INTEGER            ICST,IEND,IPAN,NSRA,LMAX,NLDAU,ISPIN
      DOUBLE COMPLEX     TR_ALPH,EZ
      LOGICAL            LDAU,TEST
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX     TMATN(LMMAXD,LMMAXD)
      DOUBLE PRECISION   CLEB(NCLEB,2),DRDI(IRMD),R(IRMD),
     +                   VINS(IRMIND:IRMD,LMPOTD),
     +                   VISP(IRMD),
     +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      INTEGER            ICLEB(NCLEB,3),IRCUT(0:IPAND),
     +                   LOFLM(LM2D),
     +                   LLDAU(LMAXD1)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX     ERYD,EK,DET
      DOUBLE PRECISION   PI
      INTEGER            LM1,LM2,L,IR,
     +                   LMLO,LMHI,MMAX,IM,ILDAU
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX     ALPHA(0:LMAXD),
     +                   FZ(IRMD,0:LMAXD),
     +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +                   PZ(IRMD,0:LMAXD),
     +                   QZ(IRMD,0:LMAXD),
     +                   SZ(IRMD,0:LMAXD),TMAT(0:LMAXD)
      DOUBLE PRECISION   RS(IRMD,0:LMAXD),S(0:LMAXD),
     +                   LDAUCUT(IRMD),
     +                   WMLDAUAV(LMAXD1)
C     ..
C     .. External Subroutines ..
      EXTERNAL CRADWF,PNSTMAT,WFMESH
C     ..
C     INTEGER MYRANK,NROFNODES
C     INTEGER MAPBLOCK
C     EXTERNAL MAPBLOCK
C     COMMON /MPI/MYRANK,NROFNODES

      PI = 4.D0*ATAN(1.D0)
C
CLDAU
C
      IF (LDAU) THEN
C
        DO ILDAU=1,NLDAU
C
          WMLDAUAV(ILDAU) = 0.0D0
          LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
          LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)
          MMAX = LMHI - LMLO + 1
          DO IM = 1,MMAX
            WMLDAUAV(ILDAU)=WMLDAUAV(ILDAU)+WMLDAU(IM,IM,ISPIN,ILDAU)
          ENDDO
          WMLDAUAV(ILDAU) = WMLDAUAV(ILDAU)/DBLE(MMAX)
C
        ENDDO
C
C -> Note: Application if WLDAU makes the potential discontinuous.
C    A cutoff can be used if necessary to make the potential continuous
C    for example (array bounds should be adjusted):
C
        IF(TEST('CUTOFF  ')) THEN
          DO IR = 1,IRMD
            LDAUCUT(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
     &                    ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
            LDAUCUT(IR) = 1D0/LDAUCUT(IR)
          ENDDO
        ELSE
          DO IR = 1,IRMD
            LDAUCUT(IR) = 1.D0
          ENDDO
        ENDIF
C
      ENDIF
C
CLDAU
C




         DO LM2 = 1,LMMAXD
            DO LM1 = 1,LMMAXD
               TMATN(LM1,LM2) = (0.0D0,0.0D0)
            END DO
         END DO
         ERYD = EZ
C
            CALL WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN),
     &                  IRMD,LMAXD)
            CALL CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S,
     &                  PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,R,ZAT,
     >                  LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT)
C-----------------------------------------------------------------------
            CALL PNSTMAT(DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,
     &                   TMATN,
     &                   VINS,IPAN,IRCUT,NSRA,CLEB,ICLEB,IEND,LOFLM,
     &                   TMAT,DET,LMAXD,ISPIN,
     >                   LDAU,NLDAU,LLDAU,
     >                   WMLDAU,WMLDAUAV,LDAUCUT)
C
            TR_ALPH = LOG(DET)
            DO L=0,LMAX 
            TR_ALPH = TR_ALPH + (L+L+1)*LOG(ALPHA(L))
            ENDDO


      RETURN
      END
