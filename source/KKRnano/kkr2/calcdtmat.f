      SUBROUTINE CALCDTMAT
     +                  (LDAU,NLDAU,ICST,
     +                   NSRA,EZ,IE,NPNT1,NPNT2,NPNT3,PI,TK,
     +                   DRDI,R,VINS,
     +                   VISP,ZAT,IPAN,
     +                   IRCUT,CLEB,LOFLM,ICLEB,IEND,
     <                   DTDE,TR_ALPH,LMAX,ISPIN,
     >                   LLDAU,WMLDAU,
C                        new input parameters after inc.p replace
     &                   nspind, ncleb, ipand, irmd, irnsd)
C
      IMPLICIT NONE


      INTEGER nspind
      INTEGER ncleb
      INTEGER ipand
      INTEGER irmd
      INTEGER irnsd

C      PARAMETER          (LMMAXD= (LMAXD+1)**2)
C      PARAMETER          (LMAXD1 = LMAXD + 1)
C      PARAMETER          (MMAXD=2*LMAXD+1)
C      PARAMETER          (LM2D= (2*LMAXD+1)**2)
C      PARAMETER          (LMPOTD= (LPOTD+1)**2)
C                                = (2*LMAX+1)**2)
C      PARAMETER          (IRMIND=IRMD-IRNSD)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION   TK,KB,PI,REALX,IMAGX
      INTEGER            LMAX,ISPIN,
     +                   NPNT1,NPNT2,NPNT3,
     +                   LM1,LM2,
     +                   IEND
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX     EZ,EZ1,EZ2
      DOUBLE PRECISION   CLEB(NCLEB,2)

C     DOUBLE PRECISION   VINS(IRMIND:IRMD,LMPOTD),
C    +                   VISP(IRMD),
C    +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      DOUBLE PRECISION   VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2),
     +                   VISP(IRMD),
     +                   WMLDAU(2*LMAX+1, 2*LMAX+1, NSPIND, LMAX + 1)
C     ..
C     DOUBLE COMPLEX     DTDE(LMMAXD,LMMAXD)
      DOUBLE COMPLEX     DTDE((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX     TR_ALPH,TR_ALPH1,TR_ALPH2,DZ

C     DOUBLE COMPLEX     TMATN1(LMMAXD,LMMAXD),TMATN2(LMMAXD,LMMAXD)
      DOUBLE COMPLEX     TMATN1((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX     TMATN2((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE PRECISION   DRDI(IRMD)
      DOUBLE PRECISION   R(IRMD)
      DOUBLE PRECISION   ZAT
C ----------------------------------------------------------------------
      INTEGER            IPAN,NLDAU

C     INTEGER            IRCUT(0:IPAND),LLDAU(LMAXD1)
      INTEGER            IRCUT(0:IPAND),LLDAU(LMAX + 1)
C     INTEGER            ICLEB(NCLEB,3),LOFLM(LM2D)
      INTEGER            ICLEB(NCLEB,3),LOFLM((2*LMAX+1)**2)

      INTEGER            ICST,NSRA,IE,IERR
      LOGICAL            LDAU
C     ..

      INTEGER             LMMAXD

      LMMAXD= (LMAX+1)**2

      KB = 0.6333659D-5
      IF (IE.LE.NPNT1 .OR. IE.GT.(NPNT1+NPNT2+NPNT3)) THEN
           DZ = DCMPLX(0.01D0*PI*KB*TK,0.0D0)
      ELSE
           DZ = DCMPLX(0.0D0,0.01D0*PI*KB*TK)
      END IF
C.. perpendicular to the contour
           EZ1 = EZ + DZ
C.. parallel to the contour
           EZ2 = EZ - DZ
 
           CALL CALCTMAT(LDAU,NLDAU,ICST,
     +                   NSRA,EZ1,
     +                   DRDI,R,VINS,
     +                   VISP,ZAT,IPAN,
     +                   IRCUT,CLEB,LOFLM,ICLEB,IEND,
     <                   TMATN1,TR_ALPH1,LMAX,ISPIN,
     >                   LLDAU,WMLDAU,
     &                   nspind, ncleb, ipand, irmd, irnsd)
 
           CALL CALCTMAT(LDAU,NLDAU,ICST,
     +                   NSRA,EZ2,
     +                   DRDI,R,VINS,
     +                   VISP,ZAT,IPAN,
     +                   IRCUT,CLEB,LOFLM,ICLEB,IEND,
     <                   TMATN2,TR_ALPH2,LMAX,ISPIN,
     >                   LLDAU,WMLDAU,
     &                   nspind, ncleb, ipand, irmd, irnsd)
C====================================================================
C
C dT(E)
C -----
C  dE
      DO LM2 = 1,LMMAXD
        DO LM1 = 1,LMMAXD
      DTDE(LM1,LM2) = (TMATN1(LM1,LM2)-TMATN2(LM1,LM2))*0.5D0/DZ
        ENDDO
      ENDDO

      TR_ALPH = -(TR_ALPH1-TR_ALPH2)*0.5D0/DZ
      RETURN
      END
 
