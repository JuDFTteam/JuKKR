C----------------------------------------------------------------------
C>    Writes local density of states.
      SUBROUTINE write_LDOS(DEN,EZ,LMAXD1,IEMXD,ITITLE,
     &                      EFERMI,E1,
     +                      E2,ALATC,TK,NSPIN,I1)

      IMPLICIT NONE

C     .. Parameters ..
      DOUBLE PRECISION KB
      PARAMETER       (KB=0.6333659D-5)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALATC,E1,E2,EFERMI,TK
      INTEGER          IELAST,IEMXD,LMAXD1
      INTEGER          NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   DEN(0:LMAXD1,IEMXD,NSPIN),
     +                 EZ(IEMXD),WEZ(IEMXD)
      DOUBLE PRECISION DOSTOT(0:LMAXD1,2),
     +                 PDOSTOT(0:LMAXD1,2)
      INTEGER          ITITLE(20,NSPIN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DOS,DOSSGN,EFCTOR,PI
      INTEGER          I1,IA,IE,IPOT,ISPIN,L
      CHARACTER(len=16) :: FNAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
      PI = 4.d0*ATAN(1.d0)
      EFCTOR = 1.d0
C     IF (TEST('EV      ')) EFCTOR = 13.6058D0
C
      IELAST = IEMXD
C
C initialize DOSTOT
C
      DO ISPIN = 1,NSPIN
        DO L = 0,LMAXD1
          DOSTOT(L,ISPIN) = 0.d0
          PDOSTOT(L,ISPIN) = 0.d0
        endDO
      endDO

C=======================================================================
C write DOS to file DOS.I1.dat - begin
C=======================================================================
C

      WRITE(UNIT=FNAME, FMT="(a,i4.4,a)") 'DOS.',I1,'.dat'

      OPEN(48, FILE=FNAME, FORM='formatted', ACTION='WRITE')
C
        DO ISPIN = 1,NSPIN
            IPOT = NSPIN * (I1-1) + ISPIN
            DOSSGN = 1.d0
            IF (ISPIN /= NSPIN) DOSSGN = -1.d0
C
            WRITE (48,FMT=9010) (ITITLE(IA,ISPIN),IA=1,19)
            WRITE (48,FMT=9020) I1
            WRITE (48,FMT=9030) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
            WRITE (48,FMT=9040) EFERMI
            WRITE (48,FMT=9050) TK,PI*KB*TK,ALATC
C
            DO IE = 1,IELAST
                DOS = 0.d0
                DO L = 0,LMAXD1
                    DOS = DOS - 2.d0 *
     &               DIMAG(DEN(L,IE,ISPIN))/PI/DBLE(NSPIN)
                    DOSTOT(L,ISPIN) = DOSTOT(L,ISPIN) +
     +                   DIMAG(WEZ(IE)*DEN(L,IE,ISPIN))
                endDO
                WRITE (48,FMT=9060) DBLE(EZ(IE))*EFCTOR,
     +               (-2.d0*DIMAG(DEN(L,IE,ISPIN))*DOSSGN/EFCTOR/PI
     &               /DBLE(NSPIN),L=0,LMAXD1),DOS*DOSSGN/EFCTOR
            endDO
C
            WRITE (48,FMT=9070) (DOSTOT(L,ISPIN)/EFCTOR/DBLE(NSPIN),
     +                                                      L=0,LMAXD1)
            IF (ISPIN /= NSPIN) WRITE (48,FMT=9000)
        endDO
        CLOSE (48)
C
C=======================================================================
C write DOS to file DOS.I1.dat - end
C=======================================================================
      RETURN
 9000 FORMAT ('&')
 9010 FORMAT ('#',19a4)
 9020 FORMAT ('# I1    :',I8)
 9030 FORMAT ('# ISPIN :',I8,'   IELAST :',I5,/,'# E1,E2 :',2f12.5,
     +       ' EFERMI :',f12.5,'   EFCTR',f10.6)
 9040 FORMAT ('# FERMI :',f12.5)
 9050 FORMAT ('# TK    =',f8.1,'   Kelvin =',3p,f8.3,' mRyd',0p,/,
     +       '# ALAT   :',f12.5)
 9060 FORMAT (1p,8e15.7)
 9065 FORMAT (1p,16d15.7)
 9070 FORMAT ('# Integrated DOS ',1p,d10.3,7d11.3)
 9080 FORMAT ('&')
      END

C> Writes complex.dos file (complex density of states).
      SUBROUTINE WRLDOS(DEN,EZ,WEZ,LMAXD1,IEMXD,NPOTD,ITITLE,EFERMI,E1,
     +                  E2,ALATC,TK,NSPIN,NAEZ,IELAST,I1,DOSTOT)

      IMPLICIT NONE

C     .. Parameters ..
      DOUBLE PRECISION KB
      PARAMETER       (KB=0.6333659D-5)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALATC,E1,E2,EFERMI,TK
      INTEGER          IELAST,IEMXD,LMAXD1,
     +                 NAEZ,NPOTD ! todo: remove NPOTD
      INTEGER          NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   DEN(0:LMAXD1,IEMXD,NSPIN),
     +                 EZ(IEMXD),WEZ(IEMXD) ! todo: remove WEZ
      DOUBLE PRECISION DOSTOT(0:LMAXD1,2),
     +                 PDOSTOT(0:LMAXD1,2)
      INTEGER          ITITLE(20,NSPIN)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   DOSCMPLX
      DOUBLE PRECISION DOSSGN,EFCTOR,PI
      INTEGER          I1,IA,IE,IPOT,ISPIN,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
      PI = 4.d0*ATAN(1.d0)
      EFCTOR = 1.d0

C initialize DOSTOT
C
      DO ISPIN = 1,NSPIN
        DO L = 0,LMAXD1
          DOSTOT(L,ISPIN) = 0.d0
          PDOSTOT(L,ISPIN) = 0.d0
        endDO
      endDO
C
C
C open file complex.dos - kept for correspondence to complexdos3.f
C
      IF(I1 == 1) THEN
      OPEN (49,FILE='complex.dos', FORM='formatted', ACTION='write')
      WRITE (49,*) NAEZ*NSPIN
      WRITE (49,*) IELAST
      WRITE (49,*) LMAXD1
      endIF

C=======================================================================
C Write complex DOS to file complex.dos - begin
C=======================================================================
C
        DO ISPIN = 1,NSPIN
            IPOT = NSPIN * (I1-1) + ISPIN
            DOSSGN = 1.d0
            IF (ISPIN /= NSPIN) DOSSGN = -1.d0
C
            WRITE (49,FMT=9010) (ITITLE(IA,ISPIN),IA=1,19)
            WRITE (49,FMT=9020) I1
            WRITE (49,FMT=9030) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
            WRITE (49,FMT=9040) EFERMI
            WRITE (49,FMT=9050) TK,PI*KB*TK,ALATC
C
            DO IE = 1,IELAST
                DOSCMPLX = DCMPLX(0.d0,0.D0)
                DO L = 0,LMAXD1
                    DOSCMPLX = DOSCMPLX - 2.d0 *
     &               DEN(L,IE,ISPIN)/PI/DBLE(NSPIN)
                endDO
                WRITE (49,FMT=9065) EZ(IE)*EFCTOR,
     +               (-2.d0*DEN(L,IE,ISPIN)*DOSSGN/EFCTOR/PI
     &               /DBLE(NSPIN),L=0,LMAXD1),DOSCMPLX*DOSSGN/EFCTOR
            endDO
C
            IF (ISPIN /= NSPIN.OR.I1 /= NAEZ) WRITE (49,FMT=9000)
        endDO
      IF(I1 == NAEZ) CLOSE (49)
C
C=======================================================================
C Write complex DOS to file complex.dos - end
C=======================================================================
      RETURN
C
 9000 FORMAT ('&')
 9010 FORMAT ('#',19a4)
 9020 FORMAT ('# I1    :',I8)
 9030 FORMAT ('# ISPIN :',I8,'   IELAST :',I5,/,'# E1,E2 :',2f12.5,
     +       ' EFERMI :',f12.5,'   EFCTR',f10.6)
 9040 FORMAT ('# FERMI :',f12.5)
 9050 FORMAT ('# TK    =',f8.1,'   Kelvin =',3p,f8.3,' mRyd',0p,/,
     +       '# ALAT   :',f12.5)
 9060 FORMAT (1p,8e15.7)
 9065 FORMAT (1p,16d15.7)
 9070 FORMAT ('# Integrated DOS ',1p,d10.3,7d11.3)
 9080 FORMAT ('&')
      END
