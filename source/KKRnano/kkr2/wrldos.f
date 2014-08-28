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
      INTEGER          I1,IA,IE,IPOT,ISPIN,L,D1,D10,D100,D1000,OFF(3)
      CHARACTER*12     FNAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
      PI = 4.0D0*ATAN(1.0D0)
      EFCTOR = 1.0D0
C     IF (TEST('EV      ')) EFCTOR = 13.6058D0
C
      IELAST = IEMXD
C
C initialize DOSTOT
C
      DO ISPIN = 1,NSPIN
        DO L = 0,LMAXD1
          DOSTOT(L,ISPIN) = 0.0D0
          PDOSTOT(L,ISPIN) = 0.0D0
        END DO
      END DO

C=======================================================================
C write DOS to file DOS.I1.dat - begin
C=======================================================================
C
      D1 = mod(I1,10)
      D10 = int( (mod(I1,100) + 0.5)/10 )
      D100 = int( (mod(I1,1000) + 0.5)/100 )
      D1000 = int( (mod(I1,10000) + 0.5)/1000 )

      OFF(1) = iachar('1')-1
      OFF(2) = iachar('1')-1
      OFF(3) = iachar('1')-1

      IF ( D10.GE.10 ) OFF(1) = iachar('7')
      IF ( D100.GE.100 ) OFF(2) = iachar('7')
      IF ( D1000.GE.1000 ) OFF(3) = iachar('7')

      FNAME='DOS.'
     + //achar(D1000+OFF(3))
     + //achar(D100+OFF(2))
     + //achar(D10+OFF(1))
     + //achar(D1+iachar('1')-1)
     + //'.dat'

      OPEN(48,FILE=FNAME,FORM='formatted')
C
        DO ISPIN = 1,NSPIN
            IPOT = NSPIN * (I1-1) + ISPIN
            DOSSGN = 1.0D0
            IF (ISPIN.NE.NSPIN) DOSSGN = -1.0D0
C
            WRITE (48,FMT=9010) (ITITLE(IA,ISPIN),IA=1,19)
            WRITE (48,FMT=9020) I1
            WRITE (48,FMT=9030) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
            WRITE (48,FMT=9040) EFERMI
            WRITE (48,FMT=9050) TK,PI*KB*TK,ALATC
C
            DO IE = 1,IELAST
                DOS = 0.0D0
                DO L = 0,LMAXD1
                    DOS = DOS - 2.0D0 *
     &               DIMAG(DEN(L,IE,ISPIN))/PI/DBLE(NSPIN)
                    DOSTOT(L,ISPIN) = DOSTOT(L,ISPIN) +
     +                   DIMAG(WEZ(IE)*DEN(L,IE,ISPIN))
                END DO
                WRITE (48,FMT=9060) DBLE(EZ(IE))*EFCTOR,
     +               (-2.0D0*DIMAG(DEN(L,IE,ISPIN))*DOSSGN/EFCTOR/PI
     &               /DBLE(NSPIN),L=0,LMAXD1),DOS*DOSSGN/EFCTOR
            END DO
C
            WRITE (48,FMT=9070) (DOSTOT(L,ISPIN)/EFCTOR/DBLE(NSPIN),
     +                                                      L=0,LMAXD1)
            IF (ISPIN.NE.NSPIN) WRITE (48,FMT=9000)
        END DO
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
     +                 NAEZ,NPOTD
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
      DOUBLE COMPLEX   DOSCMPLX
      DOUBLE PRECISION DOS,DOSSGN,EFCTOR,PI
      INTEGER          I1,IA,IE,IPOT,ISPIN,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
      PI = 4.0D0*ATAN(1.0D0)
      EFCTOR = 1.0D0

C initialize DOSTOT
C
      DO ISPIN = 1,NSPIN
        DO L = 0,LMAXD1
          DOSTOT(L,ISPIN) = 0.0D0
          PDOSTOT(L,ISPIN) = 0.0D0
        END DO
      END DO
C
C
C open file complex.dos - kept for correspondence to complexdos3.f
C
      IF(I1.EQ.1) THEN
      OPEN (49,FILE='complex.dos',FORM='formatted')
      WRITE (49,*) NAEZ*NSPIN
      WRITE (49,*) IELAST
      WRITE (49,*) LMAXD1
      END IF

C=======================================================================
C Write complex DOS to file complex.dos - begin
C=======================================================================
C
        DO ISPIN = 1,NSPIN
            IPOT = NSPIN * (I1-1) + ISPIN
            DOSSGN = 1.0D0
            IF (ISPIN.NE.NSPIN) DOSSGN = -1.0D0
C
            WRITE (49,FMT=9010) (ITITLE(IA,ISPIN),IA=1,19)
            WRITE (49,FMT=9020) I1
            WRITE (49,FMT=9030) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
            WRITE (49,FMT=9040) EFERMI
            WRITE (49,FMT=9050) TK,PI*KB*TK,ALATC
C
            DO IE = 1,IELAST
                DOSCMPLX = DCMPLX(0.0D0,0.D0)
                DO L = 0,LMAXD1
                    DOSCMPLX = DOSCMPLX - 2.0D0 *
     &               DEN(L,IE,ISPIN)/PI/DBLE(NSPIN)
                END DO
                WRITE (49,FMT=9065) EZ(IE)*EFCTOR,
     +               (-2.0D0*DEN(L,IE,ISPIN)*DOSSGN/EFCTOR/PI
     &               /DBLE(NSPIN),L=0,LMAXD1),DOSCMPLX*DOSSGN/EFCTOR
            END DO
C
            IF (ISPIN.NE.NSPIN.OR.I1.NE.NAEZ) WRITE (49,FMT=9000)
        END DO
      IF(I1.EQ.NAEZ) CLOSE (49)
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
