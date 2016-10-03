      SUBROUTINE WRLDOS(DEN,EZ,WEZ,LMAXD1,IEMXD,NPOTD,ITITLE,EFERMI,E1,
     +                  E2,ALATC,TK,NACLS1,NSPINPOT,NATYP,CONC,IELAST,
     +                  INTERVX,INTERVY,INTERVZ,DOSTOT)
      use mod_version_info
      implicit none
C     .. Parameters ..
      DOUBLE PRECISION KB
      PARAMETER (KB=0.6333659D-5)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALATC,E1,E2,EFERMI,TK
      INTEGER IELAST,IEMXD,INTERVX,INTERVY,INTERVZ,LMAXD1,NACLS1,
     +        NATYP,NPOTD
C===== uses spin-up and down also in the REL mode (KREL=1)
      INTEGER NSPINPOT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DEN(0:LMAXD1,IEMXD,NPOTD),EZ(IEMXD),WEZ(IEMXD)
      DOUBLE PRECISION DOSTOT(0:LMAXD1,2),CONC(*) ! CONC(NATYPD)
      INTEGER ITITLE(20,NPOTD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DOSCMPLX
      DOUBLE PRECISION DOS,DOSSGN,EFCTOR,PI
      INTEGER I1,IA,IE,IPOT,ISPIN,L
      CHARACTER*8 DOSFL0
      CHARACTER*11 DOSFL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
      PI = 4.0D0*ATAN(1.0D0)
      DOSFL0 = 'dos.atom'
      EFCTOR = 1.0D0
      IF (TEST('EV      ')) EFCTOR = 13.6058D0
      DO ISPIN = 1,NSPINPOT
        DO L = 0,LMAXD1
          DOSTOT(L,ISPIN) = 0.0D0
        END DO
      END DO
      DO I1 = 1,NATYP
        IF (I1<10) WRITE (dosfl,FMT='(A8,I1)') DOSFL0,I1
        IF (I1>=10 .AND. I1<100) WRITE (dosfl,FMT='(A8,I2)') DOSFL0,I1
        IF (I1>=100) WRITE (dosfl,FMT='(A8,I3)') DOSFL0,I1
        OPEN (48,FILE=trim(DOSFL),FORM='formatted')
        call version_print_header(48)
        DO ISPIN = 1,NSPINPOT
            IPOT = NSPINPOT * (I1-1) + ISPIN
            DOSSGN = 1.0D0
            IF (ISPIN.NE.NSPINPOT) DOSSGN = -1.0D0

            WRITE (48,FMT=9010) (ITITLE(IA,IPOT),IA=1,19)
            WRITE (48,FMT=9020) I1
            WRITE (48,FMT=9030) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
            WRITE (48,FMT=9040) EFERMI
            WRITE (48,FMT=9050) TK,PI*KB*TK,ALATC,INTERVX,INTERVY,
     +           INTERVZ,NACLS1
            DO IE = 1,IELAST
                DOS = 0.0D0
                DO L = 0,LMAXD1
                    DOS = DOS - 2.0D0 * 
     &               DIMAG(DEN(L,IE,IPOT))/PI/DBLE(NSPINPOT)
                    DOSTOT(L,ISPIN) = DOSTOT(L,ISPIN) +
     +                   DIMAG(WEZ(IE)*DEN(L,IE,IPOT))
                END DO
                WRITE (48,FMT=9060) DBLE(EZ(IE))*EFCTOR,
     &                DOS*DOSSGN/EFCTOR,
     &               (-2.0D0*DIMAG(DEN(L,IE,IPOT))*DOSSGN/EFCTOR/PI
     &               /DBLE(NSPINPOT),L=0,LMAXD1)
            END DO
            WRITE (48,FMT=9070) (DOSTOT(L,ISPIN)/EFCTOR/DBLE(NSPINPOT),
     +                                                      L=0,LMAXD1)
            IF (ISPIN.NE.NSPINPOT) WRITE (48,FMT=9000)
        END DO
        CLOSE (48)
      END DO

c Write complex DOS in unit 49:
      OPEN (49,FILE='complex.dos',FORM='formatted')
      call version_print_header(49)
      WRITE (49,*) NATYP*NSPINPOT
      WRITE (49,*) IELAST
      WRITE (49,*) LMAXD1
      DO I1 = 1,NATYP

        DO ISPIN = 1,NSPINPOT
            IPOT = NSPINPOT * (I1-1) + ISPIN
            DOSSGN = 1.0D0
            IF (ISPIN.NE.NSPINPOT) DOSSGN = -1.0D0

            WRITE (49,FMT=9010) (ITITLE(IA,IPOT),IA=1,19)
            WRITE (49,FMT=9020) I1
            WRITE (49,FMT=9030) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
            WRITE (49,FMT=9040) EFERMI
            WRITE (49,FMT=9050) TK,PI*KB*TK,ALATC,INTERVX,INTERVY,
     +           INTERVZ,NACLS1
            DO IE = 1,IELAST
                DOSCMPLX = DCMPLX(0.0D0,0.D0)
                DO L = 0,LMAXD1
                    DOSCMPLX = DOSCMPLX - 2.0D0 * 
     &               DEN(L,IE,IPOT)/PI/DBLE(NSPINPOT)
                END DO
                WRITE (49,FMT=9065) EZ(IE)*EFCTOR,
     +               (-2.0D0*DEN(L,IE,IPOT)*DOSSGN/EFCTOR/PI
     &               /DBLE(NSPINPOT),L=0,LMAXD1),DOSCMPLX*DOSSGN/EFCTOR
            END DO
            IF (ISPIN.NE.NSPINPOT.OR.I1.NE.NATYP) WRITE (49,FMT=9000)
        END DO
      END DO
      CLOSE (49)

! Write total DOS summed over atoms and spins(complex)
      OPEN (49,FILE='total_cmplx.dos',FORM='formatted')
      call version_print_header(49)
      WRITE(49,FMT='(4A16)') '# Real(E)','  Im(E)',' Re(DEN)',' Im(DEN)'
      DO IE = 1,IELAST
         DOSCMPLX = DCMPLX(0.0D0,0.D0)
         DO I1 = 1,NATYP
            DO ISPIN = 1,NSPINPOT
               IPOT = NSPINPOT * (I1-1) + ISPIN
               DO L = 0,LMAXD1
                  DOSCMPLX = DOSCMPLX - CONC(I1) * 2.0D0 * 
     &               DEN(L,IE,IPOT)/PI/DBLE(NSPINPOT)
               ENDDO
            ENDDO
         ENDDO
         WRITE(49,FMT='(10E16.8)') EZ(IE),DOSCMPLX
      ENDDO
      CLOSE (49)

      
      RETURN
Cccc 9000 FORMAT ('&')
 9000 FORMAT (' ')
 9010 FORMAT ('#',19a4)
 9020 FORMAT ('# I1    :',I8)
 9030 FORMAT ('# ISPIN :',I8,'   IELAST :',I5,/,'# E1,E2 :',2f12.5,
     +       ' EFERMI :',f12.5,'   EFCTR',f10.6)
 9040 FORMAT ('# FERMI :',f12.5)
 9050 FORMAT ('# TK    =',f8.1,'   Kelvin =',3p,f8.3,' mRyd',0p,/,
     +       '# ALAT   :',f12.5,/,'# INTERV X,Y,Z  :',
     +       3I5,/,'# NACLS :',I8)
 9060 FORMAT (1p,8e15.7)
 9065 FORMAT (16('(',e12.4,',',e12.4,')'))
 9070 FORMAT ('# Integrated DOS ',1p,d10.3,7d11.3)
 9080 FORMAT ('&')
      END
