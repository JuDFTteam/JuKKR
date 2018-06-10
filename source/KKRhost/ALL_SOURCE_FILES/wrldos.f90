SUBROUTINE wrldos(den,ez,wez,lmaxd1,iemxd,npotd,ititle,efermi,e1,  &
        e2,alatc,tk,nacls1,nspinpot,natyp,conc,ielast,  &
        intervx,intervy,intervz,dostot)
      use mod_version_info
      implicit none
!.. Parameters ..
      DOUBLE PRECISION KB
      PARAMETER (KB=0.6333659D-5)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ALATC,E1,E2,EFERMI,TK
      INTEGER IELAST,IEMXD,INTERVX,INTERVY,INTERVZ,LMAXD1,NACLS1, &
              NATYP,NPOTD
!===== uses spin-up and down also in the REL mode (KREL=1)
      INTEGER NSPINPOT
!..
!.. Array Arguments ..
      DOUBLE COMPLEX DEN(0:LMAXD1,IEMXD,NPOTD),EZ(IEMXD),WEZ(IEMXD)
      DOUBLE PRECISION DOSTOT(0:LMAXD1,2),CONC(*) ! CONC(NATYPD)
      INTEGER ITITLE(20,NPOTD)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX DOSCMPLX
      DOUBLE PRECISION DOS,DOSSGN,EFCTOR,PI
      INTEGER I1,IA,IE,IPOT,ISPIN,L
      CHARACTER*8 DOSFL0
      CHARACTER*11 DOSFL
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG
!..
!.. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
!     ..
pi = 4.0D0*ATAN(1.0D0)
dosfl0 = 'dos.atom'
efctor = 1.0D0
IF (test('EV      ')) efctor = 13.6058D0
DO ispin = 1,nspinpot
  DO l = 0,lmaxd1
    dostot(l,ispin) = 0.0D0
  END DO
END DO
DO i1 = 1,natyp
  IF (i1<10) WRITE (dosfl,FMT='(A8,I1)') dosfl0,i1
  IF (i1>=10 .AND. i1<100) WRITE (dosfl,FMT='(A8,I2)') dosfl0,i1
  IF (i1>=100) WRITE (dosfl,FMT='(A8,I3)') dosfl0,i1
  OPEN (48,FILE=trim(dosfl),FORM='formatted')
  CALL version_print_header(48)
  DO ispin = 1,nspinpot
    ipot = nspinpot * (i1-1) + ispin
    dossgn = 1.0D0
    IF (ispin /= nspinpot) dossgn = -1.0D0
    
    WRITE (48,FMT=9010) (ititle(ia,ipot),ia=1,19)
    WRITE (48,FMT=9020) i1
    WRITE (48,FMT=9030) ispin,ielast,e1,e2,efermi,efctor
    WRITE (48,FMT=9040) efermi
    WRITE (48,FMT=9050) tk,pi*kb*tk,alatc,intervx,intervy, intervz,nacls1
    DO ie = 1,ielast
      dos = 0.0D0
      DO l = 0,lmaxd1
        dos = dos - 2.0D0 * DIMAG(den(l,ie,ipot))/pi/DBLE(nspinpot)
        dostot(l,ispin) = dostot(l,ispin) + DIMAG(wez(ie)*den(l,ie,ipot))
      END DO
      WRITE (48,FMT=9060) DBLE(ez(ie))*efctor, dos*dossgn/efctor,  &
          (-2.0D0*DIMAG(den(l,ie,ipot))*dossgn/efctor/pi  &
          /DBLE(nspinpot),l=0,lmaxd1)
    END DO
    WRITE (48,FMT=9070) (dostot(l,ispin)/efctor/DBLE(nspinpot), l=0,lmaxd1)
    IF (ispin /= nspinpot) WRITE (48,FMT=9000)
  END DO
  CLOSE (48)
END DO

! Write complex DOS in unit 49:
OPEN (49,FILE='complex.dos',FORM='formatted')
CALL version_print_header(49)
WRITE (49,*) natyp*nspinpot
WRITE (49,*) ielast
WRITE (49,*) lmaxd1
DO i1 = 1,natyp
  
  DO ispin = 1,nspinpot
    ipot = nspinpot * (i1-1) + ispin
    dossgn = 1.0D0
    IF (ispin /= nspinpot) dossgn = -1.0D0
    
    WRITE (49,FMT=9010) (ititle(ia,ipot),ia=1,19)
    WRITE (49,FMT=9020) i1
    WRITE (49,FMT=9030) ispin,ielast,e1,e2,efermi,efctor
    WRITE (49,FMT=9040) efermi
    WRITE (49,FMT=9050) tk,pi*kb*tk,alatc,intervx,intervy, intervz,nacls1
    DO ie = 1,ielast
      doscmplx = DCMPLX(0.0D0,0.d0)
      DO l = 0,lmaxd1
        doscmplx = doscmplx - 2.0D0 * den(l,ie,ipot)/pi/DBLE(nspinpot)
      END DO
      WRITE (49,FMT=9065) ez(ie)*efctor,  &
          (-2.0D0*den(l,ie,ipot)*dossgn/efctor/pi  &
          /DBLE(nspinpot),l=0,lmaxd1),doscmplx*dossgn/efctor
    END DO
    IF (ispin /= nspinpot.OR.i1 /= natyp) WRITE (49,FMT=9000)
  END DO
END DO
CLOSE (49)

! Write total DOS summed over atoms and spins(complex)
OPEN (49,FILE='total_cmplx.dos',FORM='formatted')
CALL version_print_header(49)
WRITE(49,FMT='(4A16)') '# Real(E)','  Im(E)',' Re(DEN)',' Im(DEN)'
DO ie = 1,ielast
  doscmplx = DCMPLX(0.0D0,0.d0)
  DO i1 = 1,natyp
    DO ispin = 1,nspinpot
      ipot = nspinpot * (i1-1) + ispin
      DO l = 0,lmaxd1
        doscmplx = doscmplx - conc(i1) * 2.0D0 *  &
            den(l,ie,ipot)/pi/DBLE(nspinpot)
      END DO
    END DO
  END DO
  WRITE(49,FMT='(10E16.8)') ez(ie),doscmplx
END DO
CLOSE (49)


RETURN
!ccc 9000 FORMAT ('&')
9000 FORMAT (' ')
9010 FORMAT ('#',19A4)
9020 FORMAT ('# I1    :',i8)
9030 FORMAT ('# ISPIN :',i8,'   IELAST :',i5,/,'# E1,E2 :',2F12.5,  &
    ' EFERMI :',f12.5,'   EFCTR',f10.6)
9040 FORMAT ('# FERMI :',f12.5)
9050 FORMAT ('# TK    =',f8.1,'   Kelvin =',3P,f8.3,' mRyd',0P,/,  &
    '# ALAT   :',f12.5,/,'# INTERV X,Y,Z  :', 3I5,/,'# NACLS :',i8)
9060 FORMAT (1P,8E15.7)
9065 FORMAT (16('(',e12.4,',',e12.4,')'))
9070 FORMAT ('# Integrated DOS ',1P,d10.3,7D11.3)
9080 FORMAT ('&')
END SUBROUTINE wrldos
