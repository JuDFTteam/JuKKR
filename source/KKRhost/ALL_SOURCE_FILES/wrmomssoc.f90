SUBROUTINE wrmomssoc(krel,natyp,nspin,texts,textl,textns,  &
    charge,muorb,lmaxd,lmaxd1,korbit,ormoment)

IMPLICIT NONE

! Dummy arguments
INTEGER KREL,LMAXD,LMAXD1,NATYP,NSPIN,KORBIT
CHARACTER*5 TEXTNS
REAL*8 CHARGE(0:LMAXD1,NATYP,2)
CHARACTER*4 TEXTL(0:6)
CHARACTER*7 TEXTS(3)

! Local variables
REAL*8 CHTOT(NATYP)
DOUBLE PRECISION MUORB(0:LMAXD1+1,3,NATYP)
DOUBLE PRECISION MUSPIN(NATYP,0:LMAXD1+1), &
                 MUTOT(NATYP),SUMCH(NATYP,2), &
                 ORMOMENT(NSPIN,4,NATYP)
CHARACTER*80 FMT1,FMT2,FMT31,FMT32
INTEGER IS,ISPIN,IT,L,LF1,LF2


WRITE (1337,*)

IF ((krel == 1).OR.(nspin == 2)) THEN
  WRITE (1337,'(78(1H#))')
  WRITE (1337,99001)
  WRITE (1337,'(78(1H#))')
ELSE
  WRITE (1337,'(44(1H#))')
  WRITE (1337,99002)
  WRITE (1337,'(44(1H#))')
END IF

WRITE (1337,*)
WRITE (1337,99003)
DO it = 1,natyp
  muspin(it,lmaxd1+1) = 0D0
  sumch(it,1) = 0D0
  sumch(it,2) = 0D0
  DO l = 0,lmaxd1
    DO ispin = 1,nspin
      sumch(it,ispin) = sumch(it,ispin) + charge(l,it,ispin)
    END DO
    muspin(it,l) = charge(l,it,2) - charge(l,it,1)
    muspin(it,lmaxd1+1) = muspin(it,lmaxd1+1) + muspin(it,l)
  END DO
  chtot(it) = sumch(it,1) + sumch(it,2)
END DO

IF ( krel == 1 ) THEN
  DO it = 1,natyp
    mutot(it) = muspin(it,lmaxd1+1) + muorb(lmaxd1+1,3,it)
  END DO
END IF

is = 0
IF ( nspin == 1 ) is = is + 2
DO ispin = 1,nspin
  is = is + 1
  WRITE (1337,99004) texts(is)
END DO

IF (krel == 1) THEN
  WRITE (1337,99005)
  WRITE (1337,99006)
ELSE
  IF (nspin == 2) WRITE(1337,99005)
  WRITE(1337,*)
END IF

WRITE (1337,'(3X,26(1H=),$)')
IF ( krel == 1 ) THEN
  WRITE (1337,'(46(1H=))')
ELSE
  IF (nspin == 2) WRITE (1337,'(23(1H=),$)')
  WRITE(1337,*)
END IF

fmt1 = '(4X,I3,2X,A4,2(F12.8),2X,F8.4'
fmt2 = '(9X,A4,2(F12.8),2X,F8.4'
fmt31 = '(4X,I3,2X,A4,F12.8)'
fmt32 = '(9X,A4,F12.8)'
lf1 = 30
lf2 = 24

IF ( krel == 1.OR.korbit == 1 ) THEN
  fmt1 = fmt1(1:lf1)//',2X,3F8.4)'
  fmt2 = fmt2(1:lf2)//',2X,3F8.4)'
ELSE
  IF (nspin == 2) THEN
    fmt1 = fmt1(1:lf1)//')'
    fmt2 = fmt2(1:lf2)//')'
  ELSE
    fmt1 = fmt31
    fmt2 = fmt32
  END IF
END IF

DO it = 1,natyp
  IF ( krel == 1 ) THEN
    WRITE (1337,FMT=fmt1) it,textl(0), (charge(0,it,ispin),ispin=1,nspin),  &
        muspin(it,0),muorb(0,3,it), (muorb(0,ispin,it),ispin=1,nspin)
  ELSE
    IF (nspin == 2) THEN
      WRITE (1337,FMT=fmt1) it,textl(0),  &
          (charge(0,it,ispin),ispin=1,nspin), muspin(it,0)
    ELSE
      WRITE (1337,FMT=fmt1) it,textl(0),charge(0,it,1)
    END IF
  END IF
  
  DO l = 1,lmaxd
    IF ( krel == 1 ) THEN
      WRITE (1337,FMT=fmt2) textl(l), (charge(l,it,ispin),ispin=1,nspin),  &
          muspin(it,l),muorb(l,3,it), (muorb(l,ispin,it),ispin=1,nspin)
    ELSE
      IF (nspin == 2) THEN
        WRITE (1337,FMT=fmt2) textl(l), (charge(l,it,ispin),ispin=1,nspin),  &
            muspin(it,l)
      ELSE
        WRITE (1337,FMT=fmt2) textl(l),charge(l,it,1)
      END IF
    END IF
    
  END DO
  
  IF ( krel == 1 ) THEN
    WRITE (1337,FMT=fmt2) textns, (charge(lmaxd1,it,ispin),ispin=1,nspin),  &
        muspin(it,lmaxd1),muorb(lmaxd1,3,it),  &
        (muorb(lmaxd1,ispin,it),ispin=1,nspin)
  ELSE
    IF (nspin == 2) THEN
      WRITE (1337,FMT=fmt2) textns,  &
          (charge(lmaxd1,it,ispin),ispin=1,nspin), muspin(it,lmaxd1)
    ELSE
      WRITE (1337,FMT=fmt2) textns,charge(lmaxd1,it,1)
    END IF
  END IF
  
  WRITE (1337,'(10x,19(1H-),$)')
  IF ( krel == 1 ) THEN
    WRITE (1337,'(44(1H-))')
    WRITE (1337,FMT=fmt2) ' TOT',(sumch(it,ispin),ispin=1,nspin)  &
        ,muspin(it,lmaxd1+1),muorb(lmaxd1+1,3,it),  &
        (muorb(lmaxd1+1,ispin,it),ispin=1,nspin)
    WRITE (1337,'(25X,F12.8,12X,F8.4)') chtot(it),mutot(it)
  ELSE IF ( korbit == 1 ) THEN
    WRITE (1337,'(44(1H-))')
    WRITE (1337,FMT=fmt2) ' TOT',(sumch(it,ispin),ispin=1,nspin)  &
        ,muspin(it,lmaxd1+1), ormoment(1,4,it)+ormoment(2,4,it),  &
        (ormoment(ispin,4,it),ispin=1,nspin)
  ELSE
    IF (nspin == 2) THEN
      WRITE (1337,'(17(1H-))')
      WRITE (1337,FMT=fmt2) ' TOT',  &
          (sumch(it,ispin),ispin=1,nspin),muspin(it,lmaxd1+1)
      WRITE (1337,'(25X,F12.8)') chtot(it)
    ELSE
      WRITE (1337,*)
      WRITE (1337,FMT=fmt2) ' TOT',sumch(it,1)
    END IF
  END IF
  
  IF ( it /= natyp ) THEN
    WRITE (1337,'(3X,26(1H=),$)')
    IF ( krel == 1 ) THEN
      WRITE (1337,'(40(1H=))')
    ELSE
      IF (nspin == 2) WRITE(1337,'(17(1H=),$)')
      WRITE (1337,*)
    END IF
  END IF
END DO

WRITE (1337,*)
IF ((krel == 1).OR.(nspin == 2)) THEN
  WRITE (1337,'(78(1H#))')
ELSE
  WRITE (1337,'(44(1H#))')
END IF
WRITE (1337,*)

RETURN

99001 FORMAT (15X,'l-decomposed valence charges and magnetic moments')
99002 FORMAT (8X,'l-decomposed valence charges')
99003 FORMAT (3X,'ATOM      ',$)
99004 FORMAT (2X,'Ne ',a7,$)
99005 FORMAT ('    m_spin',$)
99006 FORMAT ('    m_orb   spin dn  spin up')
END SUBROUTINE wrmomssoc
