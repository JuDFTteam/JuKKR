C*==wrmoms.f    processed by SPAG 6.05Rc at 15:39 on  1 Mar 2002
      SUBROUTINE WRMOMSSOC(KREL,NATYP,NSPIN,TEXTS,TEXTL,TEXTNS,
     &                  CHARGE,MUORB,LMAXD,LMAXD1,KORBIT,ORMOMENT)
C
      IMPLICIT NONE
C
C Dummy arguments
C
      INTEGER KREL,LMAXD,LMAXD1,NATYP,NSPIN,KORBIT
      CHARACTER*5 TEXTNS
      REAL*8 CHARGE(0:LMAXD1,NATYP,2)
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*7 TEXTS(3)
C
C Local variables
C
      REAL*8 CHTOT(NATYP)
      DOUBLE PRECISION MUORB(0:LMAXD1+1,3,NATYP)
      DOUBLE PRECISION MUSPIN(NATYP,0:LMAXD1+1),
     &                 MUTOT(NATYP),SUMCH(NATYP,2),
     &                 ORMOMENT(NSPIN,4,NATYP)
      CHARACTER*80 FMT1,FMT2,FMT31,FMT32
      INTEGER IS,ISPIN,IT,L,LF1,LF2
C
      WRITE (1337,*)
      IF ((KREL.EQ.1).OR.(NSPIN.EQ.2)) THEN
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
      DO IT = 1,NATYP
         MUSPIN(IT,LMAXD1+1) = 0D0
         SUMCH(IT,1) = 0D0
         SUMCH(IT,2) = 0D0
         DO L = 0,LMAXD1
            DO ISPIN = 1,NSPIN
               SUMCH(IT,ISPIN) = SUMCH(IT,ISPIN) + CHARGE(L,IT,ISPIN)
            END DO
            MUSPIN(IT,L) = CHARGE(L,IT,2) - CHARGE(L,IT,1)
            MUSPIN(IT,LMAXD1+1) = MUSPIN(IT,LMAXD1+1) + MUSPIN(IT,L)
         END DO
         CHTOT(IT) = SUMCH(IT,1) + SUMCH(IT,2)
      END DO
C
      IF ( KREL.EQ.1 ) THEN
         DO IT = 1,NATYP
            MUTOT(IT) = MUSPIN(IT,LMAXD1+1) + MUORB(LMAXD1+1,3,IT)
         END DO
      END IF
C
      IS = 0
      IF ( NSPIN.EQ.1 ) IS = IS + 2
      DO ISPIN = 1,NSPIN
         IS = IS + 1
         WRITE (1337,99004) TEXTS(IS)
      END DO
C
      IF (KREL.EQ.1) THEN
         WRITE (1337,99005)
         WRITE (1337,99006)
      ELSE
         IF (NSPIN.EQ.2) WRITE(1337,99005)
         WRITE(1337,*)
      END IF
C
      WRITE (1337,'(3X,26(1H=),$)')
      IF ( KREL.EQ.1 ) THEN
         WRITE (1337,'(46(1H=))')
      ELSE
         IF (NSPIN.EQ.2) WRITE (1337,'(23(1H=),$)')
         WRITE(1337,*)
      END IF
C
      FMT1 = '(4X,I3,2X,A4,2(F12.8),2X,F8.4'
      FMT2 = '(9X,A4,2(F12.8),2X,F8.4'
      FMT31 = '(4X,I3,2X,A4,F12.8)'
      FMT32 = '(9X,A4,F12.8)'
      LF1 = 30
      LF2 = 24
C
      IF ( KREL.EQ.1.OR.KORBIT.EQ.1 ) THEN
         FMT1 = FMT1(1:LF1)//',2X,3F8.4)'
         FMT2 = FMT2(1:LF2)//',2X,3F8.4)'
      ELSE
         IF (NSPIN.EQ.2) THEN
            FMT1 = FMT1(1:LF1)//')'
            FMT2 = FMT2(1:LF2)//')'
         ELSE
            FMT1 = FMT31
            FMT2 = FMT32
         END IF
      END IF
C
      DO IT = 1,NATYP
         IF ( KREL.EQ.1 ) THEN
            WRITE (1337,FMT=FMT1) IT,TEXTL(0),
     &                         (CHARGE(0,IT,ISPIN),ISPIN=1,NSPIN),
     &                         MUSPIN(IT,0),MUORB(0,3,IT),
     &                         (MUORB(0,ISPIN,IT),ISPIN=1,NSPIN)
         ELSE 
            IF (NSPIN.EQ.2) THEN
               WRITE (1337,FMT=FMT1) IT,TEXTL(0),
     &              (CHARGE(0,IT,ISPIN),ISPIN=1,NSPIN),
     &              MUSPIN(IT,0)
            ELSE
               WRITE (1337,FMT=FMT1) IT,TEXTL(0),CHARGE(0,IT,1)
            END IF
         END IF
C
         DO L = 1,LMAXD
            IF ( KREL.EQ.1 ) THEN
               WRITE (1337,FMT=FMT2) TEXTL(L),
     &                            (CHARGE(L,IT,ISPIN),ISPIN=1,NSPIN),
     &                            MUSPIN(IT,L),MUORB(L,3,IT),
     &                            (MUORB(L,ISPIN,IT),ISPIN=1,NSPIN)
            ELSE
               IF (NSPIN.EQ.2) THEN
                  WRITE (1337,FMT=FMT2) TEXTL(L),
     &                 (CHARGE(L,IT,ISPIN),ISPIN=1,NSPIN),
     &                 MUSPIN(IT,L)
               ELSE
                  WRITE (1337,FMT=FMT2) TEXTL(L),CHARGE(L,IT,1)
               END IF
            END IF
C
         END DO
C
         IF ( KREL.EQ.1 ) THEN
            WRITE (1337,FMT=FMT2) TEXTNS,
     &                         (CHARGE(LMAXD1,IT,ISPIN),ISPIN=1,NSPIN),
     &                         MUSPIN(IT,LMAXD1),MUORB(LMAXD1,3,IT),
     &                         (MUORB(LMAXD1,ISPIN,IT),ISPIN=1,NSPIN)
         ELSE
            IF (NSPIN.EQ.2) THEN
               WRITE (1337,FMT=FMT2) TEXTNS,
     &              (CHARGE(LMAXD1,IT,ISPIN),ISPIN=1,NSPIN),
     &              MUSPIN(IT,LMAXD1)
            ELSE
               WRITE (1337,FMT=FMT2) TEXTNS,CHARGE(LMAXD1,IT,1)
            END IF
         END IF
C
         WRITE (1337,'(10x,19(1H-),$)')
         IF ( KREL.EQ.1 ) THEN
            WRITE (1337,'(44(1H-))')
            WRITE (1337,FMT=FMT2) ' TOT',(SUMCH(IT,ISPIN),ISPIN=1,NSPIN)
     &                        ,MUSPIN(IT,LMAXD1+1),MUORB(LMAXD1+1,3,IT),
     &                         (MUORB(LMAXD1+1,ISPIN,IT),ISPIN=1,NSPIN)
            WRITE (1337,'(25X,F12.8,12X,F8.4)') CHTOT(IT),MUTOT(IT)
         ELSEIF ( KORBIT.EQ.1 ) THEN
            WRITE (1337,'(44(1H-))')
            WRITE (1337,FMT=FMT2) ' TOT',(SUMCH(IT,ISPIN),ISPIN=1,NSPIN)
     &                        ,MUSPIN(IT,LMAXD1+1),
     &                         ORMOMENT(1,4,IT)+ORMOMENT(2,4,IT),
     &                         (ORMOMENT(ISPIN,4,IT),ISPIN=1,NSPIN)
         ELSE
            IF (NSPIN.EQ.2) THEN
               WRITE (1337,'(17(1H-))')
               WRITE (1337,FMT=FMT2) ' TOT',
     &              (SUMCH(IT,ISPIN),ISPIN=1,NSPIN),MUSPIN(IT,LMAXD1+1)
               WRITE (1337,'(25X,F12.8)') CHTOT(IT)
            ELSE
               WRITE (1337,*)
               WRITE (1337,FMT=FMT2) ' TOT',SUMCH(IT,1)
            END IF
         END IF
C
         IF ( IT.NE.NATYP ) THEN
            WRITE (1337,'(3X,26(1H=),$)')
            IF ( KREL.EQ.1 ) THEN
               WRITE (1337,'(40(1H=))')
            ELSE
               IF (NSPIN.EQ.2) WRITE(1337,'(17(1H=),$)')
               WRITE (1337,*)
            END IF
         END IF
      END DO
C
      WRITE (1337,*)
      IF ((KREL.EQ.1).OR.(NSPIN.EQ.2)) THEN
         WRITE (1337,'(78(1H#))')
      ELSE
         WRITE (1337,'(44(1H#))')
      END IF
      WRITE (1337,*)
C
      RETURN
C
99001 FORMAT (15X,'l-decomposed valence charges and magnetic moments')
99002 FORMAT (8X,'l-decomposed valence charges')
99003 FORMAT (3X,'ATOM      ',$)
99004 FORMAT (2X,'Ne ',A7,$)
99005 FORMAT ('    m_spin',$)
99006 FORMAT ('    m_orb   spin dn  spin up')
      END
