C*==wrmoms.f    processed by SPAG 6.05Rc at 15:39 on  1 Mar 2002
      SUBROUTINE WRMOMS(NAEZ,NSPIN,CHARGE,I1,LMAXD,LMAXD1)
C
      IMPLICIT NONE

C
C Dummy arguments
C
      INTEGER LMAXD,LMAXD1,NAEZ,NSPIN
      CHARACTER*5 TEXTNS
      REAL*8 CHARGE(0:LMAXD1,2)
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*7 TEXTS(3)
C
C Local variables
C
      DOUBLE PRECISION MUSPIN(0:LMAXD1+1),
     &                 SUMCH(2)
      CHARACTER*80 FMT1,FMT2,FMT31,FMT32
      INTEGER IS,ISPIN,I1,L,LF1,LF2
C     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/'spin dn','spin up','       '/
      DATA TEXTNS/' ns ='/
C     ..
C
      IF(I1.EQ.1) THEN

      WRITE (6,*)
      IF (NSPIN.EQ.2) THEN
         WRITE (6,'(78(1H#))')
         WRITE (6,99001)
         WRITE (6,'(78(1H#))')
      ELSE
         WRITE (6,'(44(1H#))')
         WRITE (6,99002)
         WRITE (6,'(44(1H#))')
      END IF
      
      WRITE (6,*)
      WRITE (6,99003)

C
      IS = 0
      IF ( NSPIN.EQ.1 ) IS = IS + 2
      DO ISPIN = 1,NSPIN
         IS = IS + 1
         WRITE (6,99004) TEXTS(IS)
      END DO
C
      IF (NSPIN.EQ.2) WRITE(6,99005)
         WRITE(6,*)
C
      WRITE (6,'(9X,26(1H=),$)')
      IF (NSPIN.EQ.2) WRITE (6,'(23(1H=),$)')
         WRITE(6,*)

      END IF
C
      FMT1 = '(8x,i5,2x,a4,2(f14.8),2x,f8.4'
      FMT2 = '(15x,a4,2(f14.8),2x,f8.4'
      FMT31 = '(10x,i3,2x,a4,f14.8)'
      FMT32 = '(15x,a4,f14.8)'
      LF1 = 30
      LF2 = 24
C
         IF (NSPIN.EQ.2) THEN
            FMT1 = FMT1(1:LF1)//')'
            FMT2 = FMT2(1:LF2)//')'
         ELSE
            FMT1 = FMT31
            FMT2 = FMT32
         END IF
C

         MUSPIN(LMAXD1+1) = 0D0
         SUMCH(1) = 0D0
         SUMCH(2) = 0D0
         DO L = 0,LMAXD1
            DO ISPIN = 1,NSPIN
               SUMCH(ISPIN) = SUMCH(ISPIN) + CHARGE(L,ISPIN)
            END DO
            MUSPIN(L) = CHARGE(L,2) - CHARGE(L,1)
            MUSPIN(LMAXD1+1) = MUSPIN(LMAXD1+1) + MUSPIN(L)
         END DO


            IF (NSPIN.EQ.2) THEN
               WRITE (6,FMT=FMT1) I1,TEXTL(0),
     &              (CHARGE(0,ISPIN),ISPIN=1,NSPIN),
     &              MUSPIN(0)
            ELSE
               WRITE (6,FMT=FMT1) I1,TEXTL(0),CHARGE(0,1)
            END IF
C
         DO L = 1,LMAXD
               IF (NSPIN.EQ.2) THEN
                  WRITE (6,FMT=FMT2) TEXTL(L),
     &                 (CHARGE(L,ISPIN),ISPIN=1,NSPIN),
     &                 MUSPIN(L)
               ELSE
                  WRITE (6,FMT=FMT2) TEXTL(L),CHARGE(L,1)
               END IF
C
         END DO
C
            IF (NSPIN.EQ.2) THEN
               WRITE (6,FMT=FMT2) TEXTNS,
     &              (CHARGE(LMAXD1,ISPIN),ISPIN=1,NSPIN),
     &              MUSPIN(LMAXD1)
            ELSE
               WRITE (6,FMT=FMT2) TEXTNS,CHARGE(LMAXD1,1)
            END IF
         WRITE (6,'(16x,19(1H-),$)')
            IF (NSPIN.EQ.2) THEN
               WRITE (6,'(23(1H-))')
               WRITE (6,FMT=FMT2) ' TOT',
     &              (SUMCH(ISPIN),ISPIN=1,NSPIN),MUSPIN(LMAXD1+1)
               WRITE (6,'(33X,F14.8)') (SUMCH(1)+SUMCH(2))
            ELSE
               WRITE (6,*)
               WRITE (6,FMT=FMT2) ' TOT',SUMCH(1)
            END IF
C
         IF ( I1.NE.NAEZ ) THEN
            WRITE (6,'(9X,26(1H=),$)')
               IF (NSPIN.EQ.2) WRITE(6,'(23(1H=),$)')
               WRITE (6,*)
         END IF

C
      IF(I1.EQ.NAEZ) THEN
      WRITE (6,*)
      IF (NSPIN.EQ.2) THEN
         WRITE (6,'(78(1H#))')
      ELSE
         WRITE (6,'(44(1H#))')
      END IF
      WRITE (6,*)
      END IF
C
      RETURN
C
99001 FORMAT (15X,'l-decomposed valence charges and magnetic moments')
99002 FORMAT (8X,'l-decomposed valence charges')
99003 FORMAT (8X,' ATOM      ',$)
99004 FORMAT ('  N_el ',A7,$)
99005 FORMAT ('    m_spin',$)
      END
