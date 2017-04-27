      MODULE MOD_WRMOMS

      CONTAINS

      SUBROUTINE WRMOMS(NATOM,NSPIN, &
                        density,LMAXD,LMAXD1,LMAXATOM)

      USE nrtype
      USE type_density
      IMPLICIT NONE

! Dummy arguments
      INTEGER,intent(in)           ::  NATOM
      INTEGER,intent(in)           ::  NSPIN
      TYPE(DENSITY_TYPE)           ::  DENSITY(NATOM)
!       real(kind=DP),intent(in)     ::  CHARGE(0:LMAXD1,NATOM,2)
      INTEGER,intent(in)           ::  LMAXD
      INTEGER,intent(in)           ::  LMAXD1
      INTEGER,intent(in)           ::  LMAXATOM(NATOM)
!  Local variables
      CHARACTER(len=4),parameter   ::  TEXTL(0:6) = (/ ' s =',' p =',' d =',' f =',' g =',' h =',' i =' /)
      CHARACTER(len=7),parameter   ::  TEXTS(3)   = (/ 'spin dn','spin up','       ' /)
      CHARACTER(len=5),parameter   ::  TEXTNS     =    ' ns ='
      real(kind=DP)                ::  CHTOT(NATOM)
      real(kind=DP)                ::  MUSPIN(NATOM,0:LMAXD1+1), &
                                       SUMCH(NATOM,2)
      CHARACTER(len=80)            ::  FMT1,FMT2,FMT31,FMT32
      INTEGER                      ::  IS,ISPIN,IATOM,L,LF1,LF2

      WRITE (1337,*)
      WRITE (1337,'(44(1H#))')
      WRITE (1337,99002)
      WRITE (1337,'(44(1H#))')

      WRITE (1337,*)
      WRITE (1337,99003)
      DO IATOM = 1,NATOM
         MUSPIN(IATOM,LMAXD1+1) = 0D0
         SUMCH(IATOM,1) = 0D0
         SUMCH(IATOM,2) = 0D0
         DO L = 0,LMAXATOM(IATOM)+1 !LMAXD1
            DO ISPIN = 1,NSPIN
               SUMCH(IATOM,ISPIN) = SUMCH(IATOM,ISPIN) + DENSITY(IATOM)%NCHARGE(L,ISPIN)
            END DO
            IF ( NSPIN.EQ.2 ) THEN 
              MUSPIN(IATOM,L) = DENSITY(IATOM)%NCHARGE(L,2) - DENSITY(IATOM)%NCHARGE(L,1)
              MUSPIN(IATOM,LMAXD1+1) = MUSPIN(IATOM,LMAXD1+1) + MUSPIN(IATOM,L)
            END IF
         END DO
         CHTOT(IATOM) = SUMCH(IATOM,1) + SUMCH(IATOM,2)
      END DO

      IS = 0
      IF ( NSPIN.EQ.1 ) IS = IS + 2
      DO ISPIN = 1,NSPIN
         IS = IS + 1
         WRITE (1337,99004) TEXTS(IS)
      END DO

      IF (NSPIN.EQ.2) WRITE(1337,99005)
      WRITE(1337,*)

      WRITE (1337,'(3X,26(1H=),$)')
      IF (NSPIN.EQ.2) WRITE (1337,'(23(1H=),$)')
      WRITE(1337,*)

      FMT1 = '(4X,I3,2X,A4,2(F12.8),2X,F8.4'
      FMT2 = '(9X,A4,2(F12.8),2X,F8.4'
      FMT31 = '(4X,I3,2X,A4,F12.8)'
      FMT32 = '(9X,A4,F12.8)'
      LF1 = 30
      LF2 = 24

      IF (NSPIN.EQ.2) THEN
         FMT1 = FMT1(1:LF1)//')'
         FMT2 = FMT2(1:LF2)//')'
      ELSE
         FMT1 = FMT31
         FMT2 = FMT32
      END IF

      DO IATOM = 1,NATOM
         IF (NSPIN.EQ.2) THEN
            WRITE (1337,FMT=FMT1) IATOM,TEXTL(0), &
                (DENSITY(IATOM)%NCHARGE(0,ISPIN),ISPIN=1,NSPIN), &
                 MUSPIN(IATOM,0)
         ELSE
            WRITE (1337,FMT=FMT1) IATOM,TEXTL(0),DENSITY(IATOM)%NCHARGE(0,1)
         END IF

         DO L = 1,LMAXATOM(IATOM) !LMAXD
            IF (NSPIN.EQ.2) THEN
               WRITE (1337,FMT=FMT2) TEXTL(L), &
                    (DENSITY(IATOM)%NCHARGE(L,ISPIN),ISPIN=1,NSPIN), &
                    MUSPIN(IATOM,L)
            ELSE
               WRITE (1337,FMT=FMT2) TEXTL(L),DENSITY(IATOM)%NCHARGE(L,1)
            END IF

         END DO

         IF (NSPIN.EQ.2) THEN
            WRITE (1337,FMT=FMT2) TEXTNS, &
                 (DENSITY(IATOM)%NCHARGE(LMAXATOM(IATOM)+1,ISPIN),ISPIN=1,NSPIN), &
                 MUSPIN(IATOM,LMAXD1)
         ELSE
            WRITE (1337,FMT=FMT2) TEXTNS,DENSITY(IATOM)%NCHARGE(LMAXATOM(IATOM),1)
         END IF


         WRITE (1337,'(10x,19(1H-),$)')

         IF (NSPIN.EQ.2) THEN
            WRITE (1337,'(17(1H-))')
            WRITE (1337,FMT=FMT2) ' TOT', &
!                  (SUMCH(IATOM,ISPIN),ISPIN=1,NSPIN),MUSPIN(IATOM,LMAXATOM(IATOM)+1)
                 (SUMCH(IATOM,ISPIN),ISPIN=1,NSPIN),MUSPIN(IATOM,LMAXD1+1) !Phivos
            WRITE (1337,'(25X,F12.8)') CHTOT(IATOM)
         ELSE
            WRITE (1337,*)
            WRITE (1337,FMT=FMT2) ' TOT',SUMCH(IATOM,1)
         END IF


         IF ( IATOM.NE.NATOM ) THEN
            WRITE (1337,'(3X,26(1H=),$)')
            IF (NSPIN.EQ.2) WRITE(6,'(17(1H=),$)')
            WRITE (1337,*)
         END IF
      END DO

      WRITE (1337,*)
      WRITE (1337,'(44(1H#))')
      WRITE (1337,*)


99001 FORMAT (15X,'l-decomposed valence charges and magnetic moments')
99002 FORMAT (8X,'l-decomposed valence charges')
99003 FORMAT (3X,'ATOM      ',$)
99004 FORMAT (2X,'Ne ',A7,$)
99005 FORMAT ('    m_spin',$)
99006 FORMAT ('    m_orb   spin dn  spin up')
      END SUBROUTINE WRMOMS

      END MODULE MOD_WRMOMS
