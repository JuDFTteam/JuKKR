      SUBROUTINE WRLDAUPOT(ITRUNLDAU,LOPT,UEFF,JEFF,
     &                     EREFLDAU,NATYP,WLDAU,ULDAU,PHILDAU,
     &                     IRMD,NATYPD,NSPIND,MMAXD,IRWS)
C **********************************************************************
C *                                                                    *
C * Writes out LDA+U arrays into formatted file 'ldaupot'              *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
      INTEGER IRMD,MMAXD,NATYPD,NSPIND,IRWS(NATYPD)
C     ..
C     .. Arguments ..
      INTEGER ITRUNLDAU,NATYP
      INTEGER LOPT(NATYPD)
      DOUBLE PRECISION UEFF(NATYPD),JEFF(NATYPD),EREFLDAU(NATYPD)
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
C      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD)
C     ..
C     ..  Locals 
      INTEGER IR,M1,M2,M3,M4,IT,IS
C
C ======================================================================


C      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )
      
      
      OPEN (67,FILE='ldaupot_new',FORM='FORMATTED')
      WRITE(1337,99001)
      WRITE(67,99002) ITRUNLDAU,'    ITRUNLDAU'
      WRITE(67,99002) NATYP,'    NATYP'
      WRITE(67,99003) NATYP
      WRITE(67,99004) (LOPT(IT),IT=1,NATYP)
      WRITE(67,99005) 
      DO IT = 1,NATYP
         IF ( LOPT(IT)+1.NE.0 ) WRITE(67,99006) 
     &         IT,UEFF(IT),JEFF(IT),EREFLDAU(IT)
      END DO
      DO IT = 1,NATYP
         IF ( LOPT(IT)+1.NE.0 ) THEN
            WRITE(67,99002) IT,'    WLDAU'
            DO IS = 1,NSPIND
               DO M1 = 1,MMAXD
                  WRITE (67,99007) (WLDAU(M1,M2,IS,IT),M2=1,MMAXD)
               END DO
            END DO
            WRITE(67,99002) IT,'    ULDAU'
            WRITE (67,99007) ((((ULDAU(M1,M2,M3,M4,IT),M4=1,MMAXD),
     &           M3=1,MMAXD),M2=1,MMAXD),M1=1,MMAXD)
            WRITE(67,99002) IT,'    PHILDAU'
            WRITE (67,99007) (PHILDAU(IR,IT),IR=1,IRWS(IT))
         END IF
      END DO
      CLOSE (67)
C
99001 FORMAT(/,5X,'< WRLDAUPOT > : ',
     &     'Writing out LDA+U potential (file ldaupot_new)',/)
99002 FORMAT(I6,A)
99003 FORMAT('LOPT 1..',I3)
99004 FORMAT(16I3)
99005 FORMAT('IAT',6X,'UEFF',12X,'JEFF',12X,'EREF')
99006 FORMAT(I3,3(1X,E15.8))
99007 FORMAT(5E16.8)
      END
