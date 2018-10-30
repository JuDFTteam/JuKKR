      SUBROUTINE STARTLDAU(ITRUNLDAU,IDOLDAU,KREADLDAU,LOPT,UEFF,JEFF,
     &                    EREFLDAU,NATYP,NSPIN,WLDAU,ULDAU,PHILDAU,IRWS,
     &                     NTLDAU,ITLDAU,IRMD,NATYPD,NSPIND,MMAXD)
C **********************************************************************
C *                                                                    *
C * Reads in LDA+U arrays from formatted file 'ldaupot'                *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
      INTEGER IRMD,MMAXD,NATYPD,NSPIND,IRWS(NATYPD)
C     ..
C     .. Arguments ..
      INTEGER ITRUNLDAU,IDOLDAU,KREADLDAU,NATYP,NSPIN,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION UEFF(NATYPD),JEFF(NATYPD),EREFLDAU(NATYPD)
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD)
c      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD)
C     .. 
C     .. Locals ..
      INTEGER I1,IM1,IM3,IS,IT,LL
C     ..
C ----------------------------------------------------------------------


!      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

      ITRUNLDAU = 0
      IDOLDAU = 1
      NTLDAU = 0
      DO IT = 1,NATYP
         IF ( LOPT(IT)+1.NE.0 ) THEN
            NTLDAU = NTLDAU + 1
            ITLDAU(NTLDAU) = IT
         END IF
      END DO
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      WRITE(1337,'(79(1H=),/,27X,A,/,
     &          79(1H=),/)') 'LDA+U: starting parameters'
      WRITE(1337,99001) NATYP,NTLDAU
      WRITE(1337,99002) 
      WRITE(1337,99003)
      DO IT = 1,NTLDAU
         I1 = ITLDAU(IT)
         WRITE(1337,99004) I1,UEFF(I1),JEFF(I1),EREFLDAU(I1)
      END DO
      WRITE(1337,99003)
      WRITE(1337,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
C -> read in LDA+U from file if available (KREADLDAU=1)
C
      CALL RINIT(MMAXD*MMAXD*NSPIND*NATYPD,WLDAU)
      CALL CINIT(IRMD*NATYPD,PHILDAU)
      IF ( KREADLDAU.EQ.1 ) THEN
         WRITE(1337,99005)
         CALL READLDAUPOT(ITRUNLDAU,LOPT,UEFF,JEFF,
     &                    EREFLDAU,NATYP,WLDAU,ULDAU,PHILDAU,IRWS,
     &                    NTLDAU,ITLDAU,IRMD,NATYPD,NSPIND,MMAXD)
      ELSE
         WRITE(1337,99006)
      END IF
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      IF ( ITRUNLDAU.NE.0 ) THEN
         WRITE(1337,99007) 'Coulomb matrix U(m1,m1,m3,m3)'
         DO IT = 1,NTLDAU
            I1 = ITLDAU(IT)
            LL = LOPT(I1)
            LL = MIN(3,LL)
            WRITE(1337,99008) I1
            DO IM1 = 1,2*LL+1
               WRITE(1337,99009) 
     &                   (ULDAU(IM1,IM1,IM3,IM3,I1),IM3=1,2*LL+1)
            END DO
            WRITE(1337,*)
            IF ( IT.LT.NTLDAU ) WRITE(1337,99010)
         END DO
         WRITE(1337,99007) 'Interaction potential W(m1,m2)'
         DO IT = 1,NTLDAU
            I1 = ITLDAU(IT)
            LL = LOPT(I1)
            LL = MIN(3,LL)
            DO IS = 1,NSPIN
               WRITE(1337,99011) I1,IS
               DO IM1 = 1,2*LL+1
                  WRITE(1337,99009) (WLDAU(IM1,IM3,IS,I1),IM3=1,2*LL+1)
               END DO
               WRITE(1337,*)
            END DO
            IF ( IT.LT.NTLDAU ) WRITE(1337,99010)
         END DO
         WRITE(1337,'(9X,60(1H-))')
      ELSE
         CALL RINIT(MMAXD*MMAXD*MMAXD*MMAXD*NATYPD,ULDAU)
         CALL RINIT(MMAXD*MMAXD*NSPIND*NATYPD,WLDAU)
         CALL CINIT(IRMD*NATYPD,PHILDAU)
      END IF
      WRITE(1337,*)

C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
99001 FORMAT(6X,'Number of atoms ','  in the u.c. :',I4,/,
     &       24X,'using LDA+U :',I4,/)
99002 FORMAT(9X,' IT ','   Ueff   ','   Jeff   ','   Eref   ',' (Ry)')
99003 FORMAT(9X,40(1H-))
99004 FORMAT(9X,I3,1X,3F10.6)
99005 FORMAT(9X,
     &     'Reading in LDA+U potential information (file ldaupot)',/)
99006 FORMAT(9X,
     &     'LDA+U potential initialised (set to zero)')
99007 FORMAT(9X,60(1H-),/,9X,A,/,9X,60(1H-),/)
99008 FORMAT(9X,'IT =',I3,/)
99009 FORMAT(9X,7F10.6)
99010 FORMAT(11X,58(1H~))
99011 FORMAT(9X,'IT =',I3,' ISPIN =',I2,/)
      END
