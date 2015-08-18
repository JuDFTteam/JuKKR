      SUBROUTINE READLDAUPOT(ITRUNLDAU,LOPT,UEFF,JEFF,
     &                       EREFLDAU,NATYP,WLDAU,ULDAU,PHILDAU,
     &                       NTLDAU,ITLDAU,IRMD,NATYPD,NSPIND,MMAXD)
C **********************************************************************
C *                                                                    *
C * Reads in LDA+U arrays from formatted file 'ldaupot'                *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
      INTEGER IRMD,MMAXD,NATYPD,NSPIND
C     ..
C     .. Arguments ..
      INTEGER ITRUNLDAU,NATYP,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION UEFF(NATYPD),JEFF(NATYPD),EREFLDAU(NATYPD)
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
!       DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD)
      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD)
C     ..
C     ..  Locals 
      INTEGER IOS,IR,M1,M2,M3,M4,IT,I1,I2,IS
      INTEGER IRUNLDAU,NTLOC
      INTEGER LOPTLDAU(NATYPD)
      DOUBLE PRECISION UEFF0,JEFF0,EREF0
C
C ======================================================================


      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

      OPEN (67,FILE='ldaupot',FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
      IF ( IOS.GT.0 ) THEN
         WRITE(6,99001) 'Could not find LDA+U file'
         ITRUNLDAU = 0
         RETURN
      END IF
C ======================================================================
C -> READ IN : itrunldau, natyp
C
      READ (67,*,ERR=99100) IRUNLDAU
      READ (67,*,ERR=99100) NTLOC
      IF ( NTLOC.NE.NATYP ) THEN
         CLOSE (67)
         WRITE(6,99002) 'Inconsistent NATYP value in LDA+U file'
         ITRUNLDAU = 0
         RETURN
      END IF
      READ (67,*,ERR=99100) 
C ======================================================================
C -> READ IN : lopt(1..natyp) - set NT = no. of atoms lda+u treated
C
      READ (67,*,ERR=99100) (LOPTLDAU(I2),I2=1,NATYP)
      DO I2 = 1,NATYP
         IF ( LOPTLDAU(I2).NE.LOPT(I2) ) THEN
            CLOSE (67)
            WRITE(6,99002)
     &           'Inconsistent LOPT values in LDA+U file'
            ITRUNLDAU = 0
            RETURN
         END IF
      END DO
C ======================================================================
C -> READ IN : ueff,jeff,erefldau for the NTLDAU atoms
C
      READ (67,*,ERR=99100)
      DO IT = 1,NTLDAU
         READ (67,*,ERR=99100) I2,UEFF0,JEFF0,EREF0
         I1 = 0
         DO IR = 1,NTLDAU
            IF ( I2.EQ.ITLDAU(IR) ) I1 = 1
         END DO
         IF ( I1.EQ.0 ) THEN
            CLOSE (67)
            WRITE(6,99002)
     &           'Inconsistent UEFF/JEFF/EREF values in LDA+U file'
            ITRUNLDAU = 0
            RETURN
         END IF
         UEFF0 = DABS( UEFF0-UEFF(I2) )
         JEFF0 = DABS( JEFF0-JEFF(I2) )
         EREF0 = DABS( EREF0-EREFLDAU(I2) )
         IF ( ( UEFF0.GT.1D-8 ) .OR. ( UEFF0.GT.1D-8 ) .OR.
     &        ( EREF0.GT.1D-8 ) ) THEN
            CLOSE (67)
            WRITE(6,99002)
     &           'Inconsistent UEFF/JEFF/EREF values in LDA+U file'
            ITRUNLDAU = 0
            RETURN
         END IF
      END DO
C ======================================================================
C -> READ IN : wldau,uldau for the NTLDAU atoms
C
      DO IT = 1,NTLDAU
         READ (67,*,ERR=99100) I2
         I1 = 0
         DO IR = 1,NTLDAU
            IF ( I2.EQ.ITLDAU(IR) ) I1 = 1
         END DO
         IF ( I1.EQ.0 ) THEN
            CLOSE (67)
            WRITE(6,99001) 'Inconsistent WLDAU/ULDAU in LDA+U file'
            ITRUNLDAU = 0
            RETURN
         END IF
C ---------------------------------------------------------------- WLDAU
         DO IS = 1,NSPIND
            DO M1 = 1,MMAXD
               READ(67,*,IOSTAT=IOS) (WLDAU(M1,M2,IS,I2),M2=1,MMAXD)
               IF ( IOS.NE.0 ) THEN
                  WRITE(6,99001) 'Corrupted WLDAU array in LDA+U file'
                  CLOSE (67)
                  ITRUNLDAU = 0
                  RETURN
               END IF
            END DO
         END DO
C ---------------------------------------------------------------- ULDAU
         READ (67,*,ERR=99100)

         READ(67,*,IOSTAT=IOS)
     &        ((((ULDAU(M1,M2,M3,M4,I2),M4=1,MMAXD),M3=1,MMAXD)
     &        ,M2=1,MMAXD),M1=1,MMAXD)
         IF ( IOS.NE.0 ) THEN
            WRITE(6,99001)
     &           'Corrupted ULDAU array in LDA+U file'
            CLOSE(67)
            ITRUNLDAU = 0
            RETURN
         END IF

c        DO M1 = 1,MMAXD
c           DO M2 = 1,MMAXD
c              DO M3 = 1,MMAXD 
c                 READ(67,*,IOSTAT=IOS)
c    &                 (ULDAU(M1,M2,M3,M4,I2),M4=1,MMAXD)
c                 IF ( IOS.NE.0 ) THEN
c                    WRITE(6,99001)
c    &                    'Corrupted ULDAU array in LDA+U file'
c                    CLOSE(67)
c                    ITRUNLDAU = 0
c                    RETURN
c                 END IF
c              END DO
c           END DO
c        END DO
C ----------------------------------------------------------------------
c     END DO
C ======================================================================
C -> READ IN : phildau
C 
c     DO IT = 1,NTLDAU
         READ (67,*,ERR=99100) I2
         I1 = 0
         DO IR = 1,NTLDAU
            IF ( I2.EQ.ITLDAU(IR) ) I1 = 1
         END DO
         IF ( I1.EQ.0 ) THEN
            CLOSE (67)
            WRITE(6,99001) 'Inconsistent PHILDAU values in LDA+U file'
            ITRUNLDAU = 0
            RETURN
         END IF
         READ(67,'(5E16.8)',IOSTAT=IOS) (PHILDAU(I1,I2),I1=1,IRMD)
         IF ( IOS.NE.0 ) THEN
            WRITE(6,99001) 'Corrupted PHILDAU array in LDA+U file '
            CLOSE(67)
            ITRUNLDAU = 0
            RETURN
         END IF
      END DO
C ======================================================================
C
      IF ( IRUNLDAU.EQ.0 ) WRITE(6,99002) 
     &     'ITRUNLDAU=0 found in the (otherwise consistent) LDA+U file'
      ITRUNLDAU = IRUNLDAU
      CLOSE (67)
      RETURN
C
99100 WRITE(6,99001) 'Problems reading in LDA+U file'
      ITRUNLDAU = 0
      CLOSE (67)
99001 FORMAT(9X,'WARNING: ',A,/,18X,
     &     'LDA+U potentials set to zero, iteration reinitialised')
99002 FORMAT(9X,'WARNING: ',A,/,18X,
     &     'input-card data will be used, iteration reinitialised')
      END
