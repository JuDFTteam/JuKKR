SUBROUTINE readldaupot(itrunldau,lopt,ueff,jeff,  &
        erefldau,natyp,wldau,uldau,phildau,irws,  &
        ntldau,itldau,irmd,natypd,nspind,mmaxd)
! **********************************************************************
! *                                                                    *
! * Reads in LDA+U arrays from formatted file 'ldaupot'                *
! *                                                                    *
! **********************************************************************

      use mod_version_info
      IMPLICIT NONE
!..
      INTEGER IRMD,MMAXD,NATYPD,NSPIND,IRWS(NATYPD)
!..
!.. Arguments ..
      INTEGER ITRUNLDAU,NATYP,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION UEFF(NATYPD),JEFF(NATYPD),EREFLDAU(NATYPD)
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD)
!..OUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD)
!..
!..  Locals 
      INTEGER IOS,IR,M1,M2,M3,M4,IT,I1,I2,IS
      INTEGER IRUNLDAU,NTLOC
      INTEGER LOPTLDAU(NATYPD)
      DOUBLE PRECISION UEFF0,JEFF0,EREF0
! ======================================================================


!      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

OPEN (67,FILE='ldaupot',FORM='FORMATTED',STATUS='OLD',IOSTAT=ios)
CALL version_check_header(67)
IF ( ios > 0 ) THEN
  WRITE(6,99001) 'Could not find LDA+U file'
  itrunldau = 0
  RETURN
END IF
! ======================================================================
! -> READ IN : itrunldau, natyp

READ (67,*,ERR=99100) irunldau
READ (67,*,ERR=99100) ntloc
IF ( ntloc /= natyp ) THEN
  CLOSE (67)
  WRITE(6,99002) 'Inconsistent NATYP value in LDA+U file'
  itrunldau = 0
  RETURN
END IF
READ (67,*,ERR=99100)
! ======================================================================
! -> READ IN : lopt(1..natyp) - set NT = no. of atoms lda+u treated

READ (67,*,ERR=99100) (loptldau(i2),i2=1,natyp)
DO i2 = 1,natyp
  IF ( loptldau(i2) /= lopt(i2) ) THEN
    CLOSE (67)
    WRITE(6,99002) 'Inconsistent LOPT values in LDA+U file'
    itrunldau = 0
    RETURN
  END IF
END DO
! ======================================================================
! -> READ IN : ueff,jeff,erefldau for the NTLDAU atoms

READ (67,*,ERR=99100)
DO it = 1,ntldau
  READ (67,*,ERR=99100) i2,ueff0,jeff0,eref0
  i1 = 0
  DO ir = 1,ntldau
    IF ( i2 == itldau(ir) ) i1 = 1
  END DO
  IF ( i1 == 0 ) THEN
    CLOSE (67)
    WRITE(6,99002) 'Inconsistent UEFF/JEFF/EREF values in LDA+U file'
    itrunldau = 0
    RETURN
  END IF
  ueff0 = DABS( ueff0-ueff(i2) )
  jeff0 = DABS( jeff0-jeff(i2) )
  eref0 = DABS( eref0-erefldau(i2) )
  IF ( ( ueff0 > 1D-8 ) .OR. ( ueff0 > 1D-8 ) .OR. ( eref0 > 1D-8 ) ) THEN
    CLOSE (67)
    WRITE(6,99002) 'Inconsistent UEFF/JEFF/EREF values in LDA+U file'
    itrunldau = 0
    RETURN
  END IF
END DO
! ======================================================================
! -> READ IN : wldau,uldau for the NTLDAU atoms

DO it = 1,ntldau
  READ (67,*,ERR=99100) i2
  i1 = 0
  DO ir = 1,ntldau
    IF ( i2 == itldau(ir) ) i1 = 1
  END DO
  IF ( i1 == 0 ) THEN
    CLOSE (67)
    WRITE(6,99001) 'Inconsistent WLDAU/ULDAU in LDA+U file'
    itrunldau = 0
    RETURN
  END IF
! ---------------------------------------------------------------- WLDAU
  DO is = 1,nspind
    DO m1 = 1,mmaxd
      READ(67,*,IOSTAT=ios) (wldau(m1,m2,is,i2),m2=1,mmaxd)
      IF ( ios /= 0 ) THEN
        WRITE(6,99001) 'Corrupted WLDAU array in LDA+U file'
        CLOSE (67)
        itrunldau = 0
        RETURN
      END IF
    END DO
  END DO
! ---------------------------------------------------------------- ULDAU
  READ (67,*,ERR=99100)
  
  READ(67,*,IOSTAT=ios) ((((uldau(m1,m2,m3,m4,i2),m4=1,mmaxd),m3=1,mmaxd)  &
      ,m2=1,mmaxd),m1=1,mmaxd)
  IF ( ios /= 0 ) THEN
    WRITE(6,99001) 'Corrupted ULDAU array in LDA+U file'
    CLOSE(67)
    itrunldau = 0
    RETURN
  END IF
  
!        DO M1 = 1,MMAXD
!           DO M2 = 1,MMAXD
!              DO M3 = 1,MMAXD
!                 READ(67,*,IOSTAT=IOS)
!    &                 (ULDAU(M1,M2,M3,M4,I2),M4=1,MMAXD)
!                 IF ( IOS.NE.0 ) THEN
!                    WRITE(6,99001)
!    &                    'Corrupted ULDAU array in LDA+U file'
!                    CLOSE(67)
!                    ITRUNLDAU = 0
!                    RETURN
!                 END IF
!              END DO
!           END DO
!        END DO
! ----------------------------------------------------------------------
!     END DO
! ======================================================================
! -> READ IN : phildau
  
!     DO IT = 1,NTLDAU
  READ (67,*,ERR=99100) i2
  i1 = 0
  DO ir = 1,ntldau
    IF ( i2 == itldau(ir) ) i1 = 1
  END DO
  IF ( i1 == 0 ) THEN
    CLOSE (67)
    WRITE(6,99001) 'Inconsistent PHILDAU values in LDA+U file'
    itrunldau = 0
    RETURN
  END IF
  READ(67,'(5E16.8)',IOSTAT=ios) (phildau(i1,i2),i1=1,irws(i2))
  IF ( ios /= 0 ) THEN
    WRITE(6,99001) 'Corrupted PHILDAU array in LDA+U file '
    CLOSE(67)
    itrunldau = 0
    RETURN
  END IF
END DO
! ======================================================================

IF ( irunldau == 0 ) WRITE(6,99002)  &
    'ITRUNLDAU=0 found in the (otherwise consistent) LDA+U file'
itrunldau = irunldau
CLOSE (67)
RETURN

99100 WRITE(6,99001) 'Problems reading in LDA+U file'
itrunldau = 0
CLOSE (67)
99001 FORMAT(9X,'WARNING: ',a,/,18X,  &
    'LDA+U potentials set to zero, iteration reinitialised')
99002 FORMAT(9X,'WARNING: ',a,/,18X,  &
    'input-card data will be used, iteration reinitialised')
END SUBROUTINE readldaupot
