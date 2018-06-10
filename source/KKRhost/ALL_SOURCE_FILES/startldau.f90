SUBROUTINE startldau(itrunldau,idoldau,kreadldau,lopt,ueff,jeff,  &
        erefldau,natyp,nspin,wldau,uldau,phildau,irws,  &
        ntldau,itldau,irmd,natypd,nspind,mmaxd)
! **********************************************************************
! *                                                                    *
! * Reads in LDA+U arrays from formatted file 'ldaupot'                *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!..
      INTEGER IRMD,MMAXD,NATYPD,NSPIND,IRWS(NATYPD)
!..
!.. Arguments ..
      INTEGER ITRUNLDAU,IDOLDAU,KREADLDAU,NATYP,NSPIN,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION UEFF(NATYPD),JEFF(NATYPD),EREFLDAU(NATYPD)
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD)
!      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD)
!.. 
!.. Locals ..
      INTEGER I1,IM1,IM3,IS,IT,LL
!     ..
! ----------------------------------------------------------------------


!      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

itrunldau = 0
idoldau = 1
ntldau = 0
DO it = 1,natyp
  IF ( lopt(it)+1 /= 0 ) THEN
    ntldau = ntldau + 1
    itldau(ntldau) = it
  END IF
END DO

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
WRITE(1337,'(79(1H=),/,27X,A,/, 79(1H=),/)') 'LDA+U: starting parameters'
WRITE(1337,99001) natyp,ntldau
WRITE(1337,99002)
WRITE(1337,99003)
DO it = 1,ntldau
  i1 = itldau(it)
  WRITE(1337,99004) i1,ueff(i1),jeff(i1),erefldau(i1)
END DO
WRITE(1337,99003)
WRITE(1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

! -> read in LDA+U from file if available (KREADLDAU=1)

CALL rinit(mmaxd*mmaxd*nspind*natypd,wldau)
CALL cinit(irmd*natypd,phildau)
IF ( kreadldau == 1 ) THEN
  WRITE(1337,99005)
  CALL readldaupot(itrunldau,lopt,ueff,jeff,  &
      erefldau,natyp,wldau,uldau,phildau,irws,  &
      ntldau,itldau,irmd,natypd,nspind,mmaxd)
ELSE
  WRITE(1337,99006)
END IF

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
IF ( itrunldau /= 0 ) THEN
  WRITE(1337,99007) 'Coulomb matrix U(m1,m1,m3,m3)'
  DO it = 1,ntldau
    i1 = itldau(it)
    ll = lopt(i1)
    ll = MIN(3,ll)
    WRITE(1337,99008) i1
    DO im1 = 1,2*ll+1
      WRITE(1337,99009) (uldau(im1,im1,im3,im3,i1),im3=1,2*ll+1)
    END DO
    WRITE(1337,*)
    IF ( it < ntldau ) WRITE(1337,99010)
  END DO
  WRITE(1337,99007) 'Interaction potential W(m1,m2)'
  DO it = 1,ntldau
    i1 = itldau(it)
    ll = lopt(i1)
    ll = MIN(3,ll)
    DO is = 1,nspin
      WRITE(1337,99011) i1,is
      DO im1 = 1,2*ll+1
        WRITE(1337,99009) (wldau(im1,im3,is,i1),im3=1,2*ll+1)
      END DO
      WRITE(1337,*)
    END DO
    IF ( it < ntldau ) WRITE(1337,99010)
  END DO
  WRITE(1337,'(9X,60(1H-))')
ELSE
  CALL rinit(mmaxd*mmaxd*mmaxd*mmaxd*natypd,uldau)
  CALL rinit(mmaxd*mmaxd*nspind*natypd,wldau)
  CALL cinit(irmd*natypd,phildau)
END IF
WRITE(1337,*)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

99001 FORMAT(6X,'Number of atoms ','  in the u.c. :',i4,/,  &
    24X,'using LDA+U :',i4,/)
99002 FORMAT(9X,' IT ','   Ueff   ','   Jeff   ','   Eref   ',' (Ry)')
99003 FORMAT(9X,40(1H-))
99004 FORMAT(9X,i3,1X,3F10.6)
99005 FORMAT(9X, 'Reading in LDA+U potential information (file ldaupot)',/)
99006 FORMAT(9X, 'LDA+U potential initialised (set to zero)')
99007 FORMAT(9X,60(1H-),/,9X,a,/,9X,60(1H-),/)
99008 FORMAT(9X,'IT =',i3,/)
99009 FORMAT(9X,7F10.6)
99010 FORMAT(11X,58(1H~))
99011 FORMAT(9X,'IT =',i3,' ISPIN =',i2,/)
END SUBROUTINE startldau
