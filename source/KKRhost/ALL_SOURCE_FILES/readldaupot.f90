module mod_readldaupot

contains

subroutine readldaupot(itrunldau, lopt, ueff, jeff, erefldau, natyp, wldau, &
  uldau, phildau, irws, ntldau, itldau, irmd, natypd, nspind, mmaxd)
  ! **********************************************************************
  ! *                                                                    *
  ! * Reads in LDA+U arrays from formatted file 'ldaupot'                *
  ! *                                                                    *
  ! **********************************************************************

  use :: mod_version_info
  use :: mod_datatypes, only: dp
  implicit none
  ! ..
  integer :: irmd, mmaxd, natypd, nspind, irws(natypd)
  ! ..
  ! .. Arguments ..
  integer :: itrunldau, natyp, ntldau
  integer :: lopt(natypd), itldau(natypd)
  real (kind=dp) :: ueff(natypd), jeff(natypd), erefldau(natypd)
  real (kind=dp) :: wldau(mmaxd, mmaxd, nspind, natypd)
  real (kind=dp) :: uldau(mmaxd, mmaxd, mmaxd, mmaxd, natypd)
  ! ..OUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:)
  complex (kind=dp) :: phildau(irmd, natypd)
  ! ..
  ! ..  Locals
  integer :: ios, ir, m1, m2, m3, m4, it, i1, i2, is
  integer :: irunldau, ntloc
  integer :: loptldau(natypd)
  real (kind=dp) :: ueff0, jeff0, eref0
  ! ======================================================================


  ! ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

  open (67, file='ldaupot', form='FORMATTED', status='OLD', iostat=ios)
  call version_check_header(67)
  if (ios>0) then
    write (6, 110) 'Could not find LDA+U file'
    itrunldau = 0
    return
  end if
  ! ======================================================================
  ! -> READ IN : itrunldau, natyp

  read (67, *, err=100) irunldau
  read (67, *, err=100) ntloc
  if (ntloc/=natyp) then
    close (67)
    write (6, 120) 'Inconsistent NATYP value in LDA+U file'
    itrunldau = 0
    return
  end if
  read (67, *, err=100)
  ! ======================================================================
  ! -> READ IN : lopt(1..natyp) - set NT = no. of atoms lda+u treated

  read (67, *, err=100)(loptldau(i2), i2=1, natyp)
  do i2 = 1, natyp
    if (loptldau(i2)/=lopt(i2)) then
      close (67)
      write (6, 120) 'Inconsistent LOPT values in LDA+U file'
      itrunldau = 0
      return
    end if
  end do
  ! ======================================================================
  ! -> READ IN : ueff,jeff,erefldau for the NTLDAU atoms

  read (67, *, err=100)
  do it = 1, ntldau
    read (67, *, err=100) i2, ueff0, jeff0, eref0
    i1 = 0
    do ir = 1, ntldau
      if (i2==itldau(ir)) i1 = 1
    end do
    if (i1==0) then
      close (67)
      write (6, 120) 'Inconsistent UEFF/JEFF/EREF values in LDA+U file'
      itrunldau = 0
      return
    end if
    ueff0 = abs(ueff0-ueff(i2))
    jeff0 = abs(jeff0-jeff(i2))
    eref0 = abs(eref0-erefldau(i2))
    if ((ueff0>1d-8) .or. (ueff0>1d-8) .or. (eref0>1d-8)) then
      close (67)
      write (6, 120) 'Inconsistent UEFF/JEFF/EREF values in LDA+U file'
      itrunldau = 0
      return
    end if
  end do
  ! ======================================================================
  ! -> READ IN : wldau,uldau for the NTLDAU atoms

  do it = 1, ntldau
    read (67, *, err=100) i2
    i1 = 0
    do ir = 1, ntldau
      if (i2==itldau(ir)) i1 = 1
    end do
    if (i1==0) then
      close (67)
      write (6, 110) 'Inconsistent WLDAU/ULDAU in LDA+U file'
      itrunldau = 0
      return
    end if
    ! ---------------------------------------------------------------- WLDAU
    do is = 1, nspind
      do m1 = 1, mmaxd
        read (67, *, iostat=ios)(wldau(m1,m2,is,i2), m2=1, mmaxd)
        if (ios/=0) then
          write (6, 110) 'Corrupted WLDAU array in LDA+U file'
          close (67)
          itrunldau = 0
          return
        end if
      end do
    end do
    ! ---------------------------------------------------------------- ULDAU
    read (67, *, err=100)

    read (67, *, iostat=ios)((((uldau(m1,m2,m3,m4,i2),m4=1,mmaxd),m3=1, &
      mmaxd),m2=1,mmaxd), m1=1, mmaxd)
    if (ios/=0) then
      write (6, 110) 'Corrupted ULDAU array in LDA+U file'
      close (67)
      itrunldau = 0
      return
    end if

    ! DO M1 = 1,MMAXD
    ! DO M2 = 1,MMAXD
    ! DO M3 = 1,MMAXD
    ! READ(67,*,IOSTAT=IOS)
    ! &                 (ULDAU(M1,M2,M3,M4,I2),M4=1,MMAXD)
    ! IF ( IOS.NE.0 ) THEN
    ! WRITE(6,99001)
    ! &                    'Corrupted ULDAU array in LDA+U file'
    ! CLOSE(67)
    ! ITRUNLDAU = 0
    ! RETURN
    ! END IF
    ! END DO
    ! END DO
    ! END DO
    ! ----------------------------------------------------------------------
    ! END DO
    ! ======================================================================
    ! -> READ IN : phildau

    ! DO IT = 1,NTLDAU
    read (67, *, err=100) i2
    i1 = 0
    do ir = 1, ntldau
      if (i2==itldau(ir)) i1 = 1
    end do
    if (i1==0) then
      close (67)
      write (6, 110) 'Inconsistent PHILDAU values in LDA+U file'
      itrunldau = 0
      return
    end if
    read (67, '(5E16.8)', iostat=ios)(phildau(i1,i2), i1=1, irws(i2))
    if (ios/=0) then
      write (6, 110) 'Corrupted PHILDAU array in LDA+U file '
      close (67)
      itrunldau = 0
      return
    end if
  end do
  ! ======================================================================

  if (irunldau==0) write (6, 120) &
    'ITRUNLDAU=0 found in the (otherwise consistent) LDA+U file'
  itrunldau = irunldau
  close (67)
  return

100 write (6, 110) 'Problems reading in LDA+U file'
  itrunldau = 0
  close (67)
110 format (9x, 'WARNING: ', a, /, 18x, &
    'LDA+U potentials set to zero, iteration reinitialised')
120 format (9x, 'WARNING: ', a, /, 18x, &
    'input-card data will be used, iteration reinitialised')
end subroutine readldaupot

end module mod_readldaupot
