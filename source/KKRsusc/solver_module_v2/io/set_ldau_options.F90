  subroutine set_ldau_options(itc,my_rank)
! global options read in from input file
  use global

  implicit none


! Current iteration
  integer(kind=i4b), intent(in) :: itc, my_rank
! ----------------------------------------------------------------------
! line number, column number
  integer(kind=i4b) :: iline, ipos
! was key found?
  logical           :: found
! ----------------------------------------------------------------------
  integer(kind=i4b) :: ig, istart, iend, i, ntot
  real(kind=r8b)    :: utemp(0:nlmax)


! Not an LDA+U calculation
  if (.not.lldau) then
    if (my_rank == 0) then 
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No LDA+U correction")')
        write(*,'("************************************************************"/)')
      end if 
    end if ! my_rank
    lrhomat = .false.
    return
  end if

! LDA+U options

! read density matrix from file
  call find_keyinfile('lrhomat',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lrhomat
  else 
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_ldau_options: WARNING - key lrhomat not found!")')
    end if ! my_rank 
    lrhomat = .false.
  end if

! mix factor for scf update of density matrix
  call find_keyinfile('ldaumix',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    if (lwarn .and. itc == 1) read(inputfile(iline)(ipos:nchars),*) ldaumix
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_ldau_options: WARNING - key ldaumix not found!")')
    end if ! my_rank
    ldaumix = 0.05d0
  end if

! after which scf iteration to start mixing density matrix
  call find_keyinfile('ldauitc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ldauitc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_ldau_options: WARNING - key ldauitc not found!")')
    end if ! my_rank
    ldauitc = 20
  end if

! Defaults
  ildau(1:nasusc) = 0; ntot = 0
  ueff(0:nlmax,1:nasusc) = 0.d0; jeff(0:nlmax,1:nasusc) = 0.d0

! Search for atom info lines
! first pass: which atoms to apply LDA+U to
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     which atoms are in group ig
      call find_keyinline('ildau',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) i
        if (i /= 0) then
          if (i /= 1 .and. i /= 2) stop 'set_ldau_options: unknown ildau option'
          ildau(istart:iend) = i
          ntot = ntot + igroup(ig)
        end if
      end if
    end if
  end do
! second pass: read remaining options
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     which atoms are in group ig
      if (ildau(istart) /= 0) then
!       U parameters
        call find_keyinline('ueff',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) utemp(0:nlmax)
          ueff(0:nlmax,istart:iend) = spread(utemp(0:nlmax),2,igroup(ig))
        end if
!       J parameters
        call find_keyinline('jeff',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) utemp(0:nlmax)
          jeff(0:nlmax,istart:iend) = spread(utemp(0:nlmax),2,igroup(ig))
        end if
      end if
    end if
  end do

! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" LDA+U correction  ==>  check atomic options")')
!     to start from an existing density matrix for LDA+U or not
      if (lrhomat) then
        write(*,'(" + Density matrix read in from file")')
      else
        write(*,'(" + Density matrix computed internally")')
      end if
!     From which iteration to mix density matrices, mixing factor
      write(*,'(" + Density matrix mixed after iteration=",i4," with mixing=",f8.4)') ldauitc, ldaumix
!     Atom info
      write(*,'(/" Number of atoms for LDA+U calculation:",i4)') ntot
      write(*,'(/" LDA+U configuration for each atom group:")')
      do ig=1,ngroup
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
        if (ildau(istart) /= 0) then
          write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
          if (ildau(istart) == 1) write(*,'(" + LDA+U with Dudarev U-J=",5f8.4)') ueff(:,istart) - jeff(:,istart)
          if (ildau(istart) == 2) write(*,'(" + LDA+U with full Coulomb interaction U, J=",10f8.4)') ueff(:,istart), jeff(:,istart)
        end if
      end do
      write(*,'(" WARNING: be mindful of vldaushift computed in read_rhomat and optionally added in build_vldau!!!")')
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank
! All done!
  end subroutine set_ldau_options
