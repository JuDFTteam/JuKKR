  subroutine set_global_options(lmax,natyp,nspin,nrad,itc,my_rank)
! global options read in from input file
  use global

  implicit none


! Values given in the main program
  integer(kind=i4b), intent(in) :: lmax, natyp, nspin, nrad, itc, my_rank
! -----------------------------------------------------------------------
! line number, column number
  integer(kind=i4b) :: iline, ipos
! was key found?
  logical           :: found
! -----------------------------------------------------------------------
  integer(kind=i4b) :: ig, istart, iend, ntot, iw(0:lmax)
  logical           :: lbasisdefault


! Warning about flags not being found
  call find_keyinfile('lwarn',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lwarn
  else
    lwarn = .false.
  end if

! General calculation flags
! Which KKR host program
  call find_keyinfile('ikkr',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ikkr
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key ikkr not found!")')
    end if ! my_rank
    ikkr = 2
  end if
  if (ikkr /= 1 .and. ikkr /= 2)  stop 'Unknown KKR program'

! Compute DOS?
  call find_keyinfile('ldos',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ldos
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key ldos not found!")')
    end if ! my_rank 
    ldos = .false.
  end if

! Add SOC?
  call find_keyinfile('lsoc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lsoc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lsoc not found!")')
    end if ! my_rank
    lsoc = .false.
  end if

! Add LDAU?
  call find_keyinfile('lldau',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lldau
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key ldau not found!")')
    end if ! my_rank
    lldau = .false.
  end if

! Add vector B field?
  call find_keyinfile('lbfield',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lbfield
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lbfield not found!")')
    end if ! my_rank 
   lbfield = .false.
  end if

! Compute susceptibilities?
  call find_keyinfile('lsusc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lsusc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lsusc not found!")')
    end if ! my_rank
    lsusc = .false.
  end if

! Perform spin rotations?
  call find_keyinfile('lrot',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lrot
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lrot not found!")')
    end if ! my_rank
    lrot = .false.
  end if

! Rational fit of GF?
  call find_keyinfile('lfit',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lfit
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lfit not found!")')
    end if ! my_rank
    lfit = .false.
  end if

! Compute magnetic couplings?
  call find_keyinfile('ljij',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ljij
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key ljij not found!")')
    end if ! my_rank
    ljij = .false.
  end if

! Read in model self-energy?
  call find_keyinfile('lsemodel',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lsemodel
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lsemodel not found!")')
    end if ! my_rank
    lsemodel = .false.
  end if

! Write tons of files on HDD?
  call find_keyinfile('lhdio',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lhdio
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lhdio not found!")')
    end if ! my_rank
    lhdio = .false.
  end if

! Restart from outsusc.dat
  call find_keyinfile('lrestart',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lrestart
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lrestart not found!")')
    end if ! my_rank
    lrestart = .false.
  end if

! Run in serial mode (like in the good old version without Jubas changes)
  call find_keyinfile('lserial',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) lserial
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key lserial not found!")')
    end if ! my_rank
    lserial = .true.
  end if

! ----------------------------------------------------------------------
  if (my_rank == 0) then 
    if (itc == 1) then
      write(*,'(/,"************************************************************")')
!     Which version of the KKR program
      if (ikkr == 1) then
        write(*,'(" KKR program: old")')
      else if (ikkr == 2) then
        write(*,'(" KKR program: KKRFLEX")')
      end if
!     Disk I/O option
      if (lhdio) then
        write(*,'(" Disk I/O used for wfns, pots, outsusc.dat, etc")')
      else
        write(*,'(" Disk I/O not used")')
      end if
      write(*,'("************************************************************",/)')
!     Restart mode reads from oususc.dat
      if (lrestart) then
        write(*,'(" Read from outsusc.dat the projection coeffiecients")')
      else
        write(*,'(" Nothing read outsusc.dat")')
      end if
      write(*,'("************************************************************",/)')
    end if
  end if ! my_rank
! ----------------------------------------------------------------------
! These must be given in the input file:

! How many energies
  call find_keyinfile('nesusc',nchars,nlines,inputfile,iline,ipos,found)
  if (.not.found) then
    stop 'set_global_options: key nesusc not found!'
  else
    read(inputfile(iline)(ipos:nchars),*) nesusc
  end if

! Which SRA components to use
  call find_keyinfile('isra',nchars,nlines,inputfile,iline,ipos,found)
  if (.not.found) then
    stop 'set_global_options: key isra not found!'
  else
    read(inputfile(iline)(ipos:nchars),*) isra
  end if

! Number of panels > 1
  call find_keyinfile('npanmax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) npanmax
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key npanmax not found!")')
    end if ! my_rank
!  Default value
    npanmax = 50
  endif

! These are optional, and can be set by the host KKR program:

! How many atoms to take into the projection scheme
  call find_keyinfile('nasusc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nasusc
  else 
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key nasusc not found!")')
    end if ! my_rank
    nasusc = natyp
  end if

! Angular momentum cutoff for GF
  call find_keyinfile('nlmax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nlmax
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key nlmax not found!")')
    end if ! my_rank
    nlmax = lmax
  end if

! Maximum number of radial points
  call find_keyinfile('nrmax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nrmax
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key nrmax not found!")')
    end if ! my_rank
    nrmax = nrad
  end if

! Number of spin components
  call find_keyinfile('nsmax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nsmax
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key nsmax not found!")')
    end if ! my_rank
    nsmax = nspin
  end if

! Maximum number of wfns in basis
  call find_keyinfile('nbmax',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nbmax
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_global_options: WARNING - key nbmax not found!")')
    end if ! my_rank
!   basis using default setup?
    lbasisdefault = .false.
    call find_keyinfile('lbasisdefault',nchars,nlines,inputfile,iline,ipos,found)
    if (found) read(inputfile(iline)(ipos:nchars),*) lbasisdefault
    if (lbasisdefault) then
!     default value
      nbmax = 8
    else
  !   search in input for maximum basis size
      ibasis = 2
      call find_keyinfile('ibasis',nchars,nlines,inputfile,iline,ipos,found)
      if (found) read(inputfile(iline)(ipos:nchars),*) ibasis
!     find how many atom groups there are
      call find_keyinfile('ngroup',nchars,nlines,inputfile,iline,ipos,found)
      if (.not.found) stop 'set_global_options: key ngroup not found!'
      read(inputfile(iline)(ipos:nchars),*) ngroup
!     now check basis set up
      do iline=1,nlines
        call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
        if (found) then
          read(inputfile(iline)(ipos:nchars),*) ig
          if (ig > ngroup) cycle
          call find_keyinline('iwsusc',nchars,nlines,inputfile,iline,ipos,found)
          if (found) then
            read(inputfile(iline)(ipos:nchars),*) iw(0:nlmax)
            if (ibasis == 1) nbmax = max(nbmax,maxval(iw(0:nlmax)))
            if (ibasis == 2) nbmax = max(nbmax,nsmax*maxval(iw(0:nlmax)))
            if (ibasis == 3) nbmax = max(nbmax,nsmax*sum(iw(0:nlmax)))
          end if
        end if
      end do
    end if
  end if

  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'(" nbmax=",i4," nlmax=",i2," sra=",i2," nrmax=",i6," nsmax=",i2)') nbmax, nlmax, isra, nrmax, nsmax
      write(*,'(" nesusc=",i6," nasusc=",i4)') nesusc, nasusc
    end if
  end if ! my_rank
  
! Some checks     
  if (lmax /= nlmax)  stop 'new_input: lmax /= nlmax'
  if (natyp < nasusc) stop 'new_input: natyp < nasusc'
  if (nrad > nrmax)   stop 'new_input: nrad > nrmax'
  if (nsmax /= nspin) stop 'new_input: nspin /= nsmax'
! All done!
  end subroutine set_global_options
