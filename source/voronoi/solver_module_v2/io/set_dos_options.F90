  subroutine set_dos_options(itc,my_rank)
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
  integer(kind=i4b) :: ig, istart, iend, ntot, i, il,is, i4by2(0:nlmax,nsmax)
  real(kind=r8b)    :: ere, eim


! No DOS calculation
  if (.not.ldos) then
    if (my_rank == 0) then
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No DOS calculated")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank
    return
  end if

! DOS options

! use projection energies or straigh line
  call find_keyinfile('idos',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) idos
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_dos_options: WARNING - key idos not found!")')
    end if ! my_rank
    idos = 1
  end if
! straight line option only available with GF fit
  if (.not.lfit) idos = 1
  if (idos /=1 .and. idos /= 2)  stop 'set_dos_options: unknown DOS option'

! number of energies for straight line
  call find_keyinfile('nedos',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nedos
  else
    if (my_rank == 0) then
      if (lwarn .and. idos == 2) write(*,'("set_dos_options: WARNING - key nedos not found!")')
    end if ! my_rank
    nedos = 201
  end if

! starting point for straight line
  call find_keyinfile('dose0',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ere, eim
    dose0 = cmplx(ere,eim)
  else
    if (my_rank == 0) then
      if (lwarn .and. idos == 2) write(*,'("set_dos_options: WARNING - key dose0 not found!")')
    end if ! my_rank
    dose0 = (-0.3d0,0d0)
  end if

! end point for straight line
  call find_keyinfile('dose1',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ere, eim
    dose1 = cmplx(ere,eim)
  else
    if (my_rank == 0) then
      if (lwarn .and. idos == 2) write(*,'("set_dos_options: WARNING - key dose1 not found!")')
    end if ! my_rank
    dose1 = (1.0d0,0d0)
  end if

! analytical continuation to real axis
  call find_keyinfile('ldosacon',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ldosacon
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_dos_options: WARNING - key ldosacon not found!")')
    end if ! my_rank
    ldosacon = .false.
  end if

! number of energies for rational function extrapolation
  call find_keyinfile('nedosacon',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) nedosacon
  else 
    if (my_rank == 0) then 
      if (lwarn) write(*,'("set_dos_options: WARNING - key nedosacon not found!")')
    end if ! my_rank
    nedosacon = 4
  end if

! energy zero
  call find_keyinfile('dosezero',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) dosezero
  else
    if (my_rank == 0) then
      if (lwarn) write(*,'("set_dos_options: WARNING - key dosezero not found!")')
    end if ! my_rank
    dosezero = 0.d0
  end if

! energy conversion factor
  call find_keyinfile('dosefac',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) dosefac
  else
    if (my_rank == 0) then  
      if (lwarn) write(*,'("set_dos_options: WARNING - key dosefac not found!")') 
    end if ! my_rank
    dosefac = 1.d0
  end if

! density matrix from (lm,l'm') blocks of GF
  call find_keyinfile('ldosdmat',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ldosdmat
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_dos_options: WARNING - key ldosdmat not found!")')
    end if ! my_rank
    ldosdmat = .false.
  end if

! Defaults
  iadmat(1:nasusc) = 0; ildmat(0:nlmax,1:nsmax,1:nasusc) = 1
  ntot = 0

! ******************
  if (ldosdmat) then
! ******************
! Search for atom info lines
! first pass: which atoms included in density matrix calculation
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
!     which atoms are in group ig
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     density matrix for this atom group
      call find_keyinline('iadmat',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) i
        if (i == 0) cycle
        if (i /= 1) stop 'set_ldos_options: unknown iadmat option'
        iadmat(istart:iend) = i
        ntot = ntot + igroup(ig)
      end if
      call find_keyinline('ildmat',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) i4by2(0:nlmax,1:nsmax)
        do is=1,nsmax
          do il=0,nlmax
            if (i4by2(il,is) /= 0 .and. i4by2(il,is) /= 1) stop 'set_ldos_options: unknown ildmat option'
          end do
        end do
        do i=istart,iend
          ildmat(0:nlmax,1:nsmax,i) = i4by2(0:nlmax,1:nsmax)
        end do
      end if
    end if
  end do
! ******
  end if
! ******


! ----------------------------------------------------------------------
! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      if (idos == 1) write(*,'(" DOS calculation for the projection energies")')
      if (idos == 2)  then
        write(*,'(" DOS calculation for straight line:")')
        write(*,'(" + nedos, ebottom, etop=",i6,4f10.6)') nedos, dose0, dose1
      end if
      if (ldosacon) then
        write(*,'(" Analytical continuation to real axis using rational function")')
        write(*,'(" + nedosacon=",i6)') nedosacon
      end if
      write(*,'(" Energy zero=",f16.8,"  conversion factor=",f16.8)') dosezero, dosefac
      if (ldosdmat) then
        write(*,'("************************************************************"/)')
!       Atom info
        write(*,'(/" Number of atoms for density matrix calculation:",i4)') ntot
        write(*,'(/" density matrix configuration for each atom group:")')
        do ig=1,ngroup
          istart = sum(igroup(1:ig-1)) + 1
          iend   = sum(igroup(1:ig))
          if (iadmat(istart) /= 0) then
            write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
            write(*,'(" + density matrix for il,is=",10i4)') ildmat(0:nlmax,1:nsmax,istart)
          end if
        end do
      end if
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank
! All done!
  end subroutine set_dos_options
