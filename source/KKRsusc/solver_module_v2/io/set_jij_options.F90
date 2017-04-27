  subroutine set_jij_options(itc,my_rank)
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
  real(kind=r8b)    :: re


! Not a SOC calculation
  if (.not.ljij) then 
    if (my_rank == 0) then
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No magnetic couplings")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank 
    return
  end if

! Jij options

! full tensor or scalar 
  call find_keyinfile('ljijtensor',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ljijtensor
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_soc_options: WARNING - key ljijtensor not found!")')
    end if ! my_rank
    ljijtensor = .true.
  end if

! correction from EF due to constant Ne 
  call find_keyinfile('ljijef',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ljijef
  else
    if (my_rank == 0) then 
      if (lwarn .and. itc == 1) write(*,'("set_soc_options: WARNING - key ljijef not found!")')
    end if ! my_rank
    ljijef = .false.
  end if

! onsite correction
  call find_keyinfile('ljionsite',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) ljionsite
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_soc_options: WARNING - key ljionsite not found!")')
    end if ! my_rank
    ljionsite = .false.
  end if

!
! Defaults
  ijij(1:nasusc) = 1; ntot = nasusc

! Search for atom info lines
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     which atoms are in group ig
      call find_keyinline('ijij',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) i
        if (i /= 0) then
          if (i /= 1) stop 'set_jij_options: unknown ijij option'
          ijij(istart:iend) = i
        else
          ntot = ntot - igroup(ig)
          ijij(istart:iend) = 0
        end if
      end if
    end if
  end do

! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" Magnetic couplings Jij  ==>  check atomic options")')
!     form of the tensor
      if (ljijtensor) then
        write(*,'(" + Tensor Jij using Ebert formula")')
        if (ljijef) write(*,'(" + EF correction to Jij included")')
      else
        write(*,'(" + Scalar Jij using Lichtenstein formula")')
      end if
      if (ljionsite) write(*,'("+ onsite correction to Ji included")')
!     atom info
      write(*,'(/" Number of atoms for Jij calculation:",i4)') ntot
      write(*,'(/" Jij configuration for each atom group:")')
      do ig=1,ngroup
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
        if (ijij(istart) /= 0) then
          write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
        end if
      end do
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank 
! All done!
  end subroutine set_jij_options
