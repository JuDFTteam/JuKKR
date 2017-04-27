  subroutine set_soc_options(itc,my_rank)
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
  if (.not.lsoc) then
    if (my_rank == 0) then 
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No SOC correction")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank
    return
  end if

! Defaults
  isoc(1:nasusc) = 1; socscaling(1:nasusc) = 1.d0; ntot = nasusc

! Search for atom info lines
! first pass: which atoms for SOC
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     which atoms are in group ig
      call find_keyinline('isoc',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) i
        if (i /= 0 .and. i /= 1 .and. i /= 2 .and. i /= 3) stop 'set_soc_options: unknown isoc option'
        isoc(istart:iend) = i
        if (i == 0) ntot = ntot - igroup(ig)
      end if
    end if
  end do
! second pass: atomic SOC options
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
!     which atoms are in group ig
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     scaling of SOC potential
      call find_keyinline('socscaling',nchars,nlines,inputfile,iline,ipos,found)
      if (found .and. isoc(istart) /= 0) then
        read(inputfile(iline)(ipos:nchars),*) re
        socscaling(istart:iend) = re
      end if
    end if
  end do

! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" SOC correction  ==>  check atomic options")')
!     atom info
      write(*,'(/" Number of atoms for SOC calculation:",i4)') ntot
      write(*,'(/" SOC configuration for each atom group:")')
      do ig=1,ngroup
        istart = sum(igroup(1:ig-1)) + 1
        iend   = sum(igroup(1:ig))
        if (isoc(istart) /= 0) then
          write(*,'(" ig, istart, iend=",3i6,"  iasusc=",1000i4)') ig, istart, iend, iasusc(istart:iend)
          if (isoc(istart) == 1) write(*,'(" + SOC full: L.S  =>  scaling=",f8.4)') socscaling(istart)
          if (isoc(istart) == 2) write(*,'(" + SOC longitudinal: (L.u).(S.u)  =>  scaling=",f8.4)') socscaling(istart)
          if (isoc(istart) == 3) write(*,'(" + SOC transverse: L.S - (L.u).(S.u)  =>  scaling=",f8.4)') socscaling(istart)
        end if
      end do
      write(*,'("************************************************************"/)')
    end if 
  end if ! my_rank
! All done!
  end subroutine set_soc_options
