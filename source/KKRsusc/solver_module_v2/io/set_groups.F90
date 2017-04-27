  subroutine set_groups(itc,my_rank)
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
  integer(kind=i4b) :: ig, istart, iend, ntot


! Atom groups and iasusc
  call find_keyinfile('igroup',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    call find_keyinfile('ngroup',nchars,nlines,inputfile,iline,ipos,found)
    if (.not.found) stop 'set_groups: key ngroup not found!'
    read(inputfile(iline)(ipos:nchars),*) ngroup
    call find_keyinfile('igroup',nchars,nlines,inputfile,iline,ipos,found)
    read(inputfile(iline)(ipos:nchars),*) igroup(1:ngroup)
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_groups: WARNING - key igroup not found!")')
    end if ! my_rank
    ngroup = nasusc
    igroup(1:nasusc) = 1
    do ig=1,nasusc
      iasusc(ig) = ig
    end do
  end if
! check if number of atoms matches with groups
  ntot = sum(igroup(1:ngroup))
  if (my_rank == 0) then
    write(*,'(" ngroup=",i6," ntot=",i6)') ngroup, ntot
    write(*,'(" igroup=",20i4)') igroup(1:ngroup)
    write(*,*)  ! new line
  end if ! my_rank
  if (ntot /= nasusc) stop 'set_groups: atom groups and nasusc'
! igroup not found in input, exit here
  if (.not.found) return

! search for atom info lines
  iasusc(1:nasusc) = -1000
  do iline=1,nlines
    call find_keyinline('ig',nchars,nlines,inputfile,iline,ipos,found)
    if (found) then
      read(inputfile(iline)(ipos:nchars),*) ig
      if (ig > ngroup) cycle
      istart = sum(igroup(1:ig-1)) + 1
      iend   = sum(igroup(1:ig))
!     which atoms are in group ig
      call find_keyinline('ia',nchars,nlines,inputfile,iline,ipos,found)
      if (found) then
        read(inputfile(iline)(ipos:nchars),*) iasusc(istart:iend)
      end if
    end if
  end do
  if (any(iasusc == -1000)) stop 'set_groups: iasusc not set properly with igroup'
  do ig=1,ngroup
    istart = sum(igroup(1:ig-1)) + 1
    iend   = sum(igroup(1:ig))
    if (my_rank == 0 ) then
      write(*,'(" iasusc=",20i4)') iasusc(istart:iend)
    end if ! my_rank
  end do
! All done!
  end subroutine set_groups
