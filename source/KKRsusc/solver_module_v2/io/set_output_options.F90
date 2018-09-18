  subroutine set_output_options(itc,my_rank)
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
  integer(kind=i4b) :: ig, istart, iend, i, j
  real(kind=r8b)    :: ere, eim

! Full output for the susceptibility calculation
  call find_keyinfile('loutsusc',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) loutsusc
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_output_options: WARNING - key loutsusc not found!")')
      loutsusc = .false.
    end if ! my_rank
  end if

! Full output for the preparation steps
  call find_keyinfile('loutbasis',nchars,nlines,inputfile,iline,ipos,found)
  if (found) then
    read(inputfile(iline)(ipos:nchars),*) loutbasis
  else
    if (my_rank == 0) then
      if (lwarn .and. itc == 1) write(*,'("set_output_options: WARNING - key loutbasis not found!")')
      loutbasis = .false.
    end if ! my_rank
  end if


! If warnings are turned on full output is used
  if(lwarn) then
    loutsusc = .true.
    loutbasis = .true.
  end if 

  end subroutine set_output_options

