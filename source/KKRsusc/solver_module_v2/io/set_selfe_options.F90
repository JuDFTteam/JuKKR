  subroutine set_selfe_options(itc,my_rank)
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
  integer(kind=i4b) :: ig, istart, iend
  real(kind=r8b)    :: ere, eim


! No self-energy calculation
  if (.not.lsemodel) then
    if (my_rank == 0) then
      if (itc == 1) then
        write(*,'("************************************************************")')
        write(*,'(" No self-energy read in")')
        write(*,'("************************************************************"/)')
      end if
    end if ! my_rank
    return
  end if

! ----------------------------------------------------------------------
! Print options summary
  if (my_rank == 0) then
    if (itc == 1) then
      write(*,'("************************************************************")')
      write(*,'(" Model self-energy read in  ==>  provide semodel.dat")')
      write(*,'("************************************************************"/)')
    end if
  end if ! my_rank 
! All done!
  end subroutine set_selfe_options
