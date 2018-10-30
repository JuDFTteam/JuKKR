  subroutine find_keyinfile(key,nchars,nlines,inputfile,iline,ipos,found)
! scans input file looking for key
! it looks for ' ' // keyname // '='
! this means keyname is separated on the left by space and on the right by =

  implicit none

! target key
  character(len=*),      intent(in)  :: key
! input file
  integer,               intent(in)  :: nchars, nlines
  character(len=nchars), intent(in)  :: inputfile(nlines)
! line where key is located
  integer,               intent(out) :: iline
! position of last character in key
  integer,               intent(out) :: ipos
! was the key found?
  logical,               intent(out) :: found
! ----------------------------------------------------------------------

! look for key in all lines in file
  do iline=1,nlines
    call find_keyinline(key,nchars,nlines,inputfile,iline,ipos,found)
    if (found) exit
  end do
! All done!
  end subroutine find_keyinfile
