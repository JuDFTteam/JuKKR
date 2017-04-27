  subroutine find_keyinline(key,nchars,nlines,inputfile,iline,ipos,found)
! scans chosen line of input file looking for key
! it looks for ' ' // keyname // '='
! this means keyname is separated on the left by space and on the right by =

  implicit none

! target key
  character(len=*),      intent(in)  :: key
! input file
  integer,               intent(in)  :: nchars, nlines
  character(len=nchars), intent(in)  :: inputfile(nlines)
! line where key is located
  integer,               intent(in)  :: iline
! position of last character in key
  integer,               intent(out) :: ipos
! was the key found?
  logical,               intent(out) :: found
! ----------------------------------------------------------------------
  integer :: keylen, iend

  found = .false.; ipos = -1
! size of key
  keylen = len(key) + 2
! look for key
! check for comments
  iend = index(inputfile(iline),'#')
!  write(*,'("line=",i8,"  comment mark at i=",i8)') iline, iend
! this is a comment line
  if (iend == 1) return
! no comments on this line
  if (iend == 0) iend = nchars
! check for key in current line, excluding the comments
  ipos = index(inputfile(iline)(1:iend-1),' '//key//'=')
! key found
  if (ipos > 0) then
!   this is the first position after the key in iline
    ipos  = ipos + keylen
    found = .true.
  end if
! All done!
  end subroutine find_keyinline
