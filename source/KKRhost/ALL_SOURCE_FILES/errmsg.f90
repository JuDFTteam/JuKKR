subroutine errmsg(messg, isev)
!- Write error message to message error device
! ----------------------------------------------------------------------
!i Inputs:
!i   messg :error message
!i   isev  :severity level
!r Remarks
!r   if severity level greater or equal than error tolerance
!r   program will stop.
! ----------------------------------------------------------------------
  use :: mod_types, only: t_inc
      Use mod_datatypes, Only: dp
  implicit none
! Passed parameters:                                                    
  integer :: isev
  character (len=*) :: messg
! Local parameters:                                                     
  integer :: iline, iisev, ipos(0:20), l, nline, nunit
  character (len=14) :: c(1:4)
! External calls:                                                       
! Intrinsic functions                                                   
  intrinsic :: iabs, max0, min0

  data c/'Information:', 'Warning:', 'Error:', 'Fatal error:'/

  iisev = max0(min0(isev,4), 1)

  ipos(0) = 1
  nline = 0
  l = 1
  do while (messg(l:l)/='$' .and. l<500)
    if (messg(l:l)=='|') then
      nline = nline + 1
      ipos(nline) = l
    end if
    l = l + 1
  end do
  nline = nline + 1
  ipos(nline) = l

  nunit = 1
  if ((nunit==1) .and. (t_inc%i_write>0)) write (1337, *)
  if (t_inc%i_write>0) then
    do iline = 1, nline
      write (1337, 100) c(iisev), messg(ipos(iline-1)+1:ipos(iline)-1)
    end do
  end if

  if (iabs(isev)>=3) stop

100 format (a13, a)

end subroutine
