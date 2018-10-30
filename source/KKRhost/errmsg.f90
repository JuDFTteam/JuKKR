!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_errmsg

contains

  !-------------------------------------------------------------------------------
  !> Summary: Writes error message and strops code if severity level is too high
  !> Author: 
  !> Category: KKRhost, sanity-check
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Write error message to message error device
  !> 
  !> Inputs:
  !>   messg : error message
  !>   isev  : severity level
  !> Remarks
  !>   if severity level greater or equal than error tolerance
  !>   program will stop.
  !-------------------------------------------------------------------------------
  subroutine errmsg(messg, isev)

    use :: mod_types, only: t_inc
    implicit none

    integer :: isev
    character (len=*) :: messg
    ! Local parameters:
    integer :: iline, iisev, ipos(0:20), l, nline, nunit
    character (len=14) :: c(1:4)

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

  end subroutine errmsg

end module mod_errmsg
