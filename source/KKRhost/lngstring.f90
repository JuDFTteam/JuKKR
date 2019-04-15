!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Finds the length of a string
!> Author:
!> Finds the length of a string
!------------------------------------------------------------------------------------
!> @note Jonathan Chico: Wouldn't using the fortran intrinsic functions such as `trim()`
!> and `len()` be a better solution to this?
!>
!> Jonathan Chico: This routine seems to do the same than the `length()` function. Both
!> seem to be unnecessary if one instead uses fortran intrinsic functions.
!> @endnote
!------------------------------------------------------------------------------------
module mod_lngstring

  contains
  
    !-------------------------------------------------------------------------------
    !> Summary: Finds the length of a string
    !> Author: Who wrote this function
    !> Category: input-output, KKRhost 
    !> Deprecated: False 
    !> Finds the length of a string
    !-------------------------------------------------------------------------------
    !> @note Jonathan Chico: This routine seems to do the same than the `length()` function. Both
    !> seem to be unnecessary if one instead uses fortran intrinsic functions.
    !> @endnote
    !-------------------------------------------------------------------------------
    integer function lngstring(string, lstrmax)

    implicit none
  
      ! Dummy arguments
      integer, intent(in) :: lstrmax !! Maximum length of the string
      character (len=*), intent(in) :: string !! Input string
  
      ! Local variables
      character :: c
      integer :: i
  
      lngstring = 0
      do i = lstrmax, 1, -1
        c = string(i:i)
        if (c/=' ' .and. ichar(c)>0) then
          lngstring = i
          return
        end if
      end do
    end function lngstring
  
  end module mod_lngstring
