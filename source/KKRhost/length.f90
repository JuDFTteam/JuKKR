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
!> Jonathan Chico: This seems to be used only in the subroutine `sname()` which seems to be 
!> deprecated, maybe this funtcion is unneccessary.
!> @endnote
!------------------------------------------------------------------------------------
module mod_length

  contains
    !-------------------------------------------------------------------------------
    !> Summary: Finds the length of a string
    !> Author: Who wrote this function
    !> Category: input-output, KKrhost 
    !> Deprecated: False 
    !> Finds the length of a string
    !-------------------------------------------------------------------------------
    !> @note Jonathan Chico: This seems to be used only in the subroutine `sname()` which seems to be 
    !> deprecated, maybe this funtcion is unneccessary.
    !> @endnote
    !-------------------------------------------------------------------------------
    integer function length(s, max)
      ! ************************************************************************
  
      character (len=1), intent (inout) :: s(*) !! String array
      integer, intent (in) :: max !! Maximum size of the string
  
      integer :: i
      ! ------------------------------------------------------------------------
      i = max
  
      do while (s(i)==' ')
        i = i - 1
      end do
  
      length = i
  
      return
    end function length
  
  end module mod_length
