!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Check if a string is contained in the `testc` array
!> Author: 
!> Category: undefined, KKRhost
!> Deprecated: False 
!> Check if a string is contained in the `testc` array
!-------------------------------------------------------------------------------
!> @note Jonathan Chico: This function is not contained in a module, it should be 
!> included in a tool like module where auxiliary functions and subrotuines are 
!> located.
!> @endnote
!-------------------------------------------------------------------------------
logical function test(string)
  ! ***********************************************************************

  ! TEST = 'STRING  ' IS CONTAINED IN /TESTC/.

  ! ------------------------------------------------------------------------
  use :: mod_wunfiles, only: t_params

  implicit none
  character (len=8), intent (in) :: string
  integer :: i
  character (len=8), dimension(32) :: testc

  testc = t_params%testc

  test = .false.
  do i = 1, 32
    if (string==testc(i)) test = .true.
  end do
  return
end function test