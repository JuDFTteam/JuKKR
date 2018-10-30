!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Adds two vectors together
!> Author: 
!> Adds two vectors together
!------------------------------------------------------------------------------------
module mod_vadd
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Adds two vectors together
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Adds two vectors together
  !-------------------------------------------------------------------------------
  subroutine vadd(a, b, c)

    real (kind=dp), dimension(*), intent (in) :: a !! Input vector
    real (kind=dp), dimension(*), intent (in) :: b !! Input vector
    real (kind=dp), dimension(*), intent (out) :: c !! Output vector

    integer :: i

    do i = 1, 3
      c(i) = a(i) + b(i)
    end do
    return
  end subroutine vadd

end module mod_vadd
