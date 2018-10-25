!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Multiplication of a vector by a scalar 
!> Author: 
!> Multiplication of a vector by a scalar 
!------------------------------------------------------------------------------------
!> @note Jonathan Chico: This seems unnecessary as one can just replace this by a
!> call `c(:)=a(:)*b` in the code
!> @endnote
!------------------------------------------------------------------------------------
module mod_vmul
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Multiplication of a vector by a scalar 
  !> Author: 
  !> Category: numerical-tools, KKRhost 
  !> Deprecated: False
  !> Multiplication of a vector by a scalar
  !-------------------------------------------------------------------------------
  !> @note Jonathan Chico: This seems unnecessary as one can just replace this by a
  !> call `c(:)=a(:)*b` in the code
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine vmul(a, b, c)

    real (kind=dp), intent (in) :: b  !! Input scalar
    real (kind=dp), dimension(*), intent (in) :: a !! Input vector
    real (kind=dp), dimension(*), intent (out) :: c !! Output vector

    integer :: i

    do i = 1, 3
      c(i) = b*a(i)
    end do
    return
  end subroutine vmul

end module mod_vmul
