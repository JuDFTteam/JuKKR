!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Sets two vectors to be equal
!> Author: 
!> Sets two vectors to be equal
!------------------------------------------------------------------------------------
!> @note Jonathan Chico: This seems unnecessary, it could be replace by a call like
!> `b(:)=a(:)`
!> @endnote
!------------------------------------------------------------------------------------
module mod_veq
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Sets two vectors to be equal
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False
  !> Sets two vectors to be equal
  !-------------------------------------------------------------------------------
  !> @note Jonathan Chico: This seems unnecessary, it could be replace by a call like
  !> `b(:)=a(:)`
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine veq(a, b)

    implicit none

    real (kind=dp), dimension(*), intent(in) :: a  !! Input vector
    real (kind=dp), dimension(*), intent(out) :: b !! Output vector

    integer :: i

    do i = 1, 3
      b(i) = a(i)
    end do
  end subroutine veq

end module mod_veq
