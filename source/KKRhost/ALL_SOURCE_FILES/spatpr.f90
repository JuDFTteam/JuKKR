!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Computes the triple product between three vectors to calculate a volume
!> Author:
!> Computes the triple product between three vectors to calculate a volume
!------------------------------------------------------------------------------------
module mod_spatpr
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Computes the triple product between three vectors to calculate a volume
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Computes the triple product between three vectors to calculate a volume
  !-------------------------------------------------------------------------------
  subroutine spatpr(a, b, c, v)

    implicit none

    real (kind=dp), dimension(*), intent (in) :: a !! Input vector
    real (kind=dp), dimension(*), intent (in) :: b !! Input vector
    real (kind=dp), dimension(*), intent (in) :: c !! Input vector
    real (kind=dp), intent (out) :: v !! volumne

    v = 0.0e0_dp
    v = v + c(1)*(a(2)*b(3)-a(3)*b(2))
    v = v + c(2)*(a(3)*b(1)-a(1)*b(3))
    v = v + c(3)*(a(1)*b(2)-a(2)*b(1))
    return
  end subroutine spatpr

end module mod_spatpr
