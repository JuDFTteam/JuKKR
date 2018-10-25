!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Simpson integration for a real integrand
!> Author: 
!> Simpson integration for a real integrand, \(f_x\) from 1 to `jtop` in an
!> equidistant mesh.
!> \begin{equation}
!> int =\left[F_1 + 4F_2 + 2F_3 + .... + 4F_{N-1} + F_N \right]
!> \end{equation}
!> for odd \(N\).
!------------------------------------------------------------------------------------
module mod_rintsimp
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Simpson integration for a real integrand
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> Simpson integration for a real integrand, \(f_x\) from 1 to `jtop` in an
  !> equidistant mesh.
  !> \begin{equation}
  !> int =\left[F_1 + 4F_2 + 2F_3 + .... + 4F_{N-1} + F_N \right]
  !> \end{equation}
  !> for odd \(N\).
  !-------------------------------------------------------------------------------
  subroutine rintsimp(fx, jtop, cint)

    implicit none

    ! Dummy arguments
    integer, intent(in) :: jtop !! Number of entries in the equidistant grid defining the integrand
    real (kind=dp), dimension(jtop), intent(in) :: fx !! Integrand
    real (kind=dp), intent(out) :: cint !! Integrated function

    ! Local variables
    integer :: i
    real (kind=dp) :: simp

    cint = fx(1)
    simp = -1.0e0_dp

    do i = 2, jtop - 1
      simp = -simp
      cint = cint + (3.0e0_dp+simp)*fx(i)
    end do

    cint = (cint+fx(jtop))/3.0e0_dp

  end subroutine rintsimp

end module mod_rintsimp
