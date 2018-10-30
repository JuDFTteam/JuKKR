!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: This subroutine does an integration from `istart` to `iend` of the real function `f` with an extended 3-point-simpson
!> Author: 
!> This subroutine does an integration from `istart` to `iend` of the real 
!> function `f` with an extended 3-point-simpson
!> \begin{equation}
!> f_{int}=\int_{ined}^{istart} f\left(r'\right)dr'
!> \end{equation}
!------------------------------------------------------------------------------------
module mod_simp3
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: This subroutine does an integration from `istart` to `iend` of the real function `f` with an extended 3-point-simpson
  !> Author: 
  !> Category: numerical-tools, KKRhost 
  !> Deprecated: False 
  !> This subroutine does an integration from `istart` to `iend` of the real 
  !> function `f` with an extended 3-point-simpson
  !> \begin{equation}
  !> f_{int}=\int_{ined}^{istart} f\left(r'\right)dr'
  !> \end{equation}
  !-------------------------------------------------------------------------------
  subroutine simp3(f,fint,istart,iend,drdi)
    ! .. Input variables
    integer, intent(in) :: iend
    integer, intent(in) :: istart
    real (kind=dp), dimension(*), intent(in) :: f !! Function to integrate
    real (kind=dp), dimension(*), intent(in) :: drdi !! Derivative dr/di
    ! .. Output variables
    real (kind=dp), intent(out) :: fint !! Integrated function
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: a1, a2
    integer :: i, ist
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: mod
    ! ..
    a1 = 4.0e0_dp/3.0e0_dp
    a2 = 2.0e0_dp/3.0e0_dp

    ! ---> initialize fint

    if (mod(iend-istart,2)==0) then
      fint = f(istart)*drdi(istart)/3.0e0_dp
      ist = istart + 1

    else
      fint = (f(istart+3)*drdi(istart+3)-5.0e0_dp*f(istart+2)*drdi(istart+2)+19.0e0_dp*f(istart+1)*drdi(istart+1)+9.0e0_dp*f(istart)*drdi(istart))/24.0e0_dp + &
        f(istart+1)*drdi(istart+1)/3.0e0_dp
      ist = istart + 2
    end if

    ! ---> calculate with an extended 3-point-simpson

    do i = ist, iend - 1, 2
      fint = fint + a1*f(i)*drdi(i) + a2*f(i+1)*drdi(i+1)
    end do
    fint = fint - f(iend)*drdi(iend)/3.0e0_dp

  end subroutine simp3

end module mod_simp3
