!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: This subroutine does an integration up to \(r_{cut}\) of an real function \(f\) with an extended 3-point-simpson 
!> Author: 
!> The integration of the function $$ fint =\int_{0}^{r_{cut}} f\left(r'\right)dr'$$ has been modified
!> for functions with kinks - at each kink the integration is restarted.
!------------------------------------------------------------------------------------
!> @warning Input \(f\) is destroyed
!> @endwarning
!------------------------------------------------------------------------------------
module mod_simpk

contains

  !-------------------------------------------------------------------------------
  !> Summary: This subroutine does an integration up to \(r_{cut}\) of an real function \(f\) with an extended 3-point-simpson 
  !> Author:
  !> Category: numerical-tools, KKRhost 
  !> Deprecated: False 
  !> The integration of the function $$ fint =\int_{0}^{r_{cut}} f\left(r'\right)dr'$$ has been modified
  !> for functions with kinks - at each kink the integration is restarted.
  !-------------------------------------------------------------------------------
  !> @warning Input \(f\) is destroyed 
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine simpk(f, fint, ipan, ircut, drdi)

    use :: global_variables, only: ipand
    use :: mod_datatypes, only: dp
    use :: mod_ssum, only: ssum

    integer, intent (in) :: ipan   ! < Number of panels in non-MT-region
    integer, dimension (0:ipand), intent (in) :: ircut ! < R points of panel borders
    real (kind=dp), dimension (*), intent (in) :: drdi ! < Derivative dr/di
    ! .. Output variables
    real (kind=dp), intent (out) :: fint
    ! .. In/Out variables
    real (kind=dp), dimension (*), intent (inout) :: f
    ! .. Local Scalars
    integer :: i, ien, ip, ist, n
    real (kind=dp) :: a1, a2
    ! ..
    a1 = 4.0e0_dp/3.0e0_dp
    a2 = 2.0e0_dp/3.0e0_dp
    fint = 0.0e0_dp

    do ip = 1, ipan
      ! -------------------------------------------------------------------------
      ! Loop over kinks
      ! -------------------------------------------------------------------------
      ist = ircut(ip-1) + 1
      ien = ircut(ip)

      do i = ist, ien
        f(i) = f(i)*drdi(i)
      end do                       ! I
      if (mod(ien-ist,2)==0) then
        fint = fint + (f(ist)-f(ien))/3.0e0_dp
        ist = ist + 1
        n = (ien-ist+1)/2
      else
        ! ----------------------------------------------------------------------
        ! Four point Lagrange integration for the first step
        ! ----------------------------------------------------------------------
        fint = fint + (9.0e0_dp*f(ist)+19.0e0_dp*f(ist+1)-5.0e0_dp*f(ist+2)+f(ist+3))/24.0e0_dp + (f(ist+1)-f(ien))/3.0e0_dp
        ist = ist + 2
        n = (ien-ist+1)/2
      end if
      ! -------------------------------------------------------------------------
      ! Calculate with an extended 3-point-simpson
      ! -------------------------------------------------------------------------
      fint = fint + a1*ssum(n, f(ist), 2) + a2*ssum(n, f(ist+1), 2)
    end do                         ! IP
  end subroutine simpk

end module mod_simpk
