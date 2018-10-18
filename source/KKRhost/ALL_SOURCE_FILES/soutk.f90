!------------------------------------------------------------------------------------
!> Summary: This subroutine does an outwards integration of a function with kinks
!> Author: B. Drittler
!> To integrate the function \(fint\left(r\right)=\int_{0}^{r}f\left(r'\right)dr' \)
!> at each kink the integration is restarted the starting value for this integration is determined by
!> a 4 point lagrangian integration, coefficients given by M. Abramowitz and I.A. Stegun,
!> handbook of mathematical functions, and applied mathematics series 55 (1968).
!> the weights \(drdi\) have to be multiplied before calling this subroutine.
!------------------------------------------------------------------------------------
module mod_soutk

contains

  !-------------------------------------------------------------------------------
  !> Summary: This subroutine does an outwards integration of a function with kinks
  !> Author: B. Drittler
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> To integrate the function \(fint\left(r\right)=\int_{0}^{r}f\left(r'\right)dr' \)
  !> at each kink the integration is restarted the starting value for this integration is determined by
  !> a 4 point lagrangian integration, coefficients given by M. Abramowitz and I.A. Stegun,
  !> handbook of mathematical functions, and applied mathematics series 55 (1968).
  !> the weights \(drdi\) have to be multiplied before calling this subroutine.
  !-------------------------------------------------------------------------------
  subroutine soutk(f, fint, ipan, ircut)

    use :: global_variables
    use :: mod_datatypes, only: dp

    ! .. Scalar Arguments
    integer, intent (in) :: ipan   !! Number of panels in non-MT-region
    ! .. Array Arguments
    integer, dimension (0:ipand), intent (in) :: ircut !! R points of panel borders
    real (kind=dp), dimension (*), intent (in) :: f
    ! .. Output variables
    real (kind=dp), dimension (*), intent (out) :: fint
    ! .. Local Scalars
    integer :: i, ien, ip, ist
    real (kind=dp) :: a1, a2
    ! ..
    a1 = 1.0e0_dp/3.0e0_dp
    a2 = 4.0e0_dp/3.0e0_dp

    ! ----------------------------------------------------------------------------
    ! Loop over kinks
    ! ----------------------------------------------------------------------------
    do ip = 1, ipan
      ien = ircut(ip)
      ist = ircut(ip-1) + 1

      if (ip==1) then
        fint(ist) = 0.0e0_dp
        ! ----------------------------------------------------------------------
        ! Integrate fint(ist+1) with a 4 point lagrangian
        ! ----------------------------------------------------------------------
        fint(ist+1) = (f(ist+3)-5.0e0_dp*f(ist+2)+19.0e0_dp*f(ist+1)+9.0e0_dp*f(ist))/24.0e0_dp
      else
        fint(ist) = fint(ist-1)
        ! ----------------------------------------------------------------------
        ! Integrate fint(ist+1) with a 4 point lagrangian
        ! ----------------------------------------------------------------------
        fint(ist+1) = fint(ist-1) + (f(ist+3)-5.0e0_dp*f(ist+2)+19.0e0_dp*f(ist+1)+9.0e0_dp*f(ist))/24.0e0_dp
      end if

      ! -------------------------------------------------------------------------
      ! Calculate fint with an extended 3-point-simpson
      ! -------------------------------------------------------------------------
      do i = ist + 2, ien
        fint(i) = ((fint(i-2)+f(i-2)*a1)+f(i-1)*a2) + f(i)*a1
      end do                       ! I
    end do                         ! IP

  end subroutine soutk

end module mod_soutk
