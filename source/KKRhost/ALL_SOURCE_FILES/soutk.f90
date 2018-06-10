!-------------------------------------------------------------------------------
! SUBROUTINE: SOUTK
!> @brief This subroutine does an outwards integration of a function with kinks
!> @detail To integrate the function \f$ fint\left(r\right)=\int_{0}^{r}f\left(r'\right)dr' \f$
!> at each kink the integration is restarted the starting value for this integration is determined by
!> a 4 point lagrangian integration, coefficients given by M. Abramowitz and I.A. Stegun,
!> handbook of mathematical functions, nbs applied mathematics series 55 (1968).
!> the weights \f$drdi\f$ have to be multiplied before calling this subroutine.
!> @author B. Drittler
!> @date Oct. 1989
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine soutk(f, fint, ipan, ircut)

  use :: global_variables

! .. Scalar Arguments
  integer, intent (in) :: ipan !< Number of panels in non-MT-region
! .. Array Arguments
  integer, dimension (0:ipand), intent (in) :: ircut !< R points of panel borders
  double precision, dimension (*), intent (in) :: f
! .. Output variables
  double precision, dimension (*), intent (out) :: fint
! .. Local Scalars
  integer :: i, ien, ip, ist
  double precision :: a1, a2
! ..
  a1 = 1.0d0/3.0d0
  a2 = 4.0d0/3.0d0

!----------------------------------------------------------------------------
! Loop over kinks
!----------------------------------------------------------------------------
  do ip = 1, ipan
    ien = ircut(ip)
    ist = ircut(ip-1) + 1
!
    if (ip==1) then
      fint(ist) = 0.0d0
!----------------------------------------------------------------------
! Integrate fint(ist+1) with a 4 point lagrangian
!----------------------------------------------------------------------
      fint(ist+1) = (f(ist+3)-5.0d0*f(ist+2)+19.0d0*f(ist+1)+9.0d0*f(ist))/ &
        24.0d0
    else
      fint(ist) = fint(ist-1)
!----------------------------------------------------------------------
! Integrate fint(ist+1) with a 4 point lagrangian
!----------------------------------------------------------------------
      fint(ist+1) = fint(ist-1) + (f(ist+3)-5.0d0*f(ist+2)+19.0d0*f(ist+1)+ &
        9.0d0*f(ist))/24.0d0
    end if

!-------------------------------------------------------------------------
! Calculate fint with an extended 3-point-simpson
!-------------------------------------------------------------------------
    do i = ist + 2, ien
      fint(i) = ((fint(i-2)+f(i-2)*a1)+f(i-1)*a2) + f(i)*a1
    end do ! I
  end do ! IP
!
end subroutine
