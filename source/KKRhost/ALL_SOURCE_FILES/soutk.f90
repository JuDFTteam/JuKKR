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
    Subroutine soutk(f, fint, ipan, ircut)

      Use global_variables
      Use mod_datatypes, Only: dp

! .. Scalar Arguments
      Integer, Intent (In) :: ipan !< Number of panels in non-MT-region
! .. Array Arguments
      Integer, Dimension (0:ipand), Intent (In) :: ircut !< R points of panel borders
      Real (Kind=dp), Dimension (*), Intent (In) :: f
! .. Output variables
      Real (Kind=dp), Dimension (*), Intent (Out) :: fint
! .. Local Scalars
      Integer :: i, ien, ip, ist
      Real (Kind=dp) :: a1, a2
! ..
      a1 = 1.0E0_dp/3.0E0_dp
      a2 = 4.0E0_dp/3.0E0_dp

!----------------------------------------------------------------------------
! Loop over kinks
!----------------------------------------------------------------------------
      Do ip = 1, ipan
        ien = ircut(ip)
        ist = ircut(ip-1) + 1
!
        If (ip==1) Then
          fint(ist) = 0.0E0_dp
!----------------------------------------------------------------------
! Integrate fint(ist+1) with a 4 point lagrangian
!----------------------------------------------------------------------
          fint(ist+1) = (f(ist+3)-5.0E0_dp*f(ist+2)+19.0E0_dp*f(ist+1)+ &
            9.0E0_dp*f(ist))/24.0E0_dp
        Else
          fint(ist) = fint(ist-1)
!----------------------------------------------------------------------
! Integrate fint(ist+1) with a 4 point lagrangian
!----------------------------------------------------------------------
          fint(ist+1) = fint(ist-1) + (f(ist+3)-5.0E0_dp*f(ist+2)+19.0E0_dp*f( &
            ist+1)+9.0E0_dp*f(ist))/24.0E0_dp
        End If

!-------------------------------------------------------------------------
! Calculate fint with an extended 3-point-simpson
!-------------------------------------------------------------------------
        Do i = ist + 2, ien
          fint(i) = ((fint(i-2)+f(i-2)*a1)+f(i-1)*a2) + f(i)*a1
        End Do ! I
      End Do ! IP
!
    End Subroutine
