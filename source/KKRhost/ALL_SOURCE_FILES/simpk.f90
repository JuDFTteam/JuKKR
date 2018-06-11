!-------------------------------------------------------------------------------
! SUBROUTINE: SINWK
!> @brief This subroutine does an integration up to \f$ r_{cut}\f$ of an real function
!> \f$f\f$ with an extended 3-point-simpson
!> @details The integration of the function
!> \f$ fint =\int_{0}^{r_{cut}} f\left(r'\right)dr'\f$ has been modified
!> for functions with kinks - at each kink the integration is restarted.
!> @note Attention : Input \f$f\f$ is destroyed !
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine simpk(f, fint, ipan, ircut, drdi)

      Use global_variables
      Use mod_datatypes, Only: dp

      Integer, Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (0:ipand), Intent (In) :: ircut !< R points of panel borders
      Real (Kind=dp), Dimension (*), Intent (In) :: drdi !< Derivative dr/di
! .. Output variables
      Real (Kind=dp), Intent (Out) :: fint
! .. In/Out variables
      Real (Kind=dp), Dimension (*), Intent (Inout) :: f
! .. Local Scalars
      Integer :: i, ien, ip, ist, n
      Real (Kind=dp) :: a1, a2

! .. External Functions
      Real (Kind=dp) :: ssum
      External :: ssum
!     ..
      a1 = 4.0E0_dp/3.0E0_dp
      a2 = 2.0E0_dp/3.0E0_dp
      fint = 0.0E0_dp
!
      Do ip = 1, ipan
!-------------------------------------------------------------------------
! Loop over kinks
!-------------------------------------------------------------------------
        ist = ircut(ip-1) + 1
        ien = ircut(ip)
!
        Do i = ist, ien
          f(i) = f(i)*drdi(i)
        End Do ! I
        If (mod(ien-ist,2)==0) Then
          fint = fint + (f(ist)-f(ien))/3.0E0_dp
          ist = ist + 1
          n = (ien-ist+1)/2
        Else
!----------------------------------------------------------------------
! Four point Lagrange integration for the first step
!----------------------------------------------------------------------
          fint = fint + (9.0E0_dp*f(ist)+19.0E0_dp*f(ist+1)-5.0E0_dp*f(ist+2)+ &
            f(ist+3))/24.0E0_dp + (f(ist+1)-f(ien))/3.0E0_dp
          ist = ist + 2
          n = (ien-ist+1)/2
        End If
!-------------------------------------------------------------------------
! Calculate with an extended 3-point-simpson
!-------------------------------------------------------------------------
        fint = fint + a1*ssum(n, f(ist), 2) + a2*ssum(n, f(ist+1), 2)
      End Do ! IP
    End Subroutine
