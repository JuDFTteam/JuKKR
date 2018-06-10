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
subroutine simpk(f, fint, ipan, ircut, drdi)

  use :: global_variables

  integer, intent (in) :: ipan !< Number of panels in non-MT-region
  integer, dimension (0:ipand), intent (in) :: ircut !< R points of panel borders
  double precision, dimension (*), intent (in) :: drdi !< Derivative dr/di
! .. Output variables
  double precision, intent (out) :: fint
! .. In/Out variables
  double precision, dimension (*), intent (inout) :: f
! .. Local Scalars
  integer :: i, ien, ip, ist, n
  double precision :: a1, a2

! .. External Functions
  double precision :: ssum
  external :: ssum
!     ..
  a1 = 4.0d0/3.0d0
  a2 = 2.0d0/3.0d0
  fint = 0.0d0
!
  do ip = 1, ipan
!-------------------------------------------------------------------------
! Loop over kinks
!-------------------------------------------------------------------------
    ist = ircut(ip-1) + 1
    ien = ircut(ip)
!
    do i = ist, ien
      f(i) = f(i)*drdi(i)
    end do ! I
    if (mod(ien-ist,2)==0) then
      fint = fint + (f(ist)-f(ien))/3.0d0
      ist = ist + 1
      n = (ien-ist+1)/2
    else
!----------------------------------------------------------------------
! Four point Lagrange integration for the first step
!----------------------------------------------------------------------
      fint = fint + (9.0d0*f(ist)+19.0d0*f(ist+1)-5.0d0*f(ist+2)+f(ist+3))/ &
        24.0d0 + (f(ist+1)-f(ien))/3.0d0
      ist = ist + 2
      n = (ien-ist+1)/2
    end if
!-------------------------------------------------------------------------
! Calculate with an extended 3-point-simpson
!-------------------------------------------------------------------------
    fint = fint + a1*ssum(n, f(ist), 2) + a2*ssum(n, f(ist+1), 2)
  end do ! IP
end subroutine
