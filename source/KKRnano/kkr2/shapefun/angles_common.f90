!------------------------------------------------------------------------------
!> Storage for Euler angles.
module angles_common
  use shape_constants_mod, only: DP
  implicit none
  public

  real(kind=DP), allocatable :: ALPHA(:)    !< Euler angles alpha to rotate faces perpendicular to z-axis
  real(kind=DP), allocatable :: BETA(:)     !< Euler angles beta
  real(kind=DP), allocatable :: GAMMA(:)    !< Euler angles gamma

  contains

  subroutine createangles(nfaced)
    integer, intent(in) :: nfaced
    allocate(ALPHA(nfaced))
    allocate(BETA(nfaced))
    allocate(GAMMA(nfaced))
  end subroutine

  subroutine destroyangles()
    deallocate(ALPHA)
    deallocate(BETA)
    deallocate(GAMMA)
  end subroutine

end module
