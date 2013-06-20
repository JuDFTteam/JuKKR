!------------------------------------------------------------------------------
!> Storage for Euler angles.
module angles_common
  use shape_constants_mod
  implicit none

  save

  !real(kind=DP), dimension(NFACED)::ALPHA    !< Euler angles alpha to rotate faces perpendicular to z-axis
  !real(kind=DP), dimension(NFACED)::BETA     !< Euler angles beta
  !real(kind=DP), dimension(NFACED)::GAMMA    !< Euler angles gamma

  real(kind=DP), dimension(:), allocatable ::ALPHA    !< Euler angles alpha to rotate faces perpendicular to z-axis
  real(kind=DP), dimension(:), allocatable ::BETA     !< Euler angles beta
  real(kind=DP), dimension(:), allocatable ::GAMMA    !< Euler angles gamma


  contains

  subroutine createangles(nfaced)
    implicit none
    integer, intent(in) :: nfaced

    allocate(ALPHA(nfaced))
    allocate(BETA(nfaced))
    allocate(GAMMA(nfaced))
  end subroutine

  subroutine destroyangles()
    implicit none

    deallocate(ALPHA)
    deallocate(BETA)
    deallocate(GAMMA)
  end subroutine

end module
