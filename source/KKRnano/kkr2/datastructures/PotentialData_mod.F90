!> Note: VONS is wastefully dimensioned... it is allocated for total number of
!> mesh points, but only non-spherical part is non-zero
!> also there is always space for 2 spin-directions.

! eliminates: VINS, VISP, VONS, lmpotd
! TODO: pointer to radial mesh? unnecessary - just pass atom
! ADD lmax???
! parameters are more or less the same for all atoms - given in inputcard
! write lpot, nspin, irmind, irmd for every atom??
module PotentialData_mod
  implicit none

  type PotentialData
    double precision, dimension(:,:,:), allocatable :: VINS        ! .. input potential
    double precision, dimension(:,:), allocatable :: VISP
    double precision, dimension(:,:,:), allocatable :: VONS        !     .. output potential
    integer :: nspin
    integer :: lpot
    integer :: irmind
    integer :: irmd
    integer :: irnsd
    ! derived
    integer :: lmpot

    ! ECORE, NCORE, ITITLE, LCORE ??? - separate to Core-datastructure? TODO
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createPotentialData(potential, lpot, nspin, irmind, irmd)
    implicit none
    type (PotentialData), intent(inout) :: potential
    integer, intent(in) :: lpot
    integer, intent(in) :: nspin
    integer, intent(in) :: irmind  !< number of mesh points
    integer, intent(in) :: irmd  !< number of mesh points

    ! ---- local ----------
    integer :: lmpot

    potential%nspin = nspin
    potential%lpot = lpot
    potential%irmind = irmind
    potential%irmd = irmd
    potential%irnsd = irmd - irmind

    lmpot = (lpot + 1)**2
    potential%lmpot = lmpot

    ! Unfortunately, for compatibility reasons,
    ! memory has to be allocated for both spin directions
    ! independent of nspin
    allocate(potential%VINS(IRMIND:IRMD,LMPOT,2))
    allocate(potential%VISP(IRMD,2))
    allocate(potential%VONS(IRMD,LMPOT,2))
    potential%VINS = 0.0d0
    potential%VISP = 0.0d0
    potential%VONS = 0.0d0

  end subroutine


  !----------------------------------------------------------------------------
  subroutine destroyPotentialData(potential)
    implicit none
    type (PotentialData), intent(inout) :: potential

    deallocate(potential%VINS)
    deallocate(potential%VISP)
    deallocate(potential%VONS)
  end subroutine

  !----------------------------------------------------------------------------
  !> Return the actual number of potential values.
  integer function getNumPotentialValues(potential)
    implicit none
    type (PotentialData), intent(in) :: potential
    getNumPotentialValues = (potential%irmd+(potential%irnsd+1) * &
                            (potential%lmpot-1)) &
                            * potential%nspin
  end function

end module PotentialData_mod
