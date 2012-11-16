! eliminates LCORE, NCORE, ITITLE, ECORE, LCOREMAX, QC, NPOTD

module AtomicCoreData_mod

  implicit none

  type AtomicCoreData
    integer :: LCOREMAX
    integer :: LCORE_atom(20,2)       !< for historical reasons always 2 spin directions
    integer :: NCORE_atom(2)
    double precision :: ECORE(20,2)   !< first dim 20: max. 20 core states
    integer :: ITITLE(20,2)           !< potential title as old school integer string
    double precision :: QC_corecharge !< total charge of core electrons ! remove??
    double precision, dimension(:,:), allocatable ::  RHOCAT !< radial charge density of core
    integer :: irmd
    double precision :: Z_nuclear

  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createAtomicCoreData(core, Z_nuclear, irmd)
    implicit none
    type (AtomicCoreData), intent(inout) :: core
    double precision, intent(in) :: Z_nuclear
    integer, intent(in) :: irmd

    core%irmd = irmd
    core%Z_nuclear = Z_nuclear

    ! initialise with garbage values
    core%LCOREMAX = -1
    core%LCORE_atom = -1
    core%NCORE_atom = -1
    core%ECORE = 1d9
    core%ITITLE = 0
    core%QC_corecharge = 0.0d0

    allocate(core%RHOCAT(irmd, 2)) ! allocate for both spin directions
    core%RHOCAT = 0.0d0

  end subroutine


  !----------------------------------------------------------------------------
  subroutine destroyAtomicCoreData(core)
    implicit none
    type (AtomicCoreData), intent(inout) :: core

    deallocate(core%RHOCAT)

  end subroutine

end module AtomicCoreData_mod
