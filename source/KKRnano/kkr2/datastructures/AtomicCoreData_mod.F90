! how to deal with spin???

module AtomicCoreData_mod

  implicit none

  type AtomicCoreData
    integer :: LCORE(20,2)       !< for historical reasons always 2 spin directions
    integer :: NCORE(2)
    double precision :: ECORE(20,2)   !< first dim 20: max. 20 core states
    integer :: ITITLE(20,2)           !< potential title as old school integer string
    double precision :: QC_corecharge !< total charge of core electrons
    double precision, dimension(:,:), allocatable ::  RHOCAT !< radial charge density of core
    integer :: irmd
    !!!double precision :: Z_nuclear

  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createAtomicCoreData(core, irmd)
    implicit none
    type (AtomicCoreData), intent(inout) :: core
    integer, intent(in) :: irmd

    core%irmd = irmd

    ! initialise with garbage values
    core%LCORE = -1
    core%NCORE = -1
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
