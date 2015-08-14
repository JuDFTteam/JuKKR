module AtomicCoreData_mod
  implicit none
  private
  public :: AtomicCoreData, create, destroy 
  public :: createAtomicCoreData, destroyAtomicCoreData ! deprecated

  !> Structure that contains information about core states
  type AtomicCoreData
    integer :: LCORE(20,2)       !< for historical reasons always 2 spin directions
    integer :: NCORE(2)          !< number of core states
    double precision :: ECORE(20,2)   !< first dim 20: max. 20 core states
    integer :: ITITLE(20,2)           !< potential title as old school integer string
    integer :: irmd

    ! data used for calculation - not written to disk
    double precision :: QC_corecharge !< total charge of core electrons
    double precision, dimension(:,:), allocatable ::  RHOCAT !< radial charge density of core

  end type

  interface create
    module procedure createAtomicCoreData
  endinterface
  
  interface destroy
    module procedure destroyAtomicCoreData
  endinterface
  
  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createAtomicCoreData(core, irmd)
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
    type (AtomicCoreData), intent(inout) :: core

    deallocate(core%RHOCAT)

  end subroutine

end module AtomicCoreData_mod
