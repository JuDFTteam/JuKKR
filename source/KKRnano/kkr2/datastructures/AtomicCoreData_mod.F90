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
    double precision, allocatable :: RHOCAT(:,:) !< radial charge density of core
  endtype

  interface create
    module procedure createAtomicCoreData
  endinterface
  
  interface destroy
    module procedure destroyAtomicCoreData
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  subroutine createAtomicCoreData(core, irmd)
    type(AtomicCoreData), intent(inout) :: core
    integer, intent(in) :: irmd

    core%irmd = irmd

    ! initialise with garbage values
    core%LCORE = -1
    core%NCORE = -1
    core%ECORE = 1.d9 ! much to high to be reasonable
    core%ITITLE = 0
    core%QC_corecharge = 0.d0

    allocate(core%RHOCAT(irmd, 2)) ! allocate for both spin directions
    core%RHOCAT = 0.d0

  endsubroutine ! create


  !----------------------------------------------------------------------------
  subroutine destroyAtomicCoreData(core)
    type(AtomicCoreData), intent(inout) :: core

    deallocate(core%RHOCAT)

  endsubroutine ! destroy

endmodule AtomicCoreData_mod
