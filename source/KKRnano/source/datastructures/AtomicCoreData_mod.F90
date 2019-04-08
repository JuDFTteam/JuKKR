module AtomicCoreData_mod
  implicit none
  private
  public :: AtomicCoreData, create, destroy

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
  endtype ! AtomicCoreData

  interface create
    module procedure createAtomicCoreData
  endinterface
  
  interface destroy
    module procedure destroyAtomicCoreData
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  subroutine createAtomicCoreData(self, irmd)
    type(AtomicCoreData), intent(inout) :: self
    integer, intent(in) :: irmd

    self%irmd = irmd

    ! initialise with garbage values
    self%LCORE = -1
    self%NCORE = -1
    self%ECORE = 1.d9 ! much to high to be reasonable
    self%ITITLE = 0
    self%QC_corecharge = 0.d0

    allocate(self%RHOCAT(irmd, 2)) ! allocate for both spin directions
    self%RHOCAT = 0.d0

  endsubroutine ! create


  !----------------------------------------------------------------------------
  elemental subroutine destroyAtomicCoreData(self)
    type(AtomicCoreData), intent(inout) :: self
    
    integer :: ist
    deallocate(self%RHOCAT, stat=ist)

  endsubroutine ! destroy

endmodule ! AtomicCoreData_mod
