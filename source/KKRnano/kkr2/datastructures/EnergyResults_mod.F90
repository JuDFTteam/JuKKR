! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

!> Energy contributions of one atom.
module EnergyResults_mod
  implicit none
  private
  public :: EnergyResults, create, destroy
  public :: createEnergyResults, destroyEnergyResults ! deprecated

  type EnergyResults
    double precision   :: VBC(2) !< new muffin-tin zero ???
    double precision , allocatable :: ECOU(:) !< Coulomb energies
    double precision , allocatable :: ESPC(:,:) !< core energies
    double precision , allocatable :: ESPV(:,:) !< E valence bands
    double precision , allocatable :: EXC(:) !< XC-energy
    double precision  :: EPOTIN !< kinetic energy minus sum of single particle energies
    double precision  :: VMAD !< Madelung potential - not used anymore?
    double precision , allocatable :: AC_madelung(:) !< TODO: should not be field of this structure

    double precision :: e_vxc   !< XC-part of the double counting energy
    double precision :: e_shift !< energy change due to constant potential shift ("MT-shift")
    double precision :: e_madelung  !< Madelung energy minus a part that cancels within reference sphere
    double precision :: e_total(2)

    integer :: lpot
    integer :: nspind
    integer :: lmaxd
  endtype ! EnergyResults


  interface create
    module procedure createEnergyResults
  endinterface
  
  interface destroy
    module procedure destroyEnergyResults
  endinterface

  contains

  !-----------------------------------------------------------------------------
  !> Constructs a EnergyResults object.
  !> @param[in,out] self    The EnergyResults object to construct.
  !> @param[in]     nspind
  !> @param[in]     lmaxd
  subroutine createEnergyResults(self, nspind,lmaxd)
    type(EnergyResults), intent(inout) :: self
    integer, intent(in) ::  nspind
    integer, intent(in) ::  lmaxd

    integer :: memory_stat

    self%lpot = 2 * lmaxd
    self%nspind = nspind
    self%lmaxd = lmaxd
    self%EPOTIN = 0.d0
    self%VBC = 0.d0
    self%VMAD = 0.d0

    ALLOCATECHECK(self%ECOU(0:self%lpot))
    ALLOCATECHECK(self%ESPC(0:3,nspind))
    ALLOCATECHECK(self%ESPV(0:lmaxd+1,nspind))
    ALLOCATECHECK(self%EXC(0:self%lpot))
    ALLOCATECHECK(self%AC_madelung((self%lpot+1)**2))

    self%ecou = 0.d0
    self%espc = 0.d0
    self%exc  = 0.d0
    self%espv = 0.d0

    self%e_vxc = 0.d0
    self%e_shift = 0.d0
    self%e_madelung = 0.d0
    self%e_total = 0.d0

  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a EnergyResults object.
  !> @param[in,out] self    The EnergyResults object to destroy.
  subroutine destroyEnergyResults(self)
    type(EnergyResults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%ECOU)
    DEALLOCATECHECK(self%ESPC)
    DEALLOCATECHECK(self%ESPV)
    DEALLOCATECHECK(self%EXC)
    DEALLOCATECHECK(self%AC_madelung)
  endsubroutine ! destroy

endmodule
