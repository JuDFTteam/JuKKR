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

  type EnergyResults
    double precision :: vbc(2) !< new muffin-tin zero ???
    double precision :: epotin !< kinetic energy minus sum of single particle energies
    double precision :: vmad !< Madelung potential - not used anymore?
    double precision, allocatable :: ecou(:) !< Coulomb energies
    double precision, allocatable :: espc(:,:) !< core energies
    double precision, allocatable :: espv(:,:) !< E valence bands
    double precision, allocatable :: exc(:) !< XC-energy
    double precision, allocatable :: ac_madelung(:) !< TODO: should not be field of this structure

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
  subroutine createEnergyResults(self, nspind, lmaxd)
    type(EnergyResults), intent(inout) :: self
    integer, intent(in) :: nspind
    integer, intent(in) :: lmaxd

    integer :: memory_stat

    self%lpot = 2 * lmaxd
    self%nspind = nspind
    self%lmaxd = lmaxd
    self%epotin = 0.d0
    self%vbc = 0.d0
    self%vmad = 0.d0

    ALLOCATECHECK(self%ecou(0:self%lpot))
    ALLOCATECHECK(self%espc(0:3,nspind))
    ALLOCATECHECK(self%espv(0:lmaxd+1,nspind))
    ALLOCATECHECK(self%exc(0:self%lpot))
    ALLOCATECHECK(self%ac_madelung((self%lpot+1)**2))

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

    DEALLOCATECHECK(self%ecou)
    DEALLOCATECHECK(self%espc)
    DEALLOCATECHECK(self%espv)
    DEALLOCATECHECK(self%exc)
    DEALLOCATECHECK(self%ac_madelung)
  endsubroutine ! destroy

endmodule
