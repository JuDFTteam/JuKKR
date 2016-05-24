! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module DensityResults_mod
implicit none
  private
  public :: DensityResults, create, destroy

  !> Contains densities and integrated densities (charges)
  !> for one specific atom.
  type DensityResults

    double precision  :: total_charge_neutrality     !< excess charge in system
    double precision, allocatable :: r2nef(:,:,:)    !< density of states at fermi energy
    double precision, allocatable :: rho2ns(:,:,:)   !< l-expanded density times radius**2
    double precision, allocatable :: charge(:,:)     !< charge(:,1) l-resolved charge, charge(:,2) l-resolved mag. moment
    double precision :: chrgsemicore_per_atom        !< semicore charge per atom
    double precision, allocatable  :: cminst(:)      !< charge moments in interstitial
    double precision, allocatable  :: cmom(:)        !< charge moments in muffin-tin sphere
    double precision, allocatable :: catom(:)        !< catom(1) = total charge in cell, catom(2) = total mag. moment in cell
    double complex, allocatable :: den(:,:,:)        !< complex density of states
    double precision :: force_flm(-1:1)
    double precision, allocatable :: rnorm(:,:)      !> renormalisation factors lloyd - leave here?

    integer :: irmd
    integer :: lmpotd
    integer :: lmaxd
    integer :: iemxd
    integer :: nspind
  endtype ! DensityResults

  
  interface create
    module procedure createDensityResults
  endinterface
  
  interface destroy
    module procedure destroyDensityResults
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a DensityResults object.
  !> @param[in,out] self    The DensityResults object to construct.
  !> @param[in]     irmd
  !> @param[in]     lmpotd
  !> @param[in]     lmaxd
  !> @param[in]     iemxd
  !> @param[in]     nspind
  subroutine createDensityResults(self, lmpotd, lmaxd, iemxd, nspind, irmd)
    type(DensityResults), intent(inout) :: self
    integer, intent(in) ::  lmpotd
    integer, intent(in) ::  lmaxd
    integer, intent(in) ::  iemxd
    integer, intent(in) ::  nspind
    integer, intent(in) ::  irmd

    integer :: memory_stat

    self%irmd = irmd
    self%lmpotd = lmpotd
    self%lmaxd = lmaxd
    self%iemxd = iemxd
    self%nspind = nspind

    ALLOCATECHECK(self%r2nef(irmd,lmpotd,2))
    ALLOCATECHECK(self%rho2ns(irmd,lmpotd,2))
    ALLOCATECHECK(self%charge(0:lmaxd+1,2))
    ALLOCATECHECK(self%cminst(lmpotd))
    ALLOCATECHECK(self%cmom(lmpotd))
    ALLOCATECHECK(self%catom(nspind))
    ALLOCATECHECK(self%den(0:lmaxd+1,iemxd,nspind))
    ALLOCATECHECK(self%rnorm(iemxd,2))

    !initialise to be safe
    self%rho2ns = 0.0d0
    self%r2nef = 0.0d0
    self%force_flm = 9999.9d0
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a DensityResults object.
  !> @param[in,out] self    The DensityResults object to destroy.
  subroutine destroyDensityResults(self)
    type(DensityResults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%r2nef)
    DEALLOCATECHECK(self%rho2ns)
    DEALLOCATECHECK(self%charge)
    DEALLOCATECHECK(self%cminst)
    DEALLOCATECHECK(self%cmom)
    DEALLOCATECHECK(self%catom)
    DEALLOCATECHECK(self%den)
    DEALLOCATECHECK(self%rnorm)
  endsubroutine ! destroy

endmodule ! DensityResults_mod
