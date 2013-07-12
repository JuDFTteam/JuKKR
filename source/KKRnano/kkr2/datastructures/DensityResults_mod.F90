! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module DensityResults_mod

  !> Contains densities and integrated densities (charges)
  !> for one specific atom.
  type DensityResults
    !double precision  :: DENEF
    double precision  :: total_charge_neutrality
    double precision , allocatable, dimension(:,:,:)  :: R2NEF
    double precision , allocatable, dimension(:,:,:)  :: RHO2NS
    double precision , allocatable, dimension(:,:)  :: CHARGE
    double precision , allocatable, dimension(:)  :: CMINST
    double precision , allocatable, dimension(:)  :: CMOM
    double precision , allocatable, dimension(:)  :: CATOM
    double complex , allocatable, dimension(:,:,:)  :: DEN

    double precision :: force_flm(-1:1)

    !> Renormalisation factors LLoyd - leave here?
    double precision , allocatable, dimension(:,:)  :: RNORM

    integer :: irmd
    integer :: lmpotd
    integer :: lmaxd
    integer :: iemxd
    integer :: nspind
  end type DensityResults

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a DensityResults object.
  !> @param[in,out] self    The DensityResults object to construct.
  !> @param[in]     dims
  subroutine createDensityResults(self, dims, num_radial_irmd)
    use DimParams_mod
    implicit none
    type (DensityResults), intent(inout) :: self
    type (DimParams), intent(in) :: dims
    integer, intent(in) :: num_radial_irmd

    call createDensityResultsImpl(self, num_radial_irmd, dims%lmpotd, dims%lmaxd, &
                                  dims%iemxd, dims%nspind)

  end subroutine

  !-----------------------------------------------------------------------------
  !> Constructs a DensityResults object.
  !> @param[in,out] self    The DensityResults object to construct.
  !> @param[in]     irmd
  !> @param[in]     lmpotd
  !> @param[in]     lmaxd
  !> @param[in]     iemxd
  !> @param[in]     nspind
  subroutine createDensityResultsImpl(self, irmd,lmpotd,lmaxd,iemxd,nspind)
    implicit none
    type (DensityResults), intent(inout) :: self
    integer, intent(in) ::  irmd
    integer, intent(in) ::  lmpotd
    integer, intent(in) ::  lmaxd
    integer, intent(in) ::  iemxd
    integer, intent(in) ::  nspind

    integer :: memory_stat

    self%irmd = irmd
    self%lmpotd = lmpotd
    self%lmaxd = lmaxd
    self%iemxd = iemxd
    self%nspind = nspind

    ALLOCATECHECK(self%R2NEF(irmd,lmpotd,2))
    ALLOCATECHECK(self%RHO2NS(irmd,lmpotd,2))
    ALLOCATECHECK(self%CHARGE(0:lmaxd+1,2))
    ALLOCATECHECK(self%CMINST(lmpotd))
    ALLOCATECHECK(self%CMOM(lmpotd))
    ALLOCATECHECK(self%CATOM(nspind))
    ALLOCATECHECK(self%DEN(0:lmaxd+1,iemxd,nspind))
    ALLOCATECHECK(self%RNORM(IEMXD,2))

    !initialise to be safe
    self%RHO2NS = 0.0d0
    self%R2NEF = 0.0d0
    self%force_flm = 9999.9d0
  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a DensityResults object.
  !> @param[in,out] self    The DensityResults object to destroy.
  subroutine destroyDensityResults(self)
    implicit none
    type (DensityResults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%R2NEF)
    DEALLOCATECHECK(self%RHO2NS)
    DEALLOCATECHECK(self%CHARGE)
    DEALLOCATECHECK(self%CMINST)
    DEALLOCATECHECK(self%CMOM)
    DEALLOCATECHECK(self%CATOM)
    DEALLOCATECHECK(self%DEN)
    DEALLOCATECHECK(self%RNORM)
  end subroutine

end module
