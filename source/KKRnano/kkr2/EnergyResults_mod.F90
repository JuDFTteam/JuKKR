! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

!> Energy contributions of one atom.
module EnergyResults_mod

  type EnergyResults
    double precision , dimension(2)  :: VBC
    double precision , allocatable, dimension(:)  :: ECOU
    double precision , allocatable, dimension(:,:)  :: ESPC
    double precision , allocatable, dimension(:,:)  :: ESPV
    double precision , allocatable, dimension(:)  :: EXC
    double precision  :: EPOTIN

    integer :: lpot
    integer :: nspind
    integer :: lmaxd
  end type EnergyResults

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a EnergyResults object.
  !> @param[inout] self    The EnergyResults object to construct.
  !> @param[in]    nspind
  !> @param[in]    lmaxd
  subroutine createEnergyResults(self, nspind,lmaxd)
    implicit none
    type (EnergyResults), intent(inout) :: self
    integer, intent(in) ::  nspind
    integer, intent(in) ::  lmaxd

    integer :: memory_stat

    self%lpot = 2 * lmaxd
    self%nspind = nspind
    self%lmaxd = lmaxd
    self%EPOTIN = 0.0d0
    self%VBC = 0.0d0

    ALLOCATECHECK(self%ECOU(0:self%lpot))
    ALLOCATECHECK(self%ESPC(0:3,nspind))
    ALLOCATECHECK(self%ESPV(0:lmaxd+1,nspind))
    ALLOCATECHECK(self%EXC(0:self%lpot))
  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a EnergyResults object.
  !> @param[inout] self    The EnergyResults object to destroy.
  subroutine destroyEnergyResults(self)
    implicit none
    type (EnergyResults), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%ECOU)
    DEALLOCATECHECK(self%ESPC)
    DEALLOCATECHECK(self%ESPV)
    DEALLOCATECHECK(self%EXC)
  end subroutine

end module
