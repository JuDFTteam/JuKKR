! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module LDAUData_mod

  type LDAUData
    double precision::EULDAU
    double precision::EDCLDAU
    double precision::EREFLDAU
    integer  :: lmaxd
    integer  :: mmaxd
    logical  :: ldau
    integer  :: nldau
    double complex , allocatable, dimension(:,:)  :: phildau
    double complex , allocatable, dimension(:,:,:,:)  :: dmatldau
    double precision , allocatable, dimension(:)  :: uldau
    double precision , allocatable, dimension(:)  :: jldau
    double precision , allocatable, dimension(:,:,:,:,:)  :: umldau
    double precision , allocatable, dimension(:,:,:,:)  :: wmldau
    integer , allocatable, dimension(:)  :: lldau

    integer :: irmd
    integer :: lmaxd1
    integer :: nspind
  end type LDAUData

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a LDAUData object.
  !> @param[inout] self    The LDAUData object to construct.
  !> @param[in]    irmd
  !> @param[in]    lmaxd1
  !> @param[in]    mmaxd
  !> @param[in]    nspind
  subroutine createLDAUData(self, ldau, irmd, lmaxd, nspind)
    implicit none
    type (LDAUData), intent(inout) :: self
    integer, intent(in) ::  irmd
    integer, intent(in) ::  nspind
    integer, intent(in) :: lmaxd
    logical, intent(in) :: ldau

    integer :: mmaxd
    integer :: lmaxd1

    integer :: memory_stat

    lmaxd1 = lmaxd + 1
    mmaxd = 2* lmaxd + 1

    self%irmd = irmd
    self%lmaxd1 = lmaxd1
    self%mmaxd = mmaxd
    self%nspind = nspind
    self%ldau = ldau

    ALLOCATECHECK(self%phildau(irmd,lmaxd1))
    ALLOCATECHECK(self%dmatldau(mmaxd,mmaxd,nspind,lmaxd1))
    ALLOCATECHECK(self%uldau(lmaxd1))
    ALLOCATECHECK(self%jldau(lmaxd1))
    ALLOCATECHECK(self%umldau(mmaxd,mmaxd,mmaxd,mmaxd,lmaxd1))
    ALLOCATECHECK(self%wmldau(mmaxd,mmaxd,lmaxd1,nspind))
    ALLOCATECHECK(self%lldau(lmaxd1))

    self%NLDAU = 0
    self%EULDAU = 0.0d0
    self%EDCLDAU = 0.0d0
    self%EREFLDAU = 0.0d0

  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a LDAUData object.
  !> @param[inout] self    The LDAUData object to destroy.
  subroutine destroyLDAUData(self)
    implicit none
    type (LDAUData), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%phildau)
    DEALLOCATECHECK(self%dmatldau)
    DEALLOCATECHECK(self%uldau)
    DEALLOCATECHECK(self%jldau)
    DEALLOCATECHECK(self%umldau)
    DEALLOCATECHECK(self%wmldau)
    DEALLOCATECHECK(self%lldau)
  end subroutine

end module
