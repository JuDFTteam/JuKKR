! TODO: initialisation

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module LDAUData_mod
  implicit none
  private
  public :: LDAUData, create, destroy

  type LDAUData
    double precision :: EULDAU
    double precision :: EDCLDAU
    double precision :: EREFLDAU
    integer :: lmaxd
    integer :: mmaxd
    logical :: ldau
    integer :: nldau
    double complex, allocatable :: phildau(:,:) 
    double complex, allocatable :: dmatldau(:,:,:,:) 
    double precision, allocatable :: uldau(:) 
    double precision, allocatable :: jldau(:) 
    double precision, allocatable :: umldau(:,:,:,:,:) 
    double precision, allocatable :: wmldau(:,:,:,:) 
    integer, allocatable :: lldau(:)

    integer :: irmd
    integer :: lmaxd1
    integer :: nspind
  end type LDAUData

  interface create
    module procedure createLDAUData
  endinterface
  
  interface destroy
    module procedure destroyLDAUData
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a LDAUData object.
  !> @param[inout] self    The LDAUData object to construct.
  !> @param[in]    irmd
  !> @param[in]    lmaxd1
  !> @param[in]    mmaxd
  !> @param[in]    nspind
  subroutine createLDAUData(self, ldau, irmd, lmaxd, nspind)
    type (LDAUData), intent(inout) :: self
    integer, intent(in) ::  irmd
    integer, intent(in) ::  nspind
    integer, intent(in) :: lmaxd
    logical, intent(in) :: ldau

    integer :: mmaxd
    integer :: lmaxd1

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)

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

    self%phildau = CZERO
    self%dmatldau = CZERO
    self%uldau = 0.0d0
    self%jldau = 0.0d0
    self%umldau = 0.0d0
    self%wmldau = 0.0d0
    self%lldau = 0

  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a LDAUData object.
  !> @param[inout] self    The LDAUData object to destroy.
  subroutine destroyLDAUData(self)
    type (LDAUData), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%phildau)
    DEALLOCATECHECK(self%dmatldau)
    DEALLOCATECHECK(self%uldau)
    DEALLOCATECHECK(self%jldau)
    DEALLOCATECHECK(self%umldau)
    DEALLOCATECHECK(self%wmldau)
    DEALLOCATECHECK(self%lldau)
  endsubroutine ! destroy

endmodule ! LDAUData_mod
