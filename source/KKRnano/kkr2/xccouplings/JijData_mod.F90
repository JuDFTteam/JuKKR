! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module JijData_mod

  type JijData

  !double precision, dimension(:), allocatable :: RXIJ          ! interatomic distance Ri-Rj
  !double precision, dimension(:,:), allocatable :: RXCCLS      ! position relative of j rel. to i (sorted)
  !double precision, dimension(:,:,:), allocatable :: ZKRXIJ    ! set up in clsjij, used in kkrmat01
  !integer, dimension(:), allocatable :: IXCP                   ! index to atom in elem/cell at site in cluster
  !integer, dimension(:), allocatable :: NXCP                   ! index to bravais lattice at site in cluster
  !double complex, dimension(:), allocatable ::  JXCIJINT       ! integrated Jij

    double complex , allocatable, dimension(:,:,:,:) :: GSXIJ

    logical :: do_jij_calculation
    double precision :: rcutJij
    integer :: nxij
    double precision , allocatable, dimension(:)  :: rxij
    double precision , allocatable, dimension(:,:)  :: rxccls
    double precision , allocatable, dimension(:,:,:)  :: zkrxij
    integer , allocatable, dimension(:)  :: ixcp
    integer , allocatable, dimension(:)  :: nxcp
    double complex , allocatable, dimension(:)  :: jxcijint
    double complex , allocatable, dimension(:,:,:)  :: dtixij
    double complex , allocatable, dimension(:,:,:,:)  :: gmatxij

    integer :: nxijd
    integer :: lmmaxd
    integer :: nspind

    integer :: active_spin ! spin-direction that is currently calculated

  end type JijData

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a JijData object.
  !> @param[inout] self    The JijData object to construct.
  !> @param[in]    nxijd
  !> @param[in]    lmmaxd
  !> @param[in]    nspind
  subroutine createJijData(self, do_jij_calculation, rcutjij, nxijd,lmmaxd,nspind)
    implicit none
    type (JijData), intent(inout) :: self
    logical, intent(in) ::  do_jij_calculation
    double precision, intent(in) :: rcutjij
    integer, intent(in) ::  nxijd
    integer, intent(in) ::  lmmaxd
    integer, intent(in) ::  nspind

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    integer, parameter :: NSYMAXD = 48

    integer :: memory_stat

    self%do_jij_calculation = do_jij_calculation
    self%rcutjij = rcutjij
    self%nxijd = nxijd
    self%lmmaxd = lmmaxd
    self%nspind = nspind

    ALLOCATECHECK(self%rxij(nxijd))
    ALLOCATECHECK(self%rxccls(3,nxijd))
    ALLOCATECHECK(self%zkrxij(NSYMAXD,3,nxijd))
    ALLOCATECHECK(self%ixcp(nxijd))
    ALLOCATECHECK(self%nxcp(nxijd))
    ALLOCATECHECK(self%jxcijint(nxijd))
    ALLOCATECHECK(self%dtixij(lmmaxd,lmmaxd,nspind))
    ALLOCATECHECK(self%gmatxij(lmmaxd,lmmaxd,nxijd,nspind))
    ALLOCATECHECK(self%GSXIJ(lmmaxd,lmmaxd,NSYMAXD,NXIJD))

    self%rxij = 0.0d0
    self%rxccls = 0.0d0
    self%zkrxij = 0.0d0
    self%ixcp = 0
    self%nxcp = 0
    self%jxcijint = CZERO
    self%dtixij = CZERO
    self%gmatxij = CZERO
    self%nxij = 0
    self%active_spin = 1
  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a JijData object.
  !> @param[inout] self    The JijData object to destroy.
  subroutine destroyJijData(self)
    implicit none
    type (JijData), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%rxij)
    DEALLOCATECHECK(self%rxccls)
    DEALLOCATECHECK(self%zkrxij)
    DEALLOCATECHECK(self%ixcp)
    DEALLOCATECHECK(self%nxcp)
    DEALLOCATECHECK(self%jxcijint)
    DEALLOCATECHECK(self%dtixij)
    DEALLOCATECHECK(self%gmatxij)
    DEALLOCATECHECK(self%GSXIJ)
  end subroutine
end module
