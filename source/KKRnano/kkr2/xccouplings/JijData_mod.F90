! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATEIGNORE(X) deallocate(X, stat=ist)

module JijData_mod
  implicit none
  private
  public :: JijData, create, destroy
  public :: createJijData, destroyJijData ! deprecated
  
  type JijData

    double complex, allocatable :: gsxij(:,:,:,:)

    logical :: do_jij_calculation
    double precision :: rcutJij
    integer :: nxij
    double precision, allocatable  :: rxij(:) !< interatomic distance Ri-Rj
    double precision, allocatable  :: rxccls(:,:)  !< position relative of j rel. to i (sorted)
    double precision, allocatable  :: zkrxij(:,:,:) !< set up in clsjij, used in kkrmat01
    integer, allocatable  :: ixcp(:) !< index to atom in elem/cell at site in cluster
    integer, allocatable  :: nxcp(:) !< index to bravais lattice at site in cluster
    double complex, allocatable  :: jxcijint(:) !< integrated Jij
    double complex, allocatable  :: dtixij(:,:,:)
    double complex, allocatable  :: gmatxij(:,:,:,:)

    integer :: nxijd
    integer :: lmmaxd
    integer :: nspind
    integer :: active_spin ! spin-direction that is currently calculated

  endtype

  
  interface create
    module procedure createJijData
  endinterface
  
  interface destroy
    module procedure destroyJijData
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a JijData object.
  !> @param[inout] self    The JijData object to construct.
  !> @param[in]    nxijd
  !> @param[in]    lmmaxd
  !> @param[in]    nspind
  subroutine createJijData(self, do_jij_calculation, rcutjij, nxijd, lmmaxd, nspind)
    type(JijData), intent(inout) :: self
    logical, intent(in) ::  do_jij_calculation
    double precision, intent(in) :: rcutjij
    integer, intent(in) ::  nxijd
    integer, intent(in) ::  lmmaxd
    integer, intent(in) ::  nspind

    double complex, parameter :: CZERO = (0.d0, 0.d0)
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
    ALLOCATECHECK(self%gsxij(lmmaxd,lmmaxd,NSYMAXD,NXIJD))

    self%rxij = 0.d0
    self%rxccls = 0.d0
    self%zkrxij = 0.d0
    self%ixcp = 0
    self%nxcp = 0
    self%jxcijint = CZERO
    self%dtixij = CZERO
    self%gmatxij = CZERO
    self%nxij = 0
    self%active_spin = 1
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a JijData object.
  !> @param[inout] self    The JijData object to destroy.
  elemental subroutine destroyJijData(self)
    type(JijData), intent(inout) :: self

    integer :: ist ! ignore status

    DEALLOCATEIGNORE(self%rxij)
    DEALLOCATEIGNORE(self%rxccls)
    DEALLOCATEIGNORE(self%zkrxij)
    DEALLOCATEIGNORE(self%ixcp)
    DEALLOCATEIGNORE(self%nxcp)
    DEALLOCATEIGNORE(self%jxcijint)
    DEALLOCATEIGNORE(self%dtixij)
    DEALLOCATEIGNORE(self%gmatxij)
    DEALLOCATEIGNORE(self%gsxij)
  endsubroutine ! destroy
  
endmodule
