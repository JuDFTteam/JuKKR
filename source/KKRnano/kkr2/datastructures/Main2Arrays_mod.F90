!> Datastructure that contains various arrays that
!> could not be factored out (yet).
!>
!> Includes k-mesh, basis atom positions, atomic numbers, symmetry matrices,
!> reference potential strength, Bravais vectors
!>
!> Note: Scales as O(N**2) in space and time.


! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module Main2Arrays_mod
  implicit none
  private
  public :: Main2Arrays, create, destroy
  public :: createMain2Arrays!, destroyMain2Arrays ! deprecated
  public :: readMain2Arrays, writeMain2Arrays

  type Main2Arrays
    integer :: nsymat
    integer :: lmmaxd
    integer :: isymindex(48)
    double complex,   allocatable :: dsymll(:,:,:)  !< tau symmetry matrices
    
    integer :: naez
    double precision :: bravais(3,3)
    double precision, allocatable :: rbasis(:,:)  !< basis atom positions
    double precision, allocatable :: zat(:)  !< atomic numbers
    
    integer :: kpoibz
    integer :: maxmesh
    integer :: maxmshd
    double precision, allocatable :: bzkp(:,:,:)  !< kpoints for each mesh
    double precision, allocatable :: volcub(:,:)  !< kpoint weights
    double precision, allocatable :: volbz(:)     !< bz volume?
    integer,          allocatable :: nofks(:)     !< number of k points for each mesh
    
!     integer :: iemxd
!     integer, allocatable :: kmesh(:) !< mapping of e-points to k-meshes

  endtype ! Main2Arrays

  interface create
    module procedure createMain2Arrays
  endinterface
  
  interface destroy
    module procedure destroyMain2Arrays
  endinterface

  contains

  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object.
  !> @param[in,out] self    The Main2Arrays object to construct.
  !> @param[in]    dims    Dimension parameters
  subroutine createMain2Arrays(self, dims)
    use DimParams_mod, only: DimParams
    type(Main2Arrays), intent(inout) :: self
    type(DimParams), intent(in) :: dims

    call createMain2ArraysImpl(self, &
!     dims%iemxd, &
    dims%lmmaxd, dims%naez, dims%kpoibz, dims%maxmshd)    
    
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object. (implementation, don't call directly!)
  !> @param[in,out] self    The Main2Arrays object to construct.
  !> @param[in]    iemxd
  !> @param[in]    lmmaxd
  !> @param[in]    naez
  !> @param[in]    kpoibz
  !> @param[in]    maxmshd
  subroutine createMain2ArraysImpl(self, &
!   iemxd, &
  lmmaxd, naez, kpoibz, maxmshd)
    use, intrinsic :: ieee_features
    use, intrinsic :: ieee_arithmetic
    
    type(Main2Arrays), intent(inout) :: self
!     integer, intent(in) :: iemxd
    integer, intent(in) :: lmmaxd, naez, kpoibz, maxmshd
    
    integer :: memory_stat
    double precision :: nan

    self%nsymat = 0
    self%maxmesh = 0

!     self%iemxd = iemxd
    self%lmmaxd = lmmaxd
    self%naez = naez
    self%kpoibz = kpoibz
    self%maxmshd = maxmshd

    ALLOCATECHECK(self%dsymll(lmmaxd,lmmaxd,48))
    ALLOCATECHECK(self%bzkp(3,kpoibz,maxmshd))
    ALLOCATECHECK(self%volcub(kpoibz,maxmshd))
    ALLOCATECHECK(self%volbz(maxmshd))
!     ALLOCATECHECK(self%kmesh(iemxd))
    ALLOCATECHECK(self%nofks(maxmshd))

    ALLOCATECHECK(self%rbasis(3,naez))
    ALLOCATECHECK(self%zat(naez))
    
    nan = ieee_value(nan, ieee_signaling_nan)
    self%dsymll = nan
    self%rbasis = nan
    self%bzkp = nan
    self%volcub = nan
    self%volbz = nan
!     self%kmesh = -99
    self%nofks = -99
    self%zat = nan

  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a Main2Arrays object.
  !> @param[in,out] self    The Main2Arrays object to destroy.
  subroutine destroyMain2Arrays(self)
    type(Main2Arrays), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%dsymll)
    DEALLOCATECHECK(self%rbasis)
    DEALLOCATECHECK(self%bzkp)
    DEALLOCATECHECK(self%volcub)
    DEALLOCATECHECK(self%volbz)
!     DEALLOCATECHECK(self%kmesh)
    DEALLOCATECHECK(self%nofks)
    DEALLOCATECHECK(self%zat)

  endsubroutine ! destroy

  !-----------------------------------------------------------------------------
  !> Writes Main2Arrays data to file.
  !> @param[in] self    The Main2Arrays object to write.
  subroutine writeMain2Arrays(self, filename)
    type(Main2Arrays), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 67

    open (fu, file=filename, form='unformatted')
    write(fu) self%bravais, &
              self%isymindex, &
              self%dsymll, &
              self%rbasis, &
              self%bzkp, &
              self%volcub, &
              self%volbz, &
!               self%kmesh, &
              self%nofks, &
              self%zat, &
              self%nsymat, &  ! write some scalars too
              self%maxmesh
    close(fu)

  endsubroutine ! write

  !-----------------------------------------------------------------------------
  !> Reads Main2Arrays data from file.
  !> @param[in,out] self    The Main2Arrays object to read.
  subroutine readMain2Arrays(self, filename)
    type(Main2Arrays), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 67

    open (fu, file=filename, form='unformatted')
    read (fu) self%bravais, &
              self%isymindex, &
              self%dsymll, &
              self%rbasis, &
              self%bzkp, &
              self%volcub, &
              self%volbz, &
!               self%kmesh, &
              self%nofks, &
              self%zat, &
              self%nsymat, & ! write some scalars too
              self%maxmesh
    close(fu)

  endsubroutine ! read

endmodule ! Main2Arrays_mod
