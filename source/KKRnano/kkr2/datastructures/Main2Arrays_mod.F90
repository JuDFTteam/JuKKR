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
  public :: Main2Arrays, create, destroy, load, store

  type Main2Arrays
    integer :: nsymat
    integer :: lmmaxd
    integer :: isymindex(48)
    double complex,   allocatable :: dsymll(:,:,:)  !< tau symmetry matrices
    
    integer :: naez
    double precision :: bravais(3,3)
    double precision, allocatable :: rbasis(:,:)    !< basis atom positions rbasis(3,naez)
    double precision, allocatable :: zat(:)         !< atomic numbers zat(naez)

    integer :: kpoibz
    integer :: maxmesh
    integer :: maxmshd
    double precision, allocatable :: bzkp(:,:,:)    !< kpoints for each mesh
    double precision, allocatable :: volcub(:,:)    !< kpoint weights
    double precision, allocatable :: volbz(:)       !< bz volume?
    integer,          allocatable :: nofks(:)       !< number of k points for each mesh

  endtype ! Main2Arrays

  interface create
    module procedure createMain2Arrays
  endinterface
  
  interface destroy
    module procedure destroyMain2Arrays
  endinterface

  interface load
    module procedure readMain2Arrays
  endinterface
  
  interface store
    module procedure writeMain2Arrays
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

    call createMain2ArraysImpl(self, dims%lmmaxd, dims%naez, dims%kpoibz, dims%maxmshd)    
    
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object. (implementation, don't call directly!)
  !> @param[in,out] self    The Main2Arrays object to construct.
  !> @param[in]    iemxd
  !> @param[in]    lmmaxd
  !> @param[in]    naez
  !> @param[in]    kpoibz
  !> @param[in]    maxmshd
  subroutine createMain2ArraysImpl(self, lmmaxd, naez, kpoibz, maxmshd)
    type(Main2Arrays), intent(inout) :: self
    integer, intent(in) :: lmmaxd, naez, kpoibz, maxmshd
    
    integer :: memory_stat

    self%nsymat = 0
    self%maxmesh = 0

    self%lmmaxd = lmmaxd
    self%naez = naez
    self%kpoibz = kpoibz
    self%maxmshd = maxmshd

    ALLOCATECHECK(self%dsymll(lmmaxd,lmmaxd,48))
    ALLOCATECHECK(self%bzkp(3,kpoibz,maxmshd))
    ALLOCATECHECK(self%volcub(kpoibz,maxmshd))
    ALLOCATECHECK(self%volbz(maxmshd))
    ALLOCATECHECK(self%nofks(maxmshd))

    ALLOCATECHECK(self%rbasis(3,naez))
    ALLOCATECHECK(self%zat(naez))
    
    self%dsymll = 0
    self%rbasis = 0
    self%bzkp = 0
    self%volcub = 0
    self%volbz = 0
    self%nofks = 0
    self%zat = 0
    
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

    open (fu, file=filename, form='unformatted', action='write')
    write(fu) self%bravais, &
              self%isymindex, &
              self%dsymll, &
              self%rbasis, &
              self%bzkp, &
              self%volcub, &
              self%volbz, &
              self%nofks, &
              self%zat, &
              self%nsymat, &  ! write some scalars too
              self%maxmesh
    close(fu)

  endsubroutine ! store

  !-----------------------------------------------------------------------------
  !> Reads Main2Arrays data from file.
  !> @param[in,out] self    The Main2Arrays object to read.
  subroutine readMain2Arrays(self, filename)
    type(Main2Arrays), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 67

    open (fu, file=filename, form='unformatted', action='read')
    read (fu) self%bravais, &
              self%isymindex, &
              self%dsymll, &
              self%rbasis, &
              self%bzkp, &
              self%volcub, &
              self%volbz, &
              self%nofks, &
              self%zat, &
              self%nsymat, & ! write some scalars too
              self%maxmesh
    close(fu)

  endsubroutine ! load

endmodule ! Main2Arrays_mod
