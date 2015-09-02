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
  public :: createMain2Arrays, destroyMain2Arrays ! deprecated
  public :: readMain2Arrays, writeMain2Arrays

  type Main2Arrays
    double precision :: bravais(3,3)
    integer :: isymindex(48)
    double complex, allocatable  :: DSYMLL(:,:,:)  !< tau symmetry matrices
    double precision, allocatable  :: RBASIS(:,:)  !< basis atom positions
    double precision, allocatable  :: BZKP(:,:,:)  !< kpoints for each mesh
    double precision, allocatable  :: VOLCUB(:,:)  !< kpoint weights
    double precision, allocatable  :: VOLBZ(:)     !< BZ volume?
    integer, allocatable  :: KMESH(:) !< mapping E-point to k-mesh
    integer, allocatable  :: NOFKS(:) !< number of k points for each mesh
    double precision, allocatable  :: ZAT(:)  !< atomic numbers
    double precision :: VREF !< repulsive screening pot. strength

    integer  :: NSYMAT
    integer  :: MAXMESH
    integer :: iemxd
    integer :: LMMAXD
    integer :: NAEZ
    integer :: KPOIBZ
    integer :: MAXMSHD

  endtype Main2Arrays
  
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

    call createMain2ArraysImpl(self, dims%lmaxd, dims%iemxd, dims%nspind, &
    dims%LMMAXD, dims%NAEZ, dims%LMXSPD, dims%KPOIBZ, dims%MAXMSHD, 0, &
    0, dims%nguessd, dims%ekmd, dims%smpid, dims%lpot, dims%IRMD, dims%LMPOTD)

  endsubroutine ! create



  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object. (implementation, don't call directly!)
  !> @param[in,out] self    The Main2Arrays object to construct.
  !> @param[in]    lmaxd
  !> @param[in]    iemxd
  !> @param[in]    nspind
  !> @param[in]    LMMAXD
  !> @param[in]    NAEZ
  !> @param[in]    LMXSPD
  !> @param[in]    IEMXD
  !> @param[in]    KPOIBZ
  !> @param[in]    MAXMSHD
  !> @param[in]    nrd
  !> @param[in]    NACLSD
  !> @param[in]    NSPIND
  !> @param[in]    nguessd
  !> @param[in]    lmmaxd
  !> @param[in]    ekmd
  !> @param[in]    smpid
  !> @param[in]    lpot
  !> @param[in]    IRMD
  !> @param[in]    LMPOTD
  subroutine createMain2ArraysImpl(self, lmaxd,iemxd,nspind,LMMAXD,NAEZ,LMXSPD,KPOIBZ,MAXMSHD,nrd,NACLSD,nguessd,ekmd,smpid,lpot,IRMD,LMPOTD)
    use, intrinsic :: ieee_features
    use, intrinsic :: ieee_arithmetic
    
    type(Main2Arrays), intent(inout) :: self
    integer, intent(in) ::  lmaxd
    integer, intent(in) ::  iemxd
    integer, intent(in) ::  nspind
    integer, intent(in) ::  LMMAXD
    integer, intent(in) ::  NAEZ
    integer, intent(in) ::  LMXSPD
    integer, intent(in) ::  KPOIBZ
    integer, intent(in) ::  MAXMSHD
    integer, intent(in) ::  nrd
    integer, intent(in) ::  NACLSD
    integer, intent(in) ::  nguessd
    integer, intent(in) ::  ekmd
    integer, intent(in) ::  smpid
    integer, intent(in) ::  lpot
    integer, intent(in) ::  IRMD
    integer, intent(in) ::  LMPOTD

    integer :: memory_stat
    double precision :: nan

    self%NSYMAT = 0
    self%MAXMESH = 0

!    Repulsive reference potential
!    in future: move as parameter to inputfile
    self%VREF = 8.0d0

    self%iemxd = iemxd
    self%LMMAXD = LMMAXD
    self%NAEZ = NAEZ
    self%IEMXD = IEMXD
    self%KPOIBZ = KPOIBZ
    self%MAXMSHD = MAXMSHD
    self%lmmaxd = lmmaxd

    ALLOCATECHECK(self%DSYMLL(LMMAXD,LMMAXD,48))
    ALLOCATECHECK(self%RBASIS(3,NAEZ))
    ALLOCATECHECK(self%BZKP(3,KPOIBZ,MAXMSHD))
    ALLOCATECHECK(self%VOLCUB(KPOIBZ,MAXMSHD))
    ALLOCATECHECK(self%VOLBZ(MAXMSHD))
    ALLOCATECHECK(self%KMESH(IEMXD))
    ALLOCATECHECK(self%NOFKS(MAXMSHD))
    ALLOCATECHECK(self%ZAT(NAEZ))

    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
    self%DSYMLL = nan
    self%rbasis = nan
    self%bzkp = nan
    self%volcub = nan
    self%volbz = nan
    self%kmesh = -99
    self%nofks = -99
    self%zat = nan

  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a Main2Arrays object.
  !> @param[in,out] self    The Main2Arrays object to destroy.
  subroutine destroyMain2Arrays(self)
    type(Main2Arrays), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%DSYMLL)
    DEALLOCATECHECK(self%RBASIS)
    DEALLOCATECHECK(self%BZKP)
    DEALLOCATECHECK(self%VOLCUB)
    DEALLOCATECHECK(self%VOLBZ)
    DEALLOCATECHECK(self%KMESH)
    DEALLOCATECHECK(self%NOFKS)
    DEALLOCATECHECK(self%ZAT)

  endsubroutine ! destroy

  !-----------------------------------------------------------------------------
  !> Writes Main2Arrays data to file.
  !> @param[in] self    The Main2Arrays object to write.
  subroutine writeMain2Arrays(self, filename)
    type(Main2Arrays), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 67

    open (fu, file=filename, form='unformatted')
    write(fu) self%BRAVAIS, &
              self%ISYMINDEX, &
              self%DSYMLL, &
              self%RBASIS, &
              self%BZKP, &
              self%VOLCUB, &
              self%VOLBZ, &
              self%KMESH, &
              self%NOFKS, &
              self%ZAT, &
              self%VREF, &
              self%NSYMAT, &  ! write some scalars too
              self%MAXMESH
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
    read (fu) self%BRAVAIS, &
              self%ISYMINDEX, &
              self%DSYMLL, &
              self%RBASIS, &
              self%BZKP, &
              self%VOLCUB, &
              self%VOLBZ, &
              self%KMESH, &
              self%NOFKS, &
              self%ZAT, &
              self%VREF, &
              self%NSYMAT, & ! write some scalars too
              self%MAXMESH
    close(fu)

  endsubroutine ! read

endmodule ! Main2Arrays_mod
