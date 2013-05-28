!> Datastructure that contains various arrays that could not factored out (yet).
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

  public :: createMain2Arrays
  public :: destroyMain2Arrays
  private :: createMain2ArraysImpl

  public :: Main2Arrays

  type Main2Arrays
    double precision , dimension(3,3)  :: bravais
    integer , dimension(48)  :: isymindex
    double complex , allocatable, dimension(:,:,:)  :: DSYMLL  !< tau symmetry matrices
    double precision , allocatable, dimension(:,:)  :: RBASIS  !< basis atom positions
    double precision , allocatable, dimension(:,:,:)  :: BZKP  !< kpoints for each mesh
    double precision , allocatable, dimension(:,:)  :: VOLCUB  !< kpoint weights
    double precision , allocatable, dimension(:)  :: VOLBZ     !< BZ volume?
    double precision , allocatable, dimension(:,:)  :: RR      !< lattice vectors
    integer , allocatable, dimension(:)  :: KMESH !< mapping E-point to k-mesh
    integer , allocatable, dimension(:)  :: NOFKS !< number of k points for each mesh
    integer , allocatable, dimension(:,:)  :: EZOA
    integer , allocatable, dimension(:)  :: NUMN0
    integer , allocatable, dimension(:,:)  :: INDN0
    double precision , allocatable, dimension(:)  :: ZAT  !< atomic numbers
    double precision , allocatable, dimension(:,:)  :: RCLS
    double precision , allocatable, dimension(:)  :: RMTREF
    integer , allocatable, dimension(:,:)  :: ATOM
    integer , allocatable, dimension(:)  :: CLS
    integer , allocatable, dimension(:)  :: NACLS

    double precision :: VREF !< repulsive screening pot. strength

    integer  :: NCLS
    integer  :: NSYMAT
    integer  :: MAXMESH
    integer  :: NR

    integer :: lmaxd
    integer :: iemxd
    integer :: nspind
    integer :: LMMAXD
    integer :: NAEZ
    integer :: LMXSPD
    integer :: KPOIBZ
    integer :: MAXMSHD
    integer :: nrd
    integer :: NACLSD
    integer :: nguessd !< not used?
    integer :: ekmd   !< not used, invalid
    integer :: smpid  !< not used
    integer :: lpot
    integer :: IRMD
    integer :: LMPOTD
  end type Main2Arrays

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object.
  !> @param[in,out] self    The Main2Arrays object to construct.
  !> @param[in]    dims    Dimension parameters
  subroutine createMain2Arrays(self, dims)
    use DimParams_mod
    implicit none
    type (Main2Arrays), intent(inout) :: self
    type (DimParams), intent(in) :: dims

    call createMain2ArraysImpl(self, dims%lmaxd,dims%iemxd,dims%nspind, &
    dims%LMMAXD,dims%NAEZ,dims%LMXSPD,dims%KPOIBZ,dims%MAXMSHD,dims%nrd, &
    dims%NACLSD,dims%nguessd,dims%ekmd, &
    dims%smpid,dims%lpot,dims%IRMD,dims%LMPOTD)

  end subroutine



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
    implicit none
    type (Main2Arrays), intent(inout) :: self
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

    self%NCLS = 0
    self%NSYMAT = 0
    self%MAXMESH = 0
    self%NR = 0


!    Repulsive reference potential
!    in future: move as parameter to inputfile
    self%VREF = 8.0d0

    self%lmaxd = lmaxd
    self%iemxd = iemxd
    self%nspind = nspind
    self%LMMAXD = LMMAXD
    self%NAEZ = NAEZ
    self%LMXSPD = LMXSPD
    self%IEMXD = IEMXD
    self%KPOIBZ = KPOIBZ
    self%MAXMSHD = MAXMSHD
    self%nrd = nrd
    self%NACLSD = NACLSD
    self%NSPIND = NSPIND
    self%nguessd = nguessd
    self%lmmaxd = lmmaxd
    self%ekmd = ekmd
    self%smpid = smpid
    self%lpot = lpot
    self%IRMD = IRMD
    self%LMPOTD = LMPOTD

    ALLOCATECHECK(self%DSYMLL(LMMAXD,LMMAXD,48))
    ALLOCATECHECK(self%RBASIS(3,NAEZ))
    ALLOCATECHECK(self%BZKP(3,KPOIBZ,MAXMSHD))
    ALLOCATECHECK(self%VOLCUB(KPOIBZ,MAXMSHD))
    ALLOCATECHECK(self%VOLBZ(MAXMSHD))
    ALLOCATECHECK(self%RR(3,0:NRD))
    ALLOCATECHECK(self%KMESH(IEMXD))
    ALLOCATECHECK(self%NOFKS(MAXMSHD))
    ALLOCATECHECK(self%EZOA(NACLSD,NAEZ))
    ALLOCATECHECK(self%NUMN0(NAEZ))
    ALLOCATECHECK(self%INDN0(NAEZ,NACLSD))
    ALLOCATECHECK(self%ZAT(NAEZ))
    ALLOCATECHECK(self%RCLS(3,NACLSD))
    ALLOCATECHECK(self%RMTREF(NAEZ))
    ALLOCATECHECK(self%ATOM(NACLSD,NAEZ))
    ALLOCATECHECK(self%CLS(NAEZ))
    ALLOCATECHECK(self%NACLS(1))

  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a Main2Arrays object.
  !> @param[in,out] self    The Main2Arrays object to destroy.
  subroutine destroyMain2Arrays(self)
    implicit none
    type (Main2Arrays), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%DSYMLL)
    DEALLOCATECHECK(self%RBASIS)
    DEALLOCATECHECK(self%BZKP)
    DEALLOCATECHECK(self%VOLCUB)
    DEALLOCATECHECK(self%VOLBZ)
    DEALLOCATECHECK(self%RR)
    DEALLOCATECHECK(self%KMESH)
    DEALLOCATECHECK(self%NOFKS)
    DEALLOCATECHECK(self%EZOA)
    DEALLOCATECHECK(self%NUMN0)
    DEALLOCATECHECK(self%INDN0)
    DEALLOCATECHECK(self%ZAT)
    DEALLOCATECHECK(self%RCLS)
    DEALLOCATECHECK(self%RMTREF)
    DEALLOCATECHECK(self%ATOM)
    DEALLOCATECHECK(self%CLS)
    DEALLOCATECHECK(self%NACLS)
  end subroutine

  !-----------------------------------------------------------------------------
  !> Writes Main2Arrays data to file.
  !> @param[in] self    The Main2Arrays object to write.
  subroutine writeMain2Arrays(self, filename)
    implicit none
    type (Main2Arrays), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: FILEHANDLE = 67

    open (FILEHANDLE, file=filename, form='unformatted')
    write (FILEHANDLE) self%BRAVAIS, &
                       self%ISYMINDEX, &
                       self%DSYMLL, &
                       self%RBASIS, &
                       self%BZKP, &
                       self%VOLCUB, &
                       self%VOLBZ, &
                       self%RR, &
                       self%KMESH, &
                       self%NOFKS, &
                       self%EZOA, &
                       self%NUMN0, &
                       self%INDN0, &
                       self%ZAT, &
                       self%RCLS, &
                       self%RMTREF, &
                       self%VREF, &
                       self%ATOM, &
                       self%CLS, &
                       self%NACLS, &
                       self%NCLS, &  ! write some scalars too
                       self%NSYMAT, &
                       self%MAXMESH, &
                       self%NR
    close (FILEHANDLE)

  end subroutine

  !-----------------------------------------------------------------------------
  !> Reads Main2Arrays data from file.
  !> @param[in,out] self    The Main2Arrays object to read.
  subroutine readMain2Arrays(self, filename)
    implicit none
    type (Main2Arrays), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: FILEHANDLE = 67

    open (FILEHANDLE, file=filename, form='unformatted')
    read  (FILEHANDLE) self%BRAVAIS, &
                       self%ISYMINDEX, &
                       self%DSYMLL, &
                       self%RBASIS, &
                       self%BZKP, &
                       self%VOLCUB, &
                       self%VOLBZ, &
                       self%RR, &
                       self%KMESH, &
                       self%NOFKS, &
                       self%EZOA, &
                       self%NUMN0, &
                       self%INDN0, &
                       self%ZAT, &
                       self%RCLS, &
                       self%RMTREF, &
                       self%VREF, &
                       self%ATOM, &
                       self%CLS, &
                       self%NACLS, &
                       self%NCLS, &  ! write some scalars too
                       self%NSYMAT, &
                       self%MAXMESH, &
                       self%NR
    close (FILEHANDLE)

  end subroutine

end module
