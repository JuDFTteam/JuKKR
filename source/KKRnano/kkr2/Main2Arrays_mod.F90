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
    double complex , allocatable, dimension(:,:,:)  :: DSYMLL
    double precision , allocatable, dimension(:,:)  :: RBASIS
    double precision , allocatable, dimension(:,:,:)  :: BZKP
    double precision , allocatable, dimension(:,:)  :: VOLCUB
    double precision , allocatable, dimension(:)  :: VOLBZ
    double precision , allocatable, dimension(:,:)  :: RR
    integer , allocatable, dimension(:)  :: KMESH
    integer , allocatable, dimension(:)  :: NOFKS
    integer , allocatable, dimension(:,:)  :: EZOA
    integer , allocatable, dimension(:)  :: NUMN0
    integer , allocatable, dimension(:,:)  :: INDN0
    complex , allocatable, dimension(:,:,:)  :: PRSC ! move to KKRresults?
    double precision , allocatable, dimension(:)  :: ZAT
    double precision , allocatable, dimension(:,:,:)  :: RCLS
    double precision , allocatable, dimension(:)  :: RMTREF
    double precision , allocatable, dimension(:)  :: VREF
    integer , allocatable, dimension(:,:)  :: ATOM
    integer , allocatable, dimension(:)  :: CLS
    integer , allocatable, dimension(:)  :: NACLS
    integer , allocatable, dimension(:)  :: REFPOT

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
    integer :: NREFD
    integer :: NCLSD
    integer :: nguessd
    integer :: ekmd
    integer :: smpid
    integer :: lpot
    integer :: IRMD
    integer :: LMPOTD
  end type Main2Arrays

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object.
  !> @param[inout] self    The Main2Arrays object to construct.
  !> @param[in]    dims    Dimension parameters
  subroutine createMain2Arrays(self, dims)
    use DimParams_mod
    implicit none
    type (Main2Arrays), intent(inout) :: self
    type (DimParams), intent(in) :: dims

    call createMain2ArraysImpl(self, dims%lmaxd,dims%iemxd,dims%nspind, &
    dims%LMMAXD,dims%NAEZ,dims%LMXSPD,dims%KPOIBZ,dims%MAXMSHD,dims%nrd, &
    dims%NACLSD,dims%NREFD,dims%NCLSD,dims%nguessd,dims%ekmd, &
    dims%smpid,dims%lpot,dims%IRMD,dims%LMPOTD)

  end subroutine



  !-----------------------------------------------------------------------------
  !> Constructs a Main2Arrays object. (implementation, don't call directly!)
  !> @param[inout] self    The Main2Arrays object to construct.
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
  !> @param[in]    NREFD
  !> @param[in]    NCLSD
  !> @param[in]    nguessd
  !> @param[in]    lmmaxd
  !> @param[in]    ekmd
  !> @param[in]    smpid
  !> @param[in]    lpot
  !> @param[in]    IRMD
  !> @param[in]    LMPOTD
  subroutine createMain2ArraysImpl(self, lmaxd,iemxd,nspind,LMMAXD,NAEZ,LMXSPD,KPOIBZ,MAXMSHD,nrd,NACLSD,NREFD,NCLSD,nguessd,ekmd,smpid,lpot,IRMD,LMPOTD)
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
    integer, intent(in) ::  NREFD
    integer, intent(in) ::  NCLSD
    integer, intent(in) ::  nguessd
    integer, intent(in) ::  ekmd
    integer, intent(in) ::  smpid
    integer, intent(in) ::  lpot
    integer, intent(in) ::  IRMD
    integer, intent(in) ::  LMPOTD

    integer :: memory_stat

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
    self%NREFD = NREFD
    self%NCLSD = NCLSD
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
    ALLOCATECHECK(self%PRSC(NGUESSD*LMMAXD,EKMD,NSPIND-SMPID+1))
    ALLOCATECHECK(self%ZAT(NAEZ))
    ALLOCATECHECK(self%RCLS(3,NACLSD,NCLSD))
    ALLOCATECHECK(self%RMTREF(NREFD))
    ALLOCATECHECK(self%VREF(NAEZ))
    ALLOCATECHECK(self%ATOM(NACLSD,NAEZ))
    ALLOCATECHECK(self%CLS(NAEZ))
    ALLOCATECHECK(self%NACLS(NCLSD))
    ALLOCATECHECK(self%REFPOT(NAEZ))

  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a Main2Arrays object.
  !> @param[inout] self    The Main2Arrays object to destroy.
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
    DEALLOCATECHECK(self%PRSC)
    DEALLOCATECHECK(self%ZAT)
    DEALLOCATECHECK(self%RCLS)
    DEALLOCATECHECK(self%RMTREF)
    DEALLOCATECHECK(self%VREF)
    DEALLOCATECHECK(self%ATOM)
    DEALLOCATECHECK(self%CLS)
    DEALLOCATECHECK(self%NACLS)
    DEALLOCATECHECK(self%REFPOT)
  end subroutine

end module
