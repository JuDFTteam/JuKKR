  program MAIN0

! Explanation of most variables follows below

!     ALAT                     : lattice constant (in a.u.)
!     ABASIS,BBASIS,CBASIS,    : scaling factors for rbasis
!     E1,E2,                   : energies needed in EMESHT
!     HFIELD                   : external magnetic field, for
!                              : initial potential shift in
!                              : spin polarised case
!     TK,                      : temperature
!     VCONST,                  : potential shift
!     BRAVAIS(3,3),            : bravais lattice vectors
!     RECBV(3,3),              : reciprocal basis vectors
!     CLEB(NCLEB,2),           : GAUNT coefficients (GAUNT)
!     RMTREF(NREFD),           : muffin-tin radius of reference system
!     RBASIS(3,NAEZD),         : position of atoms in the unit cell
!                              : in units of bravais vectors
!     RCLS(3,NACLSD,NCLSD),    : real space position of atom in cluster
!     RR(3,0:NRD)              : set of real space vectors (in a.u.)
!     VBC(2),                  : potential constants
!     WG(LASSLD),              : integr. weights for Legendre polynomials
!     YRG(LASSLD,0:LASSLD)     : spherical harmonics (GAUNT2)
!     ZAT(NAEZD)               : nuclear charge
!     INTERVX,INTERVY,INTERVZ, : number of intervals in x,y,z-direction
!                              : for k-net in IB of the BZ
!     ICST,                    : number of Born approximation
!     IEND,                    : number of nonzero gaunt coeffizients
!     IFILE,                   : unit specifier for potential card
!     IPE,IPF,IPFE,            : not real used, IPFE should be 0
!     KHFELD,                  : 0,1: no / yes external magnetic field
!     KVREL,                   : 0,1 : non / scalar relat. calculation
!     LMAX,                    : maximum l component in
!                              : wave function expansion
!     LPOT,                    : maximum l component in
!                              : potential expansion
!     NAEZ,                    : number of atoms in unit cell
!     NCLS,                    : number of reference clusters
!     NPNT1,NPNT2,NPNT3,       : number of E points (EMESHT)
!     NPOL,                    : number of Matsubara Pols (EMESHT)
!     NR,                      : number of real space vectors rr
!     NREF,                    : number of diff. ref. potentials
!     NSPIN,                   : counter for spin directions
!     IGUESS                   : 0,1 : no / yes (sc) initial guess, set
!                              : IGUESSD to 1 in inc.p if needed
!     BCP                      : 0,1 : no / yes bc-preconditioning, set
!                              : BCPD to 1 in inc.p if needed
!     QBOUND                   : exit condition for self-consistent
!                              : iteration
!     QMRBOUND                 : exit condition for QMR iterations
!     SCFSTEPS                 : number of scf iterations
!     INIPOL(NAEZD),           : initial spin polarisation

!     ATOM(NACLSD,NAEZD),      : atom at site in cluster
!     CLS(NAEZD),              : cluster around atom
!     NACLS(NCLSD),            : number of atoms in cluster
!     EZOA(NACLSD,NAEZD),      : EZ of atom at site in cluster
!     ICLEB(NCLEB,3),          : pointer array
!     RMT(NAEZD)               : muffin-tin radius of true system
!     RMTNEW(NAEZD)            : adapted muffin-tin radius
!     RWS(NAEZD)               : Wigner Seitz radius
!     ICST                     : the regular non spherical wavefunctions, the
!                              : alpha matrix and the t-matrix in the ICST-th. born approx
!     IMT(NAEZD),              : r point at MT radius
!     IPAN(NAEZD),             : number of panels in non-MT-region
!     IRC(NAEZD),              : r point for potential cutting
!     IRCUT(0:IPAND,NAEZD),    : r points of panel borders
!     IRMIN(NAEZD),            : max r for spherical treatment
!     IRNS(NAEZD)              : number r points for non spher. treatm.
!     IRWS(NAEZD),             : r point at WS radius
!     LMSP(NAEZD,LMXSPD)       : 0,1 : non/-vanishing lm=(l,m) component
!                              : of non-spherical potential
!     LLMSP(NAEZD,NFUND)       : lm=(l,m) of 'nfund'th nonvanishing
!                              : component of non-spherical pot.
!     JEND(LMPOTD,             : pointer array for icleb()
!     LOFLM(LM2D),             : l of lm=(l,m) (GAUNT)
!     NTCELL(NAEZD),           : index for WS cell
!     REFPOT(NAEZD)            : ref. pot. card  at position

!     A(NAEZD),B(NAEZD)        : contants for exponential r mesh
!     R(IRMD,NAEZD)            : radial mesh ( in units a Bohr)
!     DRDI(IRMD,NAEZD)         : derivative dr/di
!     THETAS(IRID,NFUND,NCELLD): shape function
!                              :         ( 0 outer space
!                              : THETA = (
!                              :         ( 1 inside WS cell
!                              : in spherical harmonics expansion
!     LCORE(20,NPOTD)          : angular momentum of core states
!     NCORE(NPOTD)             : number of core states
! ----------------------------------------------------------------------

    use Config_Reader

!     the inc.p wrapper module is used to include the inc.p and inc.cls
!     parameters - this works also for free form source files
    !use inc_p_wrapper_module
    implicit none

!   double precision KIND constant
    integer, parameter :: DP = kind(1.0D0)

!     .. Parameters ..

    integer :: NSYMAXD
    integer :: MAXMSHD

!     Maximal number of Brillouin zone symmetries, 48 is largest
!     possible number
    parameter (NSYMAXD=48)
!     Maximal number of k-meshes used
    parameter (MAXMSHD=8)

!     inc.p, inc.cls parameters

!     NAEZD
!     LPOTD
!     LMAXD
!     NREFD
!     IRID
!     BCPD
!     NACLSD
!     NCLEB
!     IRMD
!     IEMXD
!     NGSHD
!     IGUESSD
!     IPAND
!     ISHLD
!     IRNSD
!     KPOIBZ
!     KREL
!     NFUND
!     NATRCD
!     NCLSD
!     NMAXD
!     NRD
!     NSPIND
!     NUTRCD
!     NXIJD


!     ..
!     .. Local Scalars ..

    double precision :: RMAX
    double precision :: GMAX
    double precision :: RCUTJIJ
    double precision :: RCUTTRC

!     .. KKR calculation options ..
    double precision :: VCONST

    integer :: ICST
    integer :: IRM    !     eq IRMD ?
    integer :: ISHIFT
    integer :: KHFELD
    integer :: KPRE
    integer :: KTE
    integer :: KVMAD
    integer :: KVREL
    integer :: KXC
    integer :: LMAX
    integer :: NSPIN
    integer :: KFORCE

    logical :: JIJ
    logical :: LDAU

!     .. KKR calculation options options derived from others ..
    integer :: LMPOT
    integer :: LPOT
    integer :: NSRA

!     ..
!     .. Local Arrays ..
    double precision :: VBC(2)

!     .. Energy Mesh ..
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double precision :: TK

    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    integer :: IELAST
    integer :: MAXMESH

!    complex(kind=DP) :: EZ(IEMXD)
!    complex(kind=DP) :: WEZ(IEMXD)
!    integer :: KMESH(IEMXD)
    complex(kind=DP), dimension(:), allocatable :: EZ
    complex(kind=DP), dimension(:), allocatable :: WEZ
    integer, dimension(:), allocatable :: KMESH

!     .. Lattice ..
    integer :: NAEZ
    integer :: NR

    double precision :: ALAT
    double precision :: VOLUME0

    double precision :: BRAVAIS(3,3)
    double precision :: RECBV(3,3)

!    double precision :: RBASIS(3,NAEZD)
!    double precision :: RMT(NAEZD)
!    double precision :: RMTREF(NREFD)
!    double precision :: RR(3,0:NRD)
!    double precision :: RWS(NAEZD)
!    double precision :: ZAT(NAEZD)
    double precision, dimension(:,:), allocatable :: RBASIS
    double precision, dimension(:),   allocatable :: RMT
    double precision, dimension(:),   allocatable :: RMTREF
    double precision, dimension(:,:), allocatable :: RR
    double precision, dimension(:),   allocatable :: RWS
    double precision, dimension(:),   allocatable :: ZAT

!     .. Lattice aux. ..
    double precision :: ABASIS
    double precision :: BBASIS
    double precision :: CBASIS


!     .. shape functions ..
!    double precision :: THETAS(IRID,NFUND,NCELLD)
!    double precision :: GSH(NGSHD)
!    integer :: IFUNM(LMXSPD,NAEZD)
!    integer :: IPAN(NAEZD)
!    integer :: LLMSP(NFUND,NAEZD)
!    integer :: LMSP(LMXSPD,NAEZD)
!    integer :: ILM(NGSHD,3)
!    integer :: NFU(NAEZD)
!    integer :: NTCELL(NAEZD)
!    integer :: IMAXSH(0:LMPOTD)
    double precision, dimension(:,:,:), allocatable :: THETAS
    double precision, dimension(:),     allocatable :: GSH
    integer, dimension(:,:), allocatable :: IFUNM
    integer, dimension(:),   allocatable :: IPAN
    integer, dimension(:,:), allocatable :: LLMSP
    integer, dimension(:,:), allocatable :: LMSP
    integer, dimension(:,:), allocatable :: ILM
    integer, dimension(:),   allocatable :: NFU
    integer, dimension(:),   allocatable :: NTCELL
    integer, dimension(:),   allocatable :: IMAXSH

!     .. reference clusters
    double precision :: RCUTXY
    double precision :: RCUTZ

    integer :: NCLS
    integer :: NREF

!    double precision :: RCLS(3,NACLSD,NCLSD)
!    integer :: ATOM(NACLSD,NAEZD)
!    integer :: CLS(NAEZD)
!    integer :: EZOA(NACLSD,NAEZD)
!    integer :: NACLS(NCLSD)
!    integer :: NUMN0(NAEZD)
!    integer :: INDN0(NAEZD,NACLSD)
!    integer :: REFPOT(NAEZD)
    double precision, dimension(:,:,:), allocatable :: RCLS

    integer, dimension(:,:), allocatable :: ATOM
    integer, dimension(:),   allocatable :: CLS
    integer, dimension(:,:), allocatable :: EZOA
    integer, dimension(:),   allocatable :: NACLS
!     NUMN0(i) gives number of ref. cluster atoms around atom i
    integer, dimension(:),   allocatable :: NUMN0
!     INDN0(i,j) gives the atom index (from basis) of cluster atom j
!     around central atom i
    integer, dimension(:,:), allocatable :: INDN0
    integer, dimension(:),   allocatable :: REFPOT

!     .. reference system ..
!    double precision :: VREF(NAEZD)
    double precision, dimension(:), allocatable :: VREF

!     .. radial mesh(es) ..
!    double precision :: A(NAEZD)
!    double precision :: B(NAEZD)
!    double precision :: DRDI(IRMD,NAEZD)
!    double precision :: R(IRMD,NAEZD)
!    integer :: IMT(NAEZD)
!    integer :: IRC(NAEZD)
!    integer :: IRCUT(0:IPAND,NAEZD)
!    integer :: IRMIN(NAEZD)
!    integer :: IRNS(NAEZD)
!    integer :: IRWS(NAEZD)
    double precision, dimension(:),   allocatable  :: A
    double precision, dimension(:),   allocatable  :: B
    double precision, dimension(:,:), allocatable  :: DRDI
    double precision, dimension(:,:), allocatable  :: R
    integer, dimension(:),  allocatable  :: IMT
    integer, dimension(:),  allocatable  :: IRC
    integer, dimension(:,:),allocatable  :: IRCUT
    integer, dimension(:),  allocatable  :: IRMIN
    integer, dimension(:),  allocatable  :: IRNS
    integer, dimension(:),  allocatable  :: IRWS

!     .. Brillouin zone ..
!     kpoints in each direction
    integer :: INTERVX
    integer :: INTERVY
    integer :: INTERVZ
!     number of symmetries = number of symmetry matrices
    integer :: NSYMAT

!         symmetry matrices
!    complex(kind=DP) DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
!    integer :: ISYMINDEX(NSYMAXD)
    complex(kind=DP), dimension(:,:,:), allocatable :: DSYMLL
    integer, dimension(:), allocatable :: ISYMINDEX

!     .. Spherical harmonics ..
!    double precision :: WG(LASSLD) ! not passed to kkr2
!    double precision :: YRG(LASSLD,0:LASSLD,0:LASSLD) ! not passed to kkr2
    double precision, dimension(:),     allocatable :: WG ! not passed to kkr2
    double precision, dimension(:,:,:), allocatable :: YRG ! not passed to kkr2

!     .. Clebsch-Gordon coefficients
    integer :: IEND

!    double precision :: CLEB(NCLEB,2)
!    integer :: ICLEB(NCLEB,3)
!    integer :: JEND(LMPOTD,0:LMAXD,0:LMAXD)
!    integer :: LOFLM(LM2D) ! gives l from LM index

    double precision, dimension(:,:), allocatable :: CLEB
    integer, dimension(:,:),          allocatable :: ICLEB
    integer, dimension(:,:,:),        allocatable :: JEND
    integer, dimension(:),            allocatable :: LOFLM ! gives l from LM index


!     .. core states ..
!    integer :: ITITLE(20,NPOTD)
!    integer :: LCORE(20,NPOTD)
!    integer :: NCORE(NPOTD)
    integer, dimension(:,:), allocatable :: ITITLE
    integer, dimension(:,:), allocatable :: LCORE
    integer, dimension(:),   allocatable :: NCORE


!     .. Self-consistency parameters ..
    double precision :: FCM
    double precision :: MIXING
    double precision :: QBOUND

    integer :: SCFSTEPS
    integer :: IMIX
!     has no effect, use ITDBRYD in inc.p instead
    integer :: ITDBRY

!     .. Iterative solver options ..
    double precision :: QMRBOUND

    integer :: IGUESS
    integer :: BCP

!     .. auxillary variables, not passed to kkr2
    double precision :: PI
    double precision :: EREF
    double precision :: HFIELD
    double precision :: E2IN

    integer :: I1
    integer :: IPE
    integer :: IPF
    integer :: IPFE
    integer :: IFILE
    integer :: IE

    ! error code for allocations
    integer :: ierror

    logical :: EVREF
    logical :: LCARTESIAN

    character(len=40) POTENTIAL_FILENAME
    character(len=40) SHAPEFUN_FILENAME


!    complex(kind=DP) DEZ(IEMXD)
!    double precision :: RMTNEW(NAEZD)
!    integer :: INIPOL(NAEZD)

    complex(kind=DP), dimension(:), allocatable :: DEZ ! needed for EMESHT
    double precision, dimension(:), allocatable :: RMTNEW
    integer, dimension(:), allocatable :: INIPOL

    integer ::   LASSLD,LMMAXD,LMPOTD,LMXSPD,LM2D,NPOTD

! ------------- parameters derived from others or calculated
    integer ::   LPOTD
    integer ::   IEMXD
    integer ::   NCLEB
    integer, parameter :: KREL = 0

! ------------- parameters from global.conf (former inc.p, inc.cls)

    integer :: LMAXD
    integer :: NSPIND
    integer :: NAEZD
    integer :: IRNSD
    integer :: TRC
    integer :: IRMD
    integer :: NREFD
    integer :: NRD
    integer :: IRID
    integer :: NFUND
    integer :: NCELLD
    integer :: NGSHD
    integer :: NACLSD
    integer :: NCLSD
    integer :: IPAND
    integer :: NXIJD
    integer :: NATRCD
    integer :: NUTRCD
    integer :: KPOIBZ
    integer :: EKMD
    integer :: IGUESSD
    integer :: BCPD
    integer :: NMAXD
    integer :: ISHLD
    integer :: LLY

    integer :: SMPID
    integer :: EMPID
    integer :: NTHRDS
    integer :: XDIM
    integer :: YDIM
    integer :: ZDIM
    integer :: NATBLD
    integer :: ITDBRYD

    character(len=40) :: variable
    integer :: next_ptr

! ------------ end of declarations ---------------------------------

    type (ConfigReader) :: conf

    call createConfigReader(conf)
    call parseFile(conf, 'global.conf', ierror)
    if (ierror /= 0) stop

    call getValueInteger(conf, "LMAXD", LMAXD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NSPIND", NSPIND, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NAEZD", NAEZD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IRNSD", IRNSD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "TRC", TRC, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IRMD", IRMD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NREFD", NREFD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NRD", NRD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IRID", IRID, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NFUND", NFUND, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NCELLD", NCELLD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NGSHD", NGSHD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NACLSD", NACLSD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NCLSD", NCLSD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IPAND", IPAND, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NXIJD", NXIJD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NATRCD", NATRCD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NUTRCD", NUTRCD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "KPOIBZ", KPOIBZ, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IGUESSD", IGUESSD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "BCPD", BCPD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NMAXD", NMAXD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "ISHLD", ISHLD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "LLY", LLY, ierror)
    if (ierror /= 0) stop

    call getValueInteger(conf, "SMPID", SMPID, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "EMPID", EMPID, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NTHRDS", NTHRDS, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "XDIM", XDIM, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "YDIM", YDIM, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "ZDIM", ZDIM, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NATBLD", NATBLD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "ITDBRYD", ITDBRYD, ierror)
    if (ierror /= 0) stop

    write(*,*) "The following variables have not been read:"
    next_ptr = 1
    do
      call getUnreadVariable(conf, variable, next_ptr, ierror)
      if (ierror /= 0) exit
      write (*,*) variable
    end do

    call destroyConfigReader(conf)

    LPOTD = 2*LMAXD
    LMMAXD= (LMAXD+1)**2
    NPOTD= NSPIND*NAEZD
    LMPOTD= (LPOTD+1)**2
    LMXSPD= (2*LPOTD+1)**2
    LASSLD=4*LMAXD
    LM2D= (2*LMAXD+1)**2

    PI = 4.0D0*ATAN(1.0D0)
    EFERMI = 0.0d0

!-----------------------------------------------------------------------------
! Array allocations BEGIN 1
!-----------------------------------------------------------------------------

    ! Lattice
    allocate(RBASIS(3,NAEZD), stat=ierror)
    allocate(RMT(NAEZD), stat=ierror)
    allocate(RMTREF(NREFD), stat=ierror)
    allocate(RR(3,0:NRD), stat=ierror)
    allocate(RWS(NAEZD), stat=ierror)
    allocate(ZAT(NAEZD), stat=ierror)

    !Shape functions
    allocate(THETAS(IRID,NFUND,NCELLD), stat=ierror)
    allocate(GSH(NGSHD), stat=ierror)
    allocate(IFUNM(LMXSPD,NAEZD), stat=ierror)
    allocate(IPAN(NAEZD), stat=ierror)
    allocate(LLMSP(NFUND,NAEZD), stat=ierror)
    allocate(LMSP(LMXSPD,NAEZD), stat=ierror)
    allocate(ILM(NGSHD,3), stat=ierror)
    allocate(NFU(NAEZD), stat=ierror)
    allocate(NTCELL(NAEZD), stat=ierror)
    allocate(IMAXSH(0:LMPOTD), stat=ierror)

    ! Reference Cluster
    allocate(RCLS(3,NACLSD,NCLSD), stat=ierror)
    allocate(ATOM(NACLSD,NAEZD), stat=ierror)
    allocate(CLS(NAEZD), stat=ierror)
    allocate(EZOA(NACLSD,NAEZD), stat=ierror)
    allocate(NACLS(NCLSD), stat=ierror)
    allocate(NUMN0(NAEZD), stat=ierror)
    allocate(INDN0(NAEZD,NACLSD), stat=ierror)
    allocate(REFPOT(NAEZD), stat=ierror)

    ! Reference system
    allocate(VREF(NAEZD), stat=ierror)

    ! Radial mesh(es)
    allocate(A(NAEZD), stat=ierror)
    allocate(B(NAEZD), stat=ierror)
    allocate(DRDI(IRMD,NAEZD), stat=ierror)
    allocate(R(IRMD,NAEZD), stat=ierror)
    allocate(IMT(NAEZD), stat=ierror)
    allocate(IRC(NAEZD), stat=ierror)
    allocate(IRCUT(0:IPAND,NAEZD), stat=ierror)
    allocate(IRMIN(NAEZD), stat=ierror)
    allocate(IRNS(NAEZD), stat=ierror)
    allocate(IRWS(NAEZD), stat=ierror)

    !         symmetry matrices
    allocate(DSYMLL(LMMAXD,LMMAXD,NSYMAXD), stat=ierror)
    allocate(ISYMINDEX(NSYMAXD), stat=ierror)

    ! spherical harmonics
    allocate(WG(LASSLD), stat=ierror)
    allocate(YRG(LASSLD,0:LASSLD,0:LASSLD), stat=ierror)

    !     .. core states ..
    allocate(ITITLE(20,NPOTD), stat=ierror)
    allocate(LCORE(20,NPOTD), stat=ierror)
    allocate(NCORE(NPOTD), stat=ierror)

    ! auxillary
    allocate(RMTNEW(NAEZD), stat=ierror)
    allocate(INIPOL(NAEZD), stat=ierror)

!-----------------------------------------------------------------------------
! Array allocations END 1
!-----------------------------------------------------------------------------

    call RINPUT99(BRAVAIS,ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS, &
                  E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3, &
                  SCFSTEPS,IMIX,MIXING,QBOUND,FCM, &
                  ITDBRY,IRNS,NTCELL,NAEZ,IRM,ZAT, &
                  NREF,ICST,IFILE,IPE,IPF,IPFE, &
                  KHFELD,KPRE,KTE, &
                  KVMAD,KVREL,KXC,LMAX,LMPOT,LPOT, &
                  NSPIN,REFPOT, &
                  ISHIFT,INTERVX,INTERVY,INTERVZ, &
                  HFIELD,VBC,VCONST,INIPOL, &
                  POTENTIAL_FILENAME,SHAPEFUN_FILENAME, &
                  RCUTZ,RCUTXY,RCUTJIJ,JIJ,RCUTTRC, &
                  LDAU, &
                  RMTREF,KFORCE, &
                  IGUESS,BCP,QMRBOUND,LCARTESIAN,RMAX,GMAX, &
                  LMAXD, IRNSD, TRC, LPOTD, NSPIND, &
                  IRMD, NAEZD)

! unnecessary parameters, read in for compatibility: IRM, KHFELD, ...

! --------------- calculated parameters --------------------------------------
    IEMXD = NPOL + NPNT1 + NPNT2 + NPNT3
    NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2

!-----------------------------------------------------------------------------
! Array allocations BEGIN 2
!-----------------------------------------------------------------------------

    ! Energy mesh
    allocate(EZ(IEMXD), stat=ierror)
    allocate(WEZ(IEMXD), stat=ierror)
    allocate(KMESH(IEMXD), stat=ierror)

    !   auxillary
    allocate(DEZ(IEMXD), stat=ierror)

    ! Clebsch-Gordon coefficients
    allocate(CLEB(NCLEB,2), stat=ierror)
    allocate(ICLEB(NCLEB,3), stat=ierror)
    allocate(JEND(LMPOTD,0:LMAXD,0:LMAXD), stat=ierror)
    allocate(LOFLM(LM2D), stat=ierror)


!-----------------------------------------------------------------------------
! Array allocations END 2
!-----------------------------------------------------------------------------



!     in case of a LDA+U calculation - read file 'ldauinfo'
!     and write 'wldau.unf', if it does not exist already
    if (LDAU) then
      call ldauinfo_read(LMAXD, NSPIND, ZAT, NAEZD)
    end if

!===================================================================
! next short section used to search for file 'VREF'
! if present - DP-value given in that file will be used as EREF
!===================================================================
    inquire(file='VREF',exist=EVREF)
    if (EVREF) then
      open(87,file='VREF',form='formatted')
      read(87,*) EREF
      close(87)
    endif

    !     The default value for the repulsive reference potential is 8
    !     This value is used if file VREF does not exist
    !     (VREF should contain only one double precision value used as EREF)
    !     in future: move as parameter to inputfile
    do I1 = 1,NAEZD
      VREF(I1) = 8.D0
      if (EVREF) VREF(I1) = EREF
    end do
!===================================================================
!===================================================================


!     IF(NAEZD.GT.100) IPE=0


    NSRA = 1
    if (KVREL >= 1) NSRA = 2

    call TESTDIM(NSPIN,NAEZ,LMAX,IRM,NREF,IRNS, &
    LMAXD,IRMD,IRNSD,NREFD,NSPIND)

    open (19,file=SHAPEFUN_FILENAME,status='old',form='formatted')
    open (IFILE,file=POTENTIAL_FILENAME,status='old',form='formatted')

    E2IN = E2

    ! read starting potential and shapefunctions
    call STARTB1(IFILE,IPF,IPFE,IPE,KHFELD, &
                 1,NAEZ, &
                 RMTNEW,RMT,ITITLE,HFIELD,IMT,IRC,VCONST, &
                 IRNS,LPOT,NSPIN,IRMIN,NTCELL,IRCUT,IPAN, &
                 THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN, &
                 VBC,RWS,LCORE,NCORE,DRDI, &
                 R,ZAT,A,B,IRWS,INIPOL,1, &
!                New after inc.p replace
                 IPAND, IRID, NFUND, IRMD, NCELLD, &
                 NAEZD, IRNSD)

    close(IFILE)
    close(19)

! ----------------------------------------------------------------------
! update Fermi energy, adjust energy window according to running options

    if ( NPOL == 0 ) EFERMI = E2IN
!     a test if E2IN is different from E2 after call to STARTB1,
!     before it was =E2
    if ( DABS(E2IN-E2) > 1D-10 .and. NPOL /= 0 ) E2 = E2IN

! --> set up energy contour

    call EMESHT(EZ,DEZ,IELAST,E1,E2,E2IN,TK, &
    NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*DEZ(IE)
    end do


    call GAUNT2(WG,YRG,LMAX)
    call GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND,NCLEB)
!      OPEN(56,FILE='gaunt',FORM='formatted')
!      WRITE(56,FMT='(I10)') IEND
!      DO I = 1,IEND
!      WRITE(56,FMT='(3I5,1P,D25.17)')
!     +            ICLEB(I,1),ICLEB(I,2),ICLEB(I,3),CLEB(I,1)
!      END DO
!      CLOSE(56)

! --> setup of GAUNT coefficients C(l,m;l',m';l'',m'') for all
!     nonvanishing (l'',m'')-components of the shape functions THETAS

    call SHAPE(LPOT,NAEZ,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG,LMAX,NGSHD)
!      OPEN(56,FILE='gaunt_shape',FORM='formatted')
!      WRITE(56,FMT='(I10)') IMAXSH(LMPOTD)
!      DO I = 1,IMAXSH(LMPOTD)
!      WRITE(56,FMT='(3I5,1P,D25.17)')
!     +            ILM(I,1),ILM(I,2),ILM(I,3),GSH(I)
!      END DO
!      CLOSE(56)

! ================================================ deal with the lattice
    call LATTIX99(ALAT,BRAVAIS, &
                  RECBV,VOLUME0,RR,NR,NRD)

    call SCALEVEC(RBASIS,ABASIS,BBASIS,CBASIS, &
                  NAEZ,BRAVAIS,LCARTESIAN)
! ======================================================================

! ======================================================================

!   initialise arrays for clsgen99 to garbage values
    EZOA = -1
    ATOM = -1
    NACLS = -1
    CLS = -1
    NUMN0 = -1
    INDN0 = -1

    call CLSGEN99(NAEZ,RR,NR,RBASIS,CLS,NACLS,REFPOT,ATOM, &
                  EZOA, &
                  RCLS,RCUTZ,RCUTXY, &
                  NUMN0,INDN0, &
                  NRD, NCLSD, NACLSD)

! xcpl test dimensions for Jij-calculation ..
    call CLSJIJ0(NAEZ,RR,NR,RBASIS,RCUTJIJ,JIJ,NRD,NXIJD)
    ! xcpl .

    ! C2 ===================================================================
    ! C2  generate cluster 2 in order to truncate GLLKE, GLLH, etc.
    ! C2 ===================================================================

    if (TRC == 1) then
      call CLSGEN_TRC(NAEZ,RR,NR,RBASIS, &
      RCUTTRC,ALAT, &
      NRD, NATRCD, NUTRCD)
    endif
! C2 ===================================================================
! C2 ===================================================================



! ======================================================================
!     setting up kpoints
! ======================================================================
    call BZKINT0(NAEZ, &
                 RBASIS,BRAVAIS,RECBV, &
                 NSYMAT,ISYMINDEX, &
                 DSYMLL, &
                 INTERVX,INTERVY,INTERVZ, &
                 IELAST,EZ,KMESH,MAXMESH,MAXMSHD, &
                 LMAX, IEMXD, KREL, KPOIBZ, EKMD)
    ! after return from bzkint0, EKMD contains the right value

! ======================================================================
! ======================================================================


! ======================================================================
! =        write out information for the other program parts           =
! ======================================================================

    if (IGUESS .ne. IGUESSD) then
      write(*,*) 'WARNING: IGUESS from inputcard ignored. Using IGUESSD = ', IGUESSD
      IGUESS = IGUESSD
    endif
    if (BCP .ne. BCPD) then
      write(*,*) 'WARNING: BCP from inputcard ignored. Using BCPD = ', BCPD
      BCP = BCPD
    endif

!     Conversion of RMAX and GMAX to units of ALAT
    RMAX = RMAX*ALAT
    GMAX = GMAX/ALAT

    call TESTDIMLAT(ALAT,BRAVAIS,RECBV,RMAX,GMAX, &
                    NMAXD, ISHLD)


    call write_dimension_parameters( &
        LMAXD, &
        NSPIND, &
        NAEZD, &
        IRNSD, &
        TRC, &
        IRMD, &
        NREFD, &
        NRD, &
        IRID, &
        NFUND, &
        NCELLD, &
        NGSHD, &
        NACLSD, &
        NCLSD, &
        IPAND, &
        NXIJD, &
        NATRCD, &
        NUTRCD, &
        KPOIBZ, &
        IGUESSD, &
        BCPD, &
        NMAXD, &
        ISHLD, &
        LLY, &
        SMPID, &
        EMPID, &
        NTHRDS, &
        XDIM, &
        YDIM, &
        ZDIM, &
        NATBLD, &
        ITDBRYD, &
        IEMXD, &
        EKMD)

    call WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ, &
    BRAVAIS,RECBV,VOLUME0,RMAX,GMAX, &
    EFERMI,VBC, &
    SCFSTEPS,LCORE,NCORE, &
    NSRA,NAEZ,NREF,NSPIN,LMAX, &
    NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI, &
    REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB, &
    ATOM,CLS,RCLS,NACLS,LOFLM, &
    RBASIS,RR,EZOA, &
    KMESH,MAXMESH,NSYMAT, &
    DSYMLL, &
    A,B,IFUNM,ITITLE, &
    LMSP,NTCELL,THETAS, &
    LPOT,LMPOT, &
    IMIX,MIXING,QBOUND,FCM,KPRE,KTE, &
    KVMAD,KXC,ISHIFT,KFORCE, &
    LLMSP,IMT,IRC,IRMIN,IRNS,IRWS,RWS,RMT,NFU, &
    GSH,ILM,IMAXSH, &
    IEMXD,IRMD,LMPOTD,NPOTD,NAEZD, &
    LMMAXD,IPAND,NREFD,LMAXD, &
    NCLEB,NACLSD,NCLSD,LM2D,NRD, &
    NSYMAXD,IRID,NFUND, &
    NCELLD,LMXSPD,NGSHD,NUMN0,INDN0, &
    IGUESS,BCP,QMRBOUND, &
    NR,RCUTJIJ,JIJ,LDAU,ISYMINDEX)
! ======================================================================

! deallocations

    ! Energy mesh
    deallocate(EZ, stat=ierror)
    deallocate(WEZ, stat=ierror)
    deallocate(KMESH, stat=ierror)

    ! Lattice
    deallocate(RBASIS, stat=ierror)
    deallocate(RMT, stat=ierror)
    deallocate(RMTREF, stat=ierror)
    deallocate(RR, stat=ierror)
    deallocate(RWS, stat=ierror)
    deallocate(ZAT, stat=ierror)

    !Shape functions
    deallocate(THETAS, stat=ierror)
    deallocate(GSH, stat=ierror)
    deallocate(IFUNM, stat=ierror)
    deallocate(IPAN, stat=ierror)
    deallocate(LLMSP, stat=ierror)
    deallocate(LMSP, stat=ierror)
    deallocate(ILM, stat=ierror)
    deallocate(NFU, stat=ierror)
    deallocate(NTCELL, stat=ierror)
    deallocate(IMAXSH, stat=ierror)

    ! Reference Cluster
    deallocate(RCLS, stat=ierror)
    deallocate(ATOM, stat=ierror)
    deallocate(CLS, stat=ierror)
    deallocate(EZOA, stat=ierror)
    deallocate(NACLS, stat=ierror)
    deallocate(NUMN0, stat=ierror)
    deallocate(INDN0, stat=ierror)
    deallocate(REFPOT, stat=ierror)

    ! Reference system
    deallocate(VREF, stat=ierror)

    ! Radial mesh(es)
    deallocate(A, stat=ierror)
    deallocate(B, stat=ierror)
    deallocate(DRDI, stat=ierror)
    deallocate(R, stat=ierror)
    deallocate(IMT, stat=ierror)
    deallocate(IRC, stat=ierror)
    deallocate(IRCUT, stat=ierror)
    deallocate(IRMIN, stat=ierror)
    deallocate(IRNS, stat=ierror)
    deallocate(IRWS, stat=ierror)

    ! Symmetry matrices
    deallocate(DSYMLL, stat=ierror)
    deallocate(ISYMINDEX, stat=ierror)

    ! spherical harmonics
    deallocate(WG, stat=ierror)
    deallocate(YRG, stat=ierror)

    ! Clebsch-Gordon coefficients
    deallocate(CLEB, stat=ierror)
    deallocate(ICLEB, stat=ierror)
    deallocate(JEND, stat=ierror)
    deallocate(LOFLM, stat=ierror)

    !     .. core states ..
    deallocate(ITITLE, stat=ierror)
    deallocate(LCORE, stat=ierror)
    deallocate(NCORE, stat=ierror)

    !   auxillary
    deallocate(DEZ, stat=ierror)
    deallocate(RMTNEW, stat=ierror)
    deallocate(INIPOL, stat=ierror)

  end program

