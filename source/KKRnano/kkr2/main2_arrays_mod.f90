module main2_arrays_mod
  implicit none

  SAVE

  ! KKRnano-Arrays
  double complex, dimension(:), allocatable ::  EZ
  double complex, dimension(:), allocatable ::  WEZ
  double complex, dimension(:,:), allocatable ::  WEZRN
  double complex, dimension(:,:,:), allocatable ::  DEN
  double complex, dimension(:,:,:), allocatable ::  DSYMLL
  double complex, dimension(:,:), allocatable ::  PHILDAU
  double complex, dimension(:,:,:,:), allocatable ::  DMATLDAU ! LDA+U
  double precision, dimension(:), allocatable :: WG
  double precision, dimension(:,:,:), allocatable :: YRG
  double precision, dimension(:,:), allocatable :: RBASIS

  double precision, dimension(:,:), allocatable :: SMAT
  double precision, dimension(:), allocatable :: CLEB
  double precision, dimension(:,:), allocatable :: CLEB1C
  double precision, dimension(:,:), allocatable :: RNORM
  double precision, dimension(:,:), allocatable :: DFAC
  double precision, dimension(:), allocatable :: RWS
  double precision, dimension(:), allocatable :: RMT
  double precision, dimension(:,:), allocatable :: GN
  double precision, dimension(:,:), allocatable :: RM
  double precision, dimension(:,:,:), allocatable :: BZKP
  double precision, dimension(:,:), allocatable :: VOLCUB
  double precision, dimension(:), allocatable :: VOLBZ
  double precision, dimension(:,:), allocatable :: RR
  double precision, dimension(:,:,:), allocatable :: VINS        ! .. input potential
  double precision, dimension(:,:), allocatable :: VISP
  double precision, dimension(:,:,:), allocatable :: VONS        !     .. output potential
  double precision, dimension(:), allocatable :: ULDAU              ! LDA+U
  double precision, dimension(:), allocatable :: JLDAU              ! LDA+U
  double precision, dimension(:,:,:,:,:), allocatable :: UMLDAU     ! LDA+U
  double precision, dimension(:,:,:,:), allocatable :: WMLDAU

  integer, dimension(:), allocatable :: KMESH
  integer, dimension(:), allocatable :: NOFKS
  integer, dimension(:,:), allocatable :: EZOA

  integer, dimension(:), allocatable :: NUMN0
  integer, dimension(:,:), allocatable :: INDN0

  integer, dimension(:), allocatable :: NSG
  integer, dimension(:), allocatable :: NSR
  integer, dimension(:), allocatable :: LLDAU    ! LDA+U

  double complex, dimension(:,:,:), allocatable ::  TMATN
  double complex, dimension(:,:,:), allocatable ::  DTDE
  double complex, dimension(:,:,:), allocatable ::  TREFLL
  double complex, dimension(:,:,:), allocatable ::  DTREFLL
  double complex, dimension(:,:,:,:), allocatable ::  DGREFN
  double complex, dimension(:,:,:,:), allocatable ::  GREFN
  double complex, dimension(:,:,:,:), allocatable ::  GMATN
  double complex, dimension(:,:,:,:), allocatable ::  GMATN_ALL
  double complex, dimension(:,:), allocatable ::  DTIXIJ
  double complex, dimension(:,:,:,:), allocatable ::  GMATXIJ
  double complex, dimension(:,:,:,:), allocatable ::  GXIJ_ALL

  !----- Initial guess ---------------------------------------------------
  complex, dimension(:,:,:), allocatable :: PRSC
  double precision, dimension(:,:), allocatable ::  CNVFAC
  integer, dimension(:,:,:), allocatable :: SPRS

  !----- Lloyd -----------------------------------------------------------
  double complex, dimension(:,:), allocatable ::  LLY_G0TR
  double complex, dimension(:,:), allocatable ::  LLY_GRDT
  double complex, dimension(:,:), allocatable ::  LLY_GRDT_ALL
  double complex, dimension(:), allocatable ::  TR_ALPH
  double complex, dimension(:), allocatable ::  JXCIJINT
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  double precision, dimension(:), allocatable :: ECOU
  double precision, dimension(:,:), allocatable :: ESPC
  double precision, dimension(:,:), allocatable :: ESPV
  double precision, dimension(:), allocatable :: EXC

  !----------- Radial Mesh, shape functions ------------------------------
  double precision, dimension(:), allocatable :: A
  double precision, dimension(:), allocatable :: B
  double precision, dimension(:,:), allocatable :: DRDI
  double precision, dimension(:,:), allocatable :: R
  double precision, dimension(:,:,:), allocatable :: THETAS
  double precision, dimension(:), allocatable :: ZAT

  ! ----------------------------------------------------------------------
  double precision, dimension(:,:), allocatable ::  RHOCAT
  double precision, dimension(:,:,:), allocatable ::  R2NEF
  double precision, dimension(:,:,:), allocatable ::  RHO2NS
  double precision, dimension(:,:), allocatable ::  CHARGE
  !     ..
  double precision, dimension(:), allocatable :: GSH
  double precision, dimension(:), allocatable :: CMINST    ! charge moment of interstitial
  double precision, dimension(:), allocatable :: CMOM      ! LM moment of total charge
  double precision, dimension(:), allocatable :: CATOM     ! total charge per atom

  !     .. FORCES
  double precision, dimension(:,:), allocatable :: FLM
  double precision, dimension(:,:), allocatable :: FLMC
  !     .. MIXING
  double precision, dimension(:), allocatable :: SM1S
  double precision, dimension(:), allocatable :: FM1S
  double precision, dimension(:,:), allocatable :: UI2
  double precision, dimension(:,:), allocatable :: VI2
  double precision, dimension(:), allocatable :: WIT

  ! ----------------------------------------------------------------------
  integer, dimension(:), allocatable :: IMT
  integer, dimension(:), allocatable :: IPAN
  integer, dimension(:), allocatable :: IRC
  integer, dimension(:), allocatable :: IRNS
  integer, dimension(:,:), allocatable :: IRCUT
  integer, dimension(:), allocatable :: IRMIN
  integer, dimension(:), allocatable :: IRWS
  integer, dimension(:,:), allocatable :: ITITLE
  integer, dimension(:,:), allocatable :: LCORE
  integer, dimension(:,:), allocatable :: LLMSP
  integer, dimension(:), allocatable :: NCORE
  integer, dimension(:), allocatable :: NFU
  integer, dimension(:), allocatable :: NTCELL

  integer, dimension(:,:), allocatable :: ILM
  integer, dimension(:), allocatable :: IMAXSH
  integer, dimension(:,:), allocatable :: ICLEB
  integer, dimension(:), allocatable :: LOFLM
  integer, dimension(:,:), allocatable :: ICLEB1C
  integer, dimension(:), allocatable :: LOFLM1C
  integer, dimension(:,:), allocatable :: IFUNM

  integer, dimension(:,:), allocatable :: LMSP
  integer, dimension(:,:,:), allocatable :: JEND

  double precision, dimension(:,:,:), allocatable :: RCLS
  double precision, dimension(:), allocatable :: RMTREF
  double precision, dimension(:), allocatable :: VREF

  ! ------------- Jij calculation ---------------------------------------------
  double precision, dimension(:), allocatable :: RXIJ          ! interatomic distance Ri-Rj
  double precision, dimension(:,:), allocatable :: RXCCLS      ! position relative of j rel. to i (sorted)
  double precision, dimension(:,:,:), allocatable :: ZKRXIJ    ! set up in clsjij, used in kkrmat01
  integer, dimension(:), allocatable :: IXCP                   ! index to atom in elem/cell at site in cluster
  integer, dimension(:), allocatable :: NXCP                   ! index to bravais lattice at site in cluster

  integer, dimension(:,:), allocatable :: ATOM
  integer, dimension(:), allocatable :: CLS
  integer, dimension(:), allocatable :: NACLS
  integer, dimension(:), allocatable :: REFPOT

  ! -------------- MPI --------------------------------------------------------
  integer, dimension(:), allocatable :: MYLRANK
  integer, dimension(:), allocatable :: LCOMM
  integer, dimension(:), allocatable :: LGROUP
  integer, dimension(:), allocatable :: LSIZE

  !     .. LS-MPI
  integer, dimension(:,:), allocatable :: LSRANK
  integer, dimension(:,:), allocatable :: LSMYRANK

  !     .. S-MPI
  integer, dimension(:,:), allocatable :: SRANK
  integer, dimension(:,:), allocatable :: SMYRANK

  !     .. E-MPI
  real, dimension(:), allocatable :: ETIME
  integer, dimension(:,:), allocatable :: EMYRANK
  integer, dimension(:,:), allocatable :: ERANK
  integer, dimension(:), allocatable :: EPROC
  integer, dimension(:), allocatable :: EPROCO

  !==============================================================================
  ! Dimension Parameters.
  !==============================================================================

  integer :: NAEZD
  integer :: LMAXD
  integer :: NREFD
  integer :: IRID
  integer :: BCPD
  integer :: NACLSD
  integer :: NCLEB
  integer :: IRMD
  integer :: IEMXD
  integer :: NGSHD
  integer :: IGUESSD
  integer :: IPAND
  integer :: ISHLD
  integer :: IRNSD
  integer :: KPOIBZ
  integer :: NFUND
  integer :: NCLSD
  integer :: NMAXD
  integer :: NRD
  integer :: NSPIND
  integer :: NXIJD
  integer :: LLY
  integer :: EKMD
  integer :: NCELLD

  integer :: XDIM
  integer :: YDIM
  integer :: ZDIM
  integer :: NATBLD
  integer :: ITDBRYD

  ! ------------ Derived dimension variables ----------------------------------
  integer, parameter :: MAXMSHD = 8
  ! Parameters that changed to normal variables
  integer::   LMMAXD
  integer::   NPOTD
  integer::   LMAXD1
  integer::   MMAXD
  integer::   LM2D
  integer::   LMXSPD
  integer::   LASSLD
  integer::   LMPOTD
  integer::   IRMIND
  integer::   LRECPOT
  integer::   LRECRES1
  integer::   LRECRES2
  integer::   NTIRD    ! for Broyden mixing

  !derived parameters
  integer :: LPOTD
  integer :: NGUESSD

  !Parallelisation
  integer, parameter :: LMPID = 1  ! L-parallelisation not supported anymore
  integer :: SMPID
  integer :: EMPID
  integer :: NTHRDS


CONTAINS

  !---------------------------------------------------------------------
  ! Reads array dimension parameters (previously in inc.p and inc.cls)
  ! from file inp0.unf
  !---------------------------------------------------------------------
  ! missing: TRC (not used in kkr2)
  ! could be removed: KREL ?
  subroutine read_dimension_parameters()

    implicit none

    integer :: TRC ! not used anymore
    integer :: NATRCD ! not used
    integer :: NUTRCD ! not used
    integer :: FILEHANDLE = 67

    open (FILEHANDLE, FILE='inp0.unf', FORM='unformatted')

    read(FILEHANDLE) LMAXD
    read(FILEHANDLE) NSPIND
    read(FILEHANDLE) NAEZD
    read(FILEHANDLE) IRNSD
    read(FILEHANDLE) TRC
    read(FILEHANDLE) IRMD
    read(FILEHANDLE) NREFD
    read(FILEHANDLE) NRD
    read(FILEHANDLE) IRID
    read(FILEHANDLE) NFUND
    read(FILEHANDLE) NCELLD
    read(FILEHANDLE) NGSHD
    read(FILEHANDLE) NACLSD
    read(FILEHANDLE) NCLSD
    read(FILEHANDLE) IPAND
    read(FILEHANDLE) NXIJD
    read(FILEHANDLE) NATRCD
    read(FILEHANDLE) NUTRCD
    read(FILEHANDLE) KPOIBZ
    read(FILEHANDLE) IGUESSD
    read(FILEHANDLE) BCPD
    read(FILEHANDLE) NMAXD
    read(FILEHANDLE) ISHLD
    read(FILEHANDLE) LLY
    read(FILEHANDLE) SMPID
    read(FILEHANDLE) EMPID
    read(FILEHANDLE) NTHRDS
    read(FILEHANDLE) XDIM
    read(FILEHANDLE) YDIM
    read(FILEHANDLE) ZDIM
    read(FILEHANDLE) NATBLD
    read(FILEHANDLE) ITDBRYD
    read(FILEHANDLE) IEMXD
    read(FILEHANDLE) EKMD

    close(FILEHANDLE)

  end subroutine read_dimension_parameters

  !----------------------------------------------------------------------------
  subroutine allocate_main2_arrays()
    use ErrorMessages_mod
    implicit none

    integer :: memory_stat

    allocate(EZ(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(WEZ(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(WEZRN(IEMXD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DEN(0:LMAXD1,IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DSYMLL(LMMAXD,LMMAXD,48), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(PHILDAU(IRMD,LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(WG(LASSLD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(YRG(LASSLD,0:LASSLD,0:LASSLD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RBASIS(3,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(SMAT(LMXSPD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CLEB(LMXSPD*LMPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CLEB1C(NCLEB,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RNORM(IEMXD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DFAC(0:LPOTD,0:LPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RWS(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RMT(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GN(3,NMAXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RM(3,NMAXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(BZKP(3,KPOIBZ,MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VOLCUB(KPOIBZ,MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VOLBZ(MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RR(3,0:NRD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VINS(IRMIND:IRMD,LMPOTD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VISP(IRMD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VONS(IRMD,LMPOTD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ULDAU(LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(JLDAU(LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(UMLDAU(MMAXD,MMAXD,MMAXD,MMAXD,LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(KMESH(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NOFKS(MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EZOA(NACLSD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NUMN0(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(INDN0(NAEZD,NACLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NSG(ISHLD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NSR(ISHLD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLDAU(LMAXD1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(TMATN(LMMAXD,LMMAXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DTDE(LMMAXD,LMMAXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(TREFLL(LMMAXD,LMMAXD,NREFD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DTREFLL(LMMAXD,LMMAXD,NREFD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DGREFN(LMMAXD,LMMAXD,NACLSD,NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GREFN(LMMAXD,LMMAXD,NACLSD,NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GMATN_ALL(LMMAXD,LMMAXD,IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DTIXIJ(LMMAXD,LMMAXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GMATXIJ(LMMAXD,LMMAXD,NXIJD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GXIJ_ALL(LMMAXD,LMMAXD,NXIJD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(PRSC(NGUESSD*LMMAXD,EKMD,NSPIND-SMPID+1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CNVFAC(EKMD,NSPIND-SMPID+1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(SPRS(NGUESSD*LMMAXD+1,EKMD+1,NSPIND-SMPID+1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLY_G0TR(IEMXD,NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLY_GRDT(IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLY_GRDT_ALL(IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(TR_ALPH(NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(JXCIJINT(NXIJD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ECOU(0:LPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ESPC(0:3,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ESPV(0:LMAXD1,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EXC(0:LPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(A(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(B(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DRDI(IRMD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(R(IRMD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(THETAS(IRID,NFUND,NCELLD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ZAT(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RHOCAT(IRMD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(R2NEF(IRMD,LMPOTD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RHO2NS(IRMD,LMPOTD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CHARGE(0:LMAXD1,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(GSH(NGSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CMINST(LMPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CMOM(LMPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CATOM(NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(FLM(-1:1,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(FLMC(-1:1,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(SM1S(NTIRD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(FM1S(NTIRD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(UI2(NTIRD,2:ITDBRYD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VI2(NTIRD,2:ITDBRYD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(WIT(2:ITDBRYD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IMT(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IPAN(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IRC(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IRNS(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IRCUT(0:IPAND,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IRMIN(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IRWS(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ITITLE(20,NPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LCORE(20,NPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLMSP(NFUND,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NCORE(NPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NFU(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NTCELL(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ILM(NGSHD,3), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IMAXSH(0:LMPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ICLEB(LMXSPD*LMPOTD,3), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LOFLM(LMXSPD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ICLEB1C(NCLEB,3), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LOFLM1C(LM2D), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IFUNM(LMXSPD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LMSP(LMXSPD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(JEND(LMPOTD,0:LMAXD,0:LMAXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RCLS(3,NACLSD,NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RMTREF(NREFD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VREF(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RXIJ(NXIJD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RXCCLS(3,NXIJD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ZKRXIJ(48,3,NXIJD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(IXCP(NXIJD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NXCP(NXIJD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ATOM(NACLSD,NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CLS(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NACLS(NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(REFPOT(NAEZD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(MYLRANK(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LCOMM(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LGROUP(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LSIZE(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LSRANK(LMPID,NAEZD*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LSMYRANK(LMPID,NAEZD*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(SRANK(SMPID,NAEZD*LMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(SMYRANK(SMPID,NAEZD*LMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ETIME(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EMYRANK(EMPID,NAEZD*LMPID*SMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ERANK(EMPID,NAEZD*LMPID*SMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EPROC(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EPROCO(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")

  end subroutine

  !------------------------------------------------------------------------------
  subroutine deallocate_main2_arrays()
    implicit none

    integer :: memory_stat

    deallocate(EZ, stat = memory_stat)
    deallocate(WEZ, stat = memory_stat)
    deallocate(WEZRN, stat = memory_stat)
    deallocate(DEN, stat = memory_stat)
    deallocate(DSYMLL, stat = memory_stat)
    deallocate(PHILDAU, stat = memory_stat)
    deallocate(DMATLDAU, stat = memory_stat)
    deallocate(WG, stat = memory_stat)
    deallocate(YRG, stat = memory_stat)
    deallocate(RBASIS, stat = memory_stat)
    deallocate(SMAT, stat = memory_stat)
    deallocate(CLEB, stat = memory_stat)
    deallocate(CLEB1C, stat = memory_stat)
    deallocate(RNORM, stat = memory_stat)
    deallocate(DFAC, stat = memory_stat)
    deallocate(RWS, stat = memory_stat)
    deallocate(RMT, stat = memory_stat)
    deallocate(GN, stat = memory_stat)
    deallocate(RM, stat = memory_stat)
    deallocate(BZKP, stat = memory_stat)
    deallocate(VOLCUB, stat = memory_stat)
    deallocate(VOLBZ, stat = memory_stat)
    deallocate(RR, stat = memory_stat)
    deallocate(VINS, stat = memory_stat)
    deallocate(VISP, stat = memory_stat)
    deallocate(VONS, stat = memory_stat)
    deallocate(ULDAU, stat = memory_stat)
    deallocate(JLDAU, stat = memory_stat)
    deallocate(UMLDAU, stat = memory_stat)
    deallocate(WMLDAU, stat = memory_stat)
    deallocate(KMESH, stat = memory_stat)
    deallocate(NOFKS, stat = memory_stat)
    deallocate(EZOA, stat = memory_stat)
    deallocate(NUMN0, stat = memory_stat)
    deallocate(INDN0, stat = memory_stat)
    deallocate(NSG, stat = memory_stat)
    deallocate(NSR, stat = memory_stat)
    deallocate(LLDAU, stat = memory_stat)
    deallocate(TMATN, stat = memory_stat)
    deallocate(DTDE, stat = memory_stat)
    deallocate(TREFLL, stat = memory_stat)
    deallocate(DTREFLL, stat = memory_stat)
    deallocate(DGREFN, stat = memory_stat)
    deallocate(GREFN, stat = memory_stat)
    deallocate(GMATN, stat = memory_stat)
    deallocate(GMATN_ALL, stat = memory_stat)
    deallocate(DTIXIJ, stat = memory_stat)
    deallocate(GMATXIJ, stat = memory_stat)
    deallocate(GXIJ_ALL, stat = memory_stat)
    deallocate(PRSC, stat = memory_stat)
    deallocate(CNVFAC, stat = memory_stat)
    deallocate(SPRS, stat = memory_stat)
    deallocate(LLY_G0TR, stat = memory_stat)
    deallocate(LLY_GRDT, stat = memory_stat)
    deallocate(LLY_GRDT_ALL, stat = memory_stat)
    deallocate(TR_ALPH, stat = memory_stat)
    deallocate(JXCIJINT, stat = memory_stat)
    deallocate(ECOU, stat = memory_stat)
    deallocate(ESPC, stat = memory_stat)
    deallocate(ESPV, stat = memory_stat)
    deallocate(EXC, stat = memory_stat)
    deallocate(A, stat = memory_stat)
    deallocate(B, stat = memory_stat)
    deallocate(DRDI, stat = memory_stat)
    deallocate(R, stat = memory_stat)
    deallocate(THETAS, stat = memory_stat)
    deallocate(ZAT, stat = memory_stat)
    deallocate(RHOCAT, stat = memory_stat)
    deallocate(R2NEF, stat = memory_stat)
    deallocate(RHO2NS, stat = memory_stat)
    deallocate(CHARGE, stat = memory_stat)
    deallocate(GSH, stat = memory_stat)
    deallocate(CMINST, stat = memory_stat)
    deallocate(CMOM, stat = memory_stat)
    deallocate(CATOM, stat = memory_stat)
    deallocate(FLM, stat = memory_stat)
    deallocate(FLMC, stat = memory_stat)
    deallocate(SM1S, stat = memory_stat)
    deallocate(FM1S, stat = memory_stat)
    deallocate(UI2, stat = memory_stat)
    deallocate(VI2, stat = memory_stat)
    deallocate(WIT, stat = memory_stat)
    deallocate(IMT, stat = memory_stat)
    deallocate(IPAN, stat = memory_stat)
    deallocate(IRC, stat = memory_stat)
    deallocate(IRNS, stat = memory_stat)
    deallocate(IRCUT, stat = memory_stat)
    deallocate(IRMIN, stat = memory_stat)
    deallocate(IRWS, stat = memory_stat)
    deallocate(ITITLE, stat = memory_stat)
    deallocate(LCORE, stat = memory_stat)
    deallocate(LLMSP, stat = memory_stat)
    deallocate(NCORE, stat = memory_stat)
    deallocate(NFU, stat = memory_stat)
    deallocate(NTCELL, stat = memory_stat)
    deallocate(ILM, stat = memory_stat)
    deallocate(IMAXSH, stat = memory_stat)
    deallocate(ICLEB, stat = memory_stat)
    deallocate(LOFLM, stat = memory_stat)
    deallocate(ICLEB1C, stat = memory_stat)
    deallocate(LOFLM1C, stat = memory_stat)
    deallocate(IFUNM, stat = memory_stat)
    deallocate(LMSP, stat = memory_stat)
    deallocate(JEND, stat = memory_stat)
    deallocate(RCLS, stat = memory_stat)
    deallocate(RMTREF, stat = memory_stat)
    deallocate(VREF, stat = memory_stat)
    deallocate(RXIJ, stat = memory_stat)
    deallocate(RXCCLS, stat = memory_stat)
    deallocate(ZKRXIJ, stat = memory_stat)
    deallocate(IXCP, stat = memory_stat)
    deallocate(NXCP, stat = memory_stat)
    deallocate(ATOM, stat = memory_stat)
    deallocate(CLS, stat = memory_stat)
    deallocate(NACLS, stat = memory_stat)
    deallocate(REFPOT, stat = memory_stat)
    deallocate(MYLRANK, stat = memory_stat)
    deallocate(LCOMM, stat = memory_stat)
    deallocate(LGROUP, stat = memory_stat)
    deallocate(LSIZE, stat = memory_stat)
    deallocate(LSRANK, stat = memory_stat)
    deallocate(LSMYRANK, stat = memory_stat)
    deallocate(SRANK, stat = memory_stat)
    deallocate(SMYRANK, stat = memory_stat)
    deallocate(ETIME, stat = memory_stat)
    deallocate(EMYRANK, stat = memory_stat)
    deallocate(ERANK, stat = memory_stat)
    deallocate(EPROC, stat = memory_stat)
    deallocate(EPROCO, stat = memory_stat)

  end subroutine

end module main2_arrays_mod
