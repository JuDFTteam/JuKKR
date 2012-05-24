! KKRnano
! massive parallel KKR for nanoscaled systems

program MAIN2

  !use mpi
  use common_testc
  use common_optc
  use common_mpi

  use ErrorMessages_mod
  use lloyds_formula_mod
  use KKRSelfConsistency_mod

  use main2_aux_mod

  implicit none
  include 'mpif.h'

  !     .. Parameters ..

  integer::   MAXMSHD
  parameter (MAXMSHD=8)
  integer::   NSYMAXD
  parameter (NSYMAXD=48)
  double complex:: CZERO
  parameter      (CZERO=(0.0D0,0.0D0))

  !     ..
  !     .. Local Scalars ..

  double complex:: JSCAL        ! scaling factor for Jij calculation
  double precision::DENEF
  double precision::E1
  double precision::E2
  double precision::TK
  double precision::EFERMI
  double precision::EBOT
  double precision::EFOLD
  double precision::CHRGNT
  double precision::E2SHIFT
  double precision::RCUTJIJ
  double precision::DF
  double precision::ALAT
  double precision::FCM
  double precision::MIX
  double precision::MIXING
  double precision::QBOUND
  double precision::VOLUME0
  double precision::RMAX
  double precision::GMAX
  double precision::FPI
  double precision::PI
  double precision::RFPI
  double precision::RMSAVM      ! rms error magnetisation dens. (contribution of single site)
  double precision::RMSAVQ      ! rms error charge density (contribution of single site)
  double precision::EREFLDAU    ! LDA+U

  double precision::WALLCLOCK_I ! time management
  double precision::WALLCLOCK_F
  real::TIME_I
  real::TIME_S
  real::TIME_E
  real::TIME_EX

  integer::ICELL
  integer::IR
  integer::LMAX
  integer::NPNT1
  integer::NPNT2
  integer::NPNT3
  integer::NPOL
  integer::NSPIN
  integer::ITER
  integer::SCFSTEPS
  integer::IMIX
  integer::NOITER
  integer::NOITER_ALL
  integer::ISHIFT
  integer::KPRE
  integer::KTE
  integer::KVMAD
  integer::KXC
  integer::KFORCE
  integer::LPOT
  integer::LMPOT
  integer::NAEZ
  integer::IEND
  integer::NCLEBD
  integer::LM1
  integer::LM2
  integer::IEND1
  integer::IPOT
  integer::ISPIN
  integer::I1
  integer::I1BRYD
  integer::LM
  integer::NR
  integer::EKM
  integer::SYSTEM_I
  integer::SYSTEM_F
  integer::RATETIME
  integer::MAXTIME
  logical::TEST
  logical::LCORDENS
  logical::XCCPL
  logical::JIJ
  logical::STOPIT
  logical::LDORHOEF
  logical::LDAU
  logical::ERESJIJ


  ! static arrays
  double precision::BRAVAIS(3,3)
  double precision::RECBV(3,3)
  double precision::VBC(2)
  integer::ISYMINDEX(NSYMAXD)
  double precision::ECORE(20,2)
  double precision:: WORK1(2)
  double precision:: WORK2(2)

  !     .. Local Arrays ..
  double complex, dimension(:), allocatable ::  EZ
  double complex, dimension(:), allocatable ::  WEZ
  double complex, dimension(:), allocatable ::  DEZ
  double complex, dimension(:,:), allocatable ::  WEZRN
  double complex, dimension(:,:,:), allocatable ::  DEN
  double complex, dimension(:,:,:), allocatable ::  DSYMLL
  double complex, dimension(:,:), allocatable ::  PHILDAU
  double complex, dimension(:,:,:,:), allocatable ::  DMATLDAU ! LDA+U
  double precision, dimension(:), allocatable :: WG
  double precision, dimension(:,:,:), allocatable :: YRG
  double precision, dimension(:,:), allocatable :: RBASIS
  double precision::VAV0
  double precision::VOL0
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

  integer::NMESH
  integer, dimension(:), allocatable :: KMESH
  integer, dimension(:), allocatable :: NOFKS
  integer, dimension(:,:), allocatable :: EZOA
  integer::NSYMAT
  integer, dimension(:), allocatable :: NUMN0
  integer, dimension(:,:), allocatable :: INDN0
  integer::MAXMESH
  integer, dimension(:), allocatable :: NSG
  integer, dimension(:), allocatable :: NSR
  integer, dimension(:), allocatable :: LLDAU    ! LDA+U
  integer::NGMAX
  integer::NRMAX
  integer::NSHLG
  integer::NSHLR

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
  double precision:: QMRBOUND
  double precision, dimension(:,:), allocatable ::  CNVFAC
  integer::IGUESS
  integer::BCP
  integer::PRSPIN
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
  double precision::EPOTIN
  double precision, dimension(:,:), allocatable :: ESPC
  double precision, dimension(:,:), allocatable :: ESPV
  double precision, dimension(:), allocatable :: EXC
  double precision::EULDAU
  double precision::EDCLDAU
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
  double precision::QC                ! core charge
  double precision::VMAD
  !     ,,
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
  integer::LCOREMAX
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
  integer::ICST
  integer::NSRA
  integer, dimension(:,:), allocatable :: LMSP
  integer, dimension(:,:,:), allocatable :: JEND
  integer::NCLS
  integer::NREF
  integer::RF
  integer::NLDAU
  double precision, dimension(:,:,:), allocatable :: RCLS
  double precision, dimension(:), allocatable :: RMTREF
  double precision, dimension(:), allocatable :: VREF

  double precision, dimension(:), allocatable :: RXIJ          ! interatomic distance Ri-Rj
  double precision, dimension(:,:), allocatable :: RXCCLS      ! position relative of j rel. to i (sorted)
  double precision, dimension(:,:,:), allocatable :: ZKRXIJ    ! set up in clsjij, used in kkrmat01
  integer, dimension(:), allocatable :: IXCP                   ! index to atom in elem/cell at site in cluster
  integer, dimension(:), allocatable :: NXCP                   ! index to bravais lattice at site in cluster
  integer::XIJ
  integer::NXIJ
  integer, dimension(:,:), allocatable :: ATOM
  integer, dimension(:), allocatable :: CLS
  integer, dimension(:), allocatable :: NACLS
  integer, dimension(:), allocatable :: REFPOT

  !     .. L-MPI
  integer, dimension(:), allocatable :: MYLRANK
  integer, dimension(:), allocatable :: LCOMM
  integer, dimension(:), allocatable :: LGROUP
  integer, dimension(:), allocatable :: LSIZE
  integer::LMPIC

  !     .. LS-MPI
  integer, dimension(:,:), allocatable :: LSRANK
  integer, dimension(:,:), allocatable :: LSMYRANK
  integer::LSMPIC
  integer::LSMPIB
  !     .. S-MPI
  integer, dimension(:,:), allocatable :: SRANK
  integer, dimension(:,:), allocatable :: SMYRANK
  integer::SMPIC
  integer::SMPIB
  integer::MAPSPIN

  !     .. E-MPI
  real, dimension(:), allocatable :: ETIME
  integer, dimension(:,:), allocatable :: EMYRANK
  integer, dimension(:,:), allocatable :: ERANK
  integer::EMPIC
  integer::EMPIB
  integer::IE
  integer::IELAST
  integer, dimension(:), allocatable :: EPROC
  integer, dimension(:), allocatable :: EPROCO

  !     .. ACTV-MPI
  integer::MYACTVRANK
  integer::ACTVCOMM
  integer::ACTVGROUP
  integer::ACTVSIZE
  integer::MYBCRANK
  integer::BCRANK

  integer::   IERR
  integer::   MAPBLOCK
  external     MAPBLOCK

  ! Array allocations
  integer:: memory_stat
    !logical :: memory_fail

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

  ! dimension parameters
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

  !Parallelisation
  integer, parameter :: LMPID = 1  ! L-parallelisation not supported anymore
  integer :: SMPID
  integer :: EMPID
  integer :: NTHRDS

  !derived parameters
  integer :: LPOTD
  integer :: NGUESSD

!============================================================= CONSTANTS
  LCORDENS=.true.
  PI = 4.0D0*ATAN(1.0D0)
  FPI = 4.0D0*PI
  RFPI = SQRT(FPI)
!=============================================================

  call read_dimension_parameters( LMAXD, NSPIND, NAEZD, IRNSD, &
                                  IRMD, NREFD, NRD, IRID, NFUND, NCELLD, &
                                  NGSHD, NACLSD, NCLSD, IPAND, NXIJD, KPOIBZ, &
                                  IGUESSD, BCPD, NMAXD, ISHLD, &
                                  LLY, SMPID, EMPID, NTHRDS, XDIM, YDIM, ZDIM, &
                                  NATBLD, ITDBRYD, IEMXD, EKMD)

  ! from dimension parameters - calculate some derived parameters
  call getDerivedParameters(IGUESSD, IRMD, IRMIND, IRNSD, LASSLD, LM2D, LMAXD, &
                            LMAXD1, LMMAXD, LMPOTD, LMXSPD, LPOTD, LRECPOT, &
                            LRECRES2, MMAXD, NAEZD, NCLEB, NGUESSD, NPOTD, NSPIND, NTIRD)

  ! consistency check of some dimension parameters
  call consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)

!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  allocate(EZ(IEMXD), stat = memory_stat)
  if(memory_stat /= 0) call fatalMemoryError("main2")
  allocate(WEZ(IEMXD), stat = memory_stat)
  if(memory_stat /= 0) call fatalMemoryError("main2")
  allocate(DEZ(IEMXD), stat = memory_stat)
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

  !-----------------------------------------------------------------------------
  ! Array allocations END
  !-----------------------------------------------------------------------------

  call readKKR0Input      (NSYMAXD, A, ALAT, ATOM, B, BCP, BRAVAIS, CLEB1C, &
                           CLS, DRDI, DSYMLL, EZOA, FCM, GMAX, GSH, ICLEB1C, ICST, &
                           IEND1, IFUNM, IGUESS, ILM, IMAXSH, IMIX, IMT, INDN0, IPAN, &
                           IRC, IRCUT, IRMIN, IRNS, IRWS, ISHIFT, ISYMINDEX, ITITLE, &
                           JEND, JIJ, KFORCE, KMESH, KPRE, KTE, KVMAD, KXC, LCORE, &
                           LDAU, LLMSP, LMAX, LMPOT, LMSP, LOFLM1C, LPOT, MAXMESH, &
                           MIXING, NACLS, NAEZ, NCLS, NCORE, NFU, NR, NREF, NSPIN, &
                           NSRA, NSYMAT, NTCELL, NUMN0, OPTC, QBOUND, QMRBOUND, R, &
                           RBASIS, RCLS, RCUTJIJ, RECBV, REFPOT, RMAX, RMT, RMTREF, &
                           RR, RWS, SCFSTEPS, TESTC, THETAS, VOLUME0, VREF, ZAT)

  ! ---------------------------------------------------------- k_mesh
  call readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)

  call readEnergyMesh(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)

  if (KFORCE==1) open (54,file='force',form='formatted')   ! every process opens file 'force' !!!

 ! ======================================================================
 ! =                     End read in variables                          =
 ! ======================================================================

  call consistencyCheck02(IELAST, IEMXD, IGUESS, IGUESSD, LMAX, LMAXD, NAEZ, NAEZD, &
                          NPNT1, NPNT2, NPNT3, NPOL, NR, NRD, NSPIN, NSPIND)

  call consistencyCheck03(ATOM, CLS, EZOA, INDN0, NACLS, NACLSD, NAEZ, NCLSD, NR, NUMN0)

! ------------------------------------------------------------------

  ! TODO: get rid of LSMYRANK, LSRANK, LSMPIB, LSMPIC
  call IMPI(NAEZ,MYRANK,NROFNODES, &
            LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE, &
            LSMPIB,LSMPIC,LSRANK,LSMYRANK, &
            SMPIB,SMPIC,SRANK,SMYRANK, &
            EMPIB,EMPIC,ERANK,EMYRANK, &
            MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE, &
            lmpid, smpid, empid, nthrds)

!====================================================================

!=====================================================================
!     processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ...
! ... and wait after SC-ITER loop
!=====================================================================


  if (LMPIC/=0.or.LSMPIC/=0) then   !     ACTVGROUP could also test EMPIC

    MYBCRANK = 0

! ========= TIMING ======================================================
    if (MYLRANK(1) == 0) then
      RATETIME = 100
      MAXTIME  = 100000

      call SYSTEM_CLOCK(SYSTEM_I,RATETIME,MAXTIME)

      WALLCLOCK_I = MPI_WTIME()

      call CPU_TIME(TIME_I)

      open (2,file='time-info',form='formatted')
    endif
!========= TIMING END ======================================================

    CNVFAC = 1000.0D0

! initialise the arrays for (gen. Anderson/Broyden) potential mixing
    UI2 = 0.00
    VI2 = 0.00
    SM1S = 0.00
    FM1S = 0.00

! ######################################################################
! ######################################################################
    do ITER = 1, SCFSTEPS
! ######################################################################
! ######################################################################

      call CPU_TIME(TIME_S)

      EKM    = 0
      NOITER = 0

      if (MYLRANK(1)==0) then
        write(2,'(79(1H=))')
        call OUTTIME(MYLRANK(1),'started at ..........',TIME_I,ITER)
        write(2,'(79(1H=))')
      endif

      call GAUNT2(WG,YRG,LMAX)

      call MADELUNG3D(LPOT,YRG,WG,ALAT, &
      RMAX,GMAX,BRAVAIS,RECBV, &
      LMXSPD,LASSLD,LPOTD,LMPOTD, &
      NMAXD,ISHLD, &
      LMPOT,CLEB,ICLEB,IEND, &
      NCLEBD,LOFLM,DFAC, &
      NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM) ! does it have to be in SCF-loop?

      do LM = 1,LMPOTD
        CMOM(LM) = 0.0D0
        CMINST(LM) = 0.0D0
      end do

      CHRGNT = 0.0D0

      LRECRES1 = 8*43 + 16*(LMAXD+2)
      if (NPOL==0 .or. TEST('DOS     ')) then
        LRECRES1 = LRECRES1 + 32*(LMAXD+2)*IEMXD
      end if

      ! needed for results.f - find better solution - unnecessary I/O
      open (71,access='direct',recl=LRECRES1,file='results1', &
      form='unformatted')
      open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
      form='unformatted')


!N ====================================================================
!     BEGIN do loop over atoms (NMPID-parallel)
!N ====================================================================

      do I1 = 1,NAEZ
        if(MYLRANK(LMPIC)==MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

!=======================================================================
! xccpl

          XCCPL = .false.
          ERESJIJ = .false.

          ! calculate exchange couplings only at last self-consistency step and when Jij=true
          if ((ITER==SCFSTEPS).and.JIJ) XCCPL = .true.

          if (XCCPL) then

            inquire(file='ERESJIJ',exist=ERESJIJ)

            call CLSJIJ(I1,NAEZ,RR,NR,RBASIS,RCUTJIJ,NSYMAT,ISYMINDEX, &
                        IXCP,NXCP,NXIJ,RXIJ,RXCCLS,ZKRXIJ, &
                        nrd, nxijd)

            do ISPIN = 1, NSPIN
              do XIJ = 1, NXIJ
                do LM1 = 1,LMMAXD
                  do LM2 = 1,LMMAXD
                    GMATXIJ(LM1,LM2,XIJ,ISPIN) = CZERO
                  enddo
                enddo
              enddo
            enddo

          endif

! xccpl
!=======================================================================

          ! This read is probably NOT NECESSARY (except for 1st iteration) !!!
          ! TURNS out it is necessary - otherwise wrong results - why?
          read(66,rec=I1) VINS,VISP,ECORE  ! Read potential from file!!!

! LDA+U
          if (LDAU) then

            EREFLDAU = EFERMI
            EREFLDAU = 0.48       ! FIXME: hardcoded

            call LDAUINIT(I1,ITER,NSRA,NLDAU,LLDAU,ULDAU,JLDAU,EREFLDAU, &
                          VISP,NSPIN,R(1,I1),DRDI(1,I1), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          PHILDAU,UMLDAU,WMLDAU, &
                          lmax, irmd, ipand)

          endif
! LDA+U

! TIME
          call OUTTIME(MYLRANK(1),'initialized .........',TIME_I,ITER)
! TIME

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================

          do IE = 1,IELAST

            call CPU_TIME(TIME_E)

            ETIME(IE) = 0.0D0

            if (ITER==1.and.IE==1) then

              call EBALANCE('I',ITER,SCFSTEPS, &
                            IELAST,NPNT1, &
                            MYACTVRANK,ACTVCOMM, &
                            ETIME,EPROC,EPROCO, &
                            empid, iemxd)

            endif

! IE ====================================================================
            if (EMPIB==EPROC(IE)) then
! IE ====================================================================


              do RF = 1,NREF

                call TREF(EZ(IE),VREF(RF),LMAX,RMTREF(RF), &
                          TREFLL(1,1,RF),DTREFLL(1,1,RF), LLY)

              end do

              call GREF(EZ(IE),ALAT,IEND1,NCLS,NAEZ, &
                        CLEB1C,RCLS,ATOM,CLS,ICLEB1C,LOFLM1C,NACLS, &
                        REFPOT, &
                        TREFLL(1,1,1),DTREFLL(1,1,1),GREFN,DGREFN, &
                        IE, &
                        LLY_G0TR, &
                        LMPIC,MYLRANK,LCOMM,LSIZE, &
                        naez, lmax, naclsd, ncleb, nrefd, iemxd, nclsd, &
                        LLY, LMPID*SMPID*EMPID)

spinloop:     do ISPIN = 1,NSPIN

                call CALCTMAT(LDAU,NLDAU,ICST, &
                              NSRA,EZ(IE), &
                              DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                              VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                              IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                              TMATN(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN, &
                              LLDAU,WMLDAU, &
                              nspind, ncleb, ipand, irmd, irnsd)

                if(LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
                  call CALCDTMAT(LDAU,NLDAU,ICST, &
                                NSRA,EZ(IE),IE,NPNT1,NPNT2,NPNT3,PI,TK, &
                                DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                                DTDE(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN, &
                                LLDAU,WMLDAU, &
                                nspind, ncleb, ipand, irmd, irnsd)
                end if

                ! calculate DTIXIJ = T_down - T_up
                if (XCCPL) then
                  if (ISPIN==1) then
                    do LM1 = 1,LMMAXD
                      do LM2 = 1,LMMAXD
                        DTIXIJ(LM1,LM2) = TMATN(LM1,LM2,ISPIN)
                      enddo
                    enddo
                  else
                    do LM1 = 1,LMMAXD
                      do LM2 = 1,LMMAXD
                        DTIXIJ(LM1,LM2) = DTIXIJ(LM1,LM2)-TMATN(LM1,LM2,ISPIN)
                      enddo
                    enddo
                  endif
                endif


                RF = REFPOT(I1)
                do LM1 = 1,LMMAXD
                  TMATN(LM1,LM1,ISPIN) =  TMATN(LM1,LM1,ISPIN) &
                  - TREFLL(LM1,LM1,RF)
                  DTDE(LM1,LM1,ISPIN) =  DTDE(LM1,LM1,ISPIN) &
                  - DTREFLL(LM1,LM1,RF)
                end do

                ! TMATN now contains Delta t = t - t_ref !!!
                ! DTDE now contains Delta dt !!!


! SPIN ==================================================================
!     BEGIN do loop over spins (SMPID-parallel)
! SPIN===================================================================

                ! renormalize TR_ALPH
                TR_ALPH(ISPIN) = TR_ALPH(ISPIN) - LLY_G0TR(IE,CLS(I1))
                LLY_GRDT(IE,ISPIN) = CZERO ! initialize LLY_GRDT, shouldn't be necessary

                if (SMPID==1) then
                  MAPSPIN = 0
                  PRSPIN   = ISPIN
                else
                  MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
                  PRSPIN   = 1
                endif

!         true beginning of SMPID-parallel section

                if(SRANK(SMPIB,SMPIC)==MAPSPIN) then

                  NMESH = KMESH(IE)

                  if( MYLRANK(LMPIC)==0 ) then
                    write (6,'(A,I3,A,2(1X,F10.6),A,I3,A,I3)')  &
                    ' ** IE = ',IE,' ENERGY =',EZ(IE), &
                    ' KMESH = ', NMESH,' ISPIN = ',ISPIN
                  end if

! <<>>

                  call KLOOPZ1( &
                  GMATN(1,1,1,ISPIN), &
                  ALAT,IE,ITER,NAEZ, &
                  NOFKS(NMESH),VOLBZ(NMESH), &
                  BZKP(1,1,NMESH),VOLCUB(1,NMESH), &
                  CLS,NACLS,RR, &
                  EZOA,ATOM,GREFN,DGREFN, &
                  NSYMAT,DSYMLL, &
                  TMATN(:,:,ISPIN),DTDE(:,:,ISPIN), &
                  NUMN0,INDN0,I1, &
                  SPRS(1,1,PRSPIN),PRSC(1,1,PRSPIN), &
                  EKM,NOITER, &
                  QMRBOUND,IGUESS,BCP,CNVFAC(1,PRSPIN), &
                  NXIJ,XCCPL,IXCP,ZKRXIJ, &
                  LLY_GRDT(IE,ISPIN),TR_ALPH(ISPIN), &
                  GMATXIJ(1,1,1,ISPIN), &
                  LCOMM(LMPIC),LSIZE(LMPIC), &
                  iemxd, nthrds, &
                  lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
                  nxijd, nguessd, kpoibz, nrd, ekmd)

                endif

              end do spinloop                          ! ISPIN = 1,NSPIN

! SPIN ==================================================================
!     END do loop over spins (SMPID-parallel)
! SPIN===================================================================

! =====================================================================
! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
! xccpl

              if (XCCPL) then

                call SREDGX( ISPIN,NSPIN, &
                             MYRANK, &
                             SMPIC,SMYRANK, &
                             GMATXIJ, &
                             GXIJ_ALL, &
                             naez, lmax, lmpid, empid, smpid, nxijd)

                JSCAL = WEZ(IE)/DBLE(NSPIN)

                call XCCPLJIJ_START(I1,IE,JSCAL, &
                               RXIJ,NXIJ,IXCP,RXCCLS, &
                               GXIJ_ALL,DTIXIJ, &
                               LMPIC,LCOMM, &
                               JXCIJINT,ERESJIJ, &
                               naez, lmmaxd, nxijd, nspind, &
                               lmpid*smpid*empid)

              end if

! xccpl
! End of Jij calculation
! =====================================================================

              call CPU_TIME(TIME_EX)
              ETIME(IE) = TIME_EX-TIME_E

! IE ====================================================================
            endif
! IE ====================================================================

! for initial guess calculate sparse indices combining IE.KPT
            EKM = EKM + NOFKS(KMESH(IE))

          end do                   ! IE = 1,IELAST

! IE ====================================================================
!     END do loop over energies (EMPID-parallel) to be implemented
! IE ====================================================================


!=======================================================================
!     "allreduce" information of 1 .. EMPID and 1 .. SMPID processors
!=======================================================================
          call SREDGM(NSPIN,IELAST, &
                      MYRANK, &
                      SMPIC,SMYRANK, &
                      EMPIC,EMYRANK,EPROC, &
                      GMATN,LLY_GRDT, &
                      GMATN_ALL,LLY_GRDT_ALL, &
                      naez, lmax, lmpid, smpid, empid, iemxd)
!=======================================================================
!=======================================================================

! TIME
          call OUTTIME(MYLRANK(1),'G obtained ..........',TIME_I,ITER)
! TIME

!=======================================================================
!     output of Jij's .. calling xccpljij with flag 'F'
!=======================================================================
          if (XCCPL) then

            call XCCPLJIJ_OUT(I1, &  ! I1 needed for filenames
                          RXIJ,NXIJ,IXCP,RXCCLS, &
                          LMPIC, &
                          MYRANK,EMPIC,EMYRANK, &
                          JXCIJINT, &
                          naez, nxijd, &
                          lmpid, smpid, empid)
          endif
!=======================================================================
!=======================================================================

!=======================================================================
!     on the basis of new timings determine now new distribution of
!     work to 1 .. EMPID processors
!=======================================================================
          call EBALANCE('R',ITER,SCFSTEPS, &
                        IELAST,NPNT1, &
                        MYACTVRANK,ACTVCOMM, &
                        ETIME,EPROC,EPROCO, &
                        empid, iemxd)
!=======================================================================
!=======================================================================


!=======================================================================
!     in case of IGUESS and EMPID > 1 initial guess arrays might
!     have to be adjusted to new distributions
!=======================================================================
          if ((IGUESS==1).and.(EMPID>1)) then

            do ISPIN = 1,NSPIN

              if (SMPID==1) then
                MAPSPIN = 0
                PRSPIN   = ISPIN
              else
                MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
                PRSPIN   = 1
              endif

!       true beginning of SMPID-parallel section

              if(SRANK(SMPIB,SMPIC) == MAPSPIN) then

                call EPRDIST(IELAST,KMESH,NOFKS, &
                             PRSC(1,1,PRSPIN), &
                             SPRS(1,1,PRSPIN), &
                             CNVFAC(1,PRSPIN), &
                             MYRANK,EMPIC,EMYRANK, &
                             EPROC,EPROCO, &
                             lmpid, smpid, empid, naez, lmax, nguessd, ekmd, iemxd)

              endif
            enddo

          endif  ! IGUESS == 1 .and. EMPID > 1
!=======================================================================
!=======================================================================


!----------------------------------------------------------------------
! BEGIN only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
          if (LMPIC==1) then

            if (LLY==1) then
              ! get WEZRN and RNORM, the important input from previous
              ! calculations is LLY_GRDT_ALL
              ! TODO: all THETAS passed, but only 1 needed
              ! here atom processes communicate with each other
              call LLOYD0(EZ,WEZ,CLEB1C,DRDI,R,IRMIN,VINS,VISP, &
                          THETAS,ZAT,ICLEB1C, &
                          IFUNM,IPAN,IRCUT,LMSP,JEND,LOFLM1C, &
                          NTCELL,ICST, &
                          IELAST,IEND1,NAEZ,NSPIN,NSRA, &
                          WEZRN,RNORM, &
                          GMATN_ALL, &
                          LLY_GRDT_ALL, &
                          LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                          DMATLDAU, &
                          LMPIC,MYLRANK, &
                          LCOMM,LSIZE, &
                          lmpid*smpid*empid, lmax, irmd, irnsd, iemxd, &
                          irid, nfund, ncelld, ipand, ncleb)

! IME
              call OUTTIME(MYLRANK(1),'Lloyd processed......',TIME_I,ITER)
! IME
            else ! no Lloyd

              do IE=1,IELAST
                WEZRN(IE,1) = WEZ(IE)
                WEZRN(IE,2) = WEZ(IE)
              enddo
            endif

            ! now WEZRN stores the weights for E-integration

            call CINIT(IEMXD*(LMAXD+2)*NSPIND,DEN)
            DENEF = 0.0D0

            if (LDAU) then
              call CINIT(MMAXD*MMAXD*NSPIND*LMAXD1,DMATLDAU(1,1,1,1))
            endif

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN ==================================================================

            do ISPIN = 1,NSPIN
              ICELL = NTCELL(I1)
              IPOT = (I1-1) * NSPIN + ISPIN

              LDORHOEF = NPOL/=0  ! needed in RHOVAL
              call RHOVAL(LDORHOEF,ICST,IELAST, &
                          NSRA,ISPIN,NSPIN,EZ,WEZRN(1,ISPIN), &
                          DRDI(1,I1),R(1,I1),IRMIN(I1), &
                          VINS(IRMIND,1,ISPIN),VISP(1,ISPIN), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          THETAS(1,1,ICELL),IFUNM(1,ICELL),LMSP(1,ICELL), &
                          RHO2NS,R2NEF, &
                          DEN(0,1,ISPIN),ESPV(0,ISPIN), &
                          CLEB1C,LOFLM1C,ICLEB1C,IEND1,JEND, &
                          GMATN_ALL, &
                          LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                          DMATLDAU, &
                          iemxd, &
                          lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

              EBOT = E1

              call RHOCORE(EBOT,NSRA,ISPIN,NSPIN,I1, &  ! I1 is used only for debugging output
                           DRDI(1,I1),R(1,I1),VISP(1,ISPIN), &
                           A(I1),B(I1),ZAT(I1), &
                           IRCUT(0,I1),RHOCAT,QC, &
                           ECORE(1,ISPIN),NCORE(IPOT),LCORE(1,IPOT), &
                           irmd, ipand)

            end do

! SPIN ==================================================================
!      END do loop over spins
! SPIN ===================================================================

            if (LLY == 1) then
              call renormalizeDOS(DEN,RNORM,LMAXD1,IELAST,NSPIN,IEMXD)
            end if

            ! calculate DOS at Fermi level
            DENEF = calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIN)

            ! ---> l/m_s/atom-resolved charges, output -> CHARGE
            ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
            ! CHARGE -> written to result file, only for informative purposes
            call calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIN, WEZ, IEMXD)

! LDAU

            EULDAU = 0.0D0
            EDCLDAU = 0.0D0

            if (LDAU.and.NLDAU>=1) then

              call LDAUWMAT(I1,NSPIN,ITER,MIXING,DMATLDAU,NLDAU,LLDAU, &
                            ULDAU,JLDAU,UMLDAU,WMLDAU,EULDAU,EDCLDAU, &
                            lmaxd)

            endif

! LDAU

! ----------------------------------------------------------------------
! -->   determine total charge density expanded in spherical harmonics
! -------------------------------------------------------------- density

            call RHOTOTB_NEW(NSPIN,RHO2NS,RHOCAT, &
                         DRDI(:,I1),IRCUT(:,I1), &
                         LPOT,NFU(ICELL),LLMSP(1,ICELL),THETAS(:,:,ICELL),IPAN(I1), &
                         CATOM, &
                         irmd, irid, ipand, nfund)

            CHRGNT = CHRGNT + CATOM(1) - ZAT(I1)

            ! write to 'results1' - only to be read in in results.f - very unnecessary
            if (NPOL==0 .or. TEST('DOS     ')) then
              write(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN  ! write density of states (DEN) only when certain options set
            else
              write(71,rec=I1) QC,CATOM,CHARGE,ECORE
            end if

          endif
!----------------------------------------------------------------------
! END L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
        end if
      end do

!N ====================================================================
!     END do loop over atoms (NMPID-parallel)
!N ====================================================================

      close(66)
      close(71)

      call OUTTIME(MYLRANK(1),'density calculated ..',TIME_I,ITER)

!----------------------------------------------------------------------
! BEGIN L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
      if (LMPIC==1) then

!****************************************************** MPI COLLECT DATA

        call MPI_ALLREDUCE(CHRGNT,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        call DCOPY(1,WORK1,1,CHRGNT,1)

        call MPI_ALLREDUCE(DENEF,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        call DCOPY(1,WORK1,1,DENEF,1)

!****************************************************** MPI COLLECT DATA

! ----------------------------------------------------------------------

! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),0.03D0)*DSIGN(1.0D0,E2SHIFT) !FIXME: hardcoded maximal shift of 0.03
        EFOLD = E2

        if (ISHIFT < 2) E2 = E2 - E2SHIFT

        if( MYLRANK(LMPIC) == 0 ) then
          write (6,fmt=9020) EFOLD,E2SHIFT

! --> divided by NAEZ because the weight of each atom has been already
!     taken into account in 1c

          write (6,fmt=9030) E2,DENEF/DBLE(NAEZ)
          write(6,'(79(1H+),/)')
        end if

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIN)
! ----------------------------------------------------------------------

        open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
        form='unformatted')
        open (72,access='direct',recl=LRECRES2,file='results2', &
        form='unformatted')

! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================
        do I1 = 1,NAEZ
          if(MYLRANK(LMPIC) == MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then
            ICELL = NTCELL(I1)

            do ISPIN = 1,NSPIN

! -->     get correct density and valence band energies

              ESPV(0,ISPIN) = ESPV(0,ISPIN) - &
              EFOLD*CHRGNT/DBLE(NSPIN*NAEZ)
              if (LCORDENS) then

                do LM = 1,LMPOT
                  call DAXPY(IRC(I1),DF,R2NEF(1,LM,ISPIN),1, &
                  RHO2NS(1,LM,ISPIN),1)
                end do

              end if
! ----------------------------------------------------------------------
            end do

            call RHOMOM_NEW(CMOM,CMINST,LPOT,RHO2NS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ipand, ngshd)

            call OUTTIME(MYLRANK(1),'RHOMOM ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

            !call VINTRAS(LPOT,NSPIN,I1,RHO2NS,VONS, &
            !R,DRDI,IRCUT,IPAN,ICELL,ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            !THETAS,LMSP(1,ICELL), &
            !irmd, irid, nfund, ngshd, ipand)
            call VINTRAS_NEW(LPOT,NSPIN,RHO2NS,VONS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ngshd, ipand)

            call OUTTIME(MYLRANK(1),'VINTRAS ......',TIME_I,ITER)

            call STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM, &
            RBASIS,SMAT,VOLUME0,LASSLD,LMXSPD,NAEZD,I1)

            call OUTTIME(MYLRANK(1),'STRMAT ......',TIME_I,ITER)

            call VMADELBLK(CMOM,CMINST,LPOT,NSPIN, &
            NAEZ,VONS,ZAT,R,IRCUT,IPAN, &
            VMAD, &
            LMPOT,SMAT,CLEB,ICLEB,IEND, &
            LMXSPD,NCLEBD,LOFLM,DFAC,I1, &
            LMPIC,MYLRANK, &
            LCOMM,LSIZE, &
            irmd, ipand, lmpid*smpid*empid)

            call OUTTIME(MYLRANK(1),'VMADELBLK ......',TIME_I,ITER)

! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

            if (KFORCE==1 .and. ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCEH(CMOM,FLM,LPOT,I1,RHO2NS,VONS, &
              R,DRDI,IMT,ZAT,irmd)
              call FORCE(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R, &
              DRDI,IMT,naez,irmd)
! ---------------------------------------------------------------------
            end if

! Force Calculation stops here look after VXCDRV

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

            if (KTE==1) then
              call ESPCB(ESPC,NSPIN,I1,ECORE,LCORE,LCOREMAX,NCORE) !TODO

              call EPOTINB_NEW(EPOTIN,NSPIN,RHO2NS,VISP,R(:,I1),DRDI(:,I1), &
              IRMIN(I1),IRWS(I1),LPOT,VINS,IRCUT(:,I1),IPAN(I1),ZAT(I1), &
              irmd, irnsd, ipand)


              call ECOUB_NEW(CMOM,ECOU,LPOT,NSPIN,RHO2NS, &
              VONS,ZAT(I1),R(:,I1), &
              DRDI(:,I1),KVMAD,IRCUT(:,I1),IPAN(I1),IMAXSH,IFUNM(1,ICELL), &
              ILM,GSH,THETAS(:,:,ICELL),LMSP(1,ICELL), &
              irmd, irid, nfund, ipand, ngshd)

            end if
            call OUTTIME(MYLRANK(1),'KTE ......',TIME_I,ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
            call VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,I1,RHO2NS, &
            VONS,R,DRDI,A, &
            IRWS,IRCUT,IPAN,ICELL,GSH,ILM,IMAXSH,IFUNM(1,ICELL), &
            THETAS,LMSP(1,ICELL), &
            naez, irmd, irid, nfund, ngshd, ipand)

            call OUTTIME(MYLRANK(1),'VXCDRV ......',TIME_I,ITER)
! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

! Force calculation continues here

            if (KFORCE==1.and.ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCXC(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R, &
              ALAT,DRDI,IMT,ZAT, &
              LMPIC,MYLRANK, &
              LCOMM, &
              naez, irmd, lmpid*smpid*empid)
! ---------------------------------------------------------------------
            end if

            write(72,rec=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX, &
            EULDAU,EDCLDAU

! Force calculation ends
! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! =====================================================================
            VAV0 = 0.0D0
            VOL0 = 0.0D0

            call MTZERO(LMPOT,I1,NSPIN,VONS,ZAT,R,DRDI,IMT,IRCUT, &
                        IPAN,ICELL,LMSP(1,ICELL),IFUNM(1,ICELL), &
                        THETAS,IRWS,VAV0,VOL0, &
                        irmd, irid, nfund, ipand)

            call OUTTIME(MYLRANK(1),'MTZERO ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
          end if
        end do
        call OUTTIME(MYLRANK(1),'calculated pot ......',TIME_I,ITER)
! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================

! Calculate muffin-tin potential shift
!****************************************************** MPI COLLECT DATA
        WORK1(1) = VAV0
        WORK1(2) = VOL0
        call MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        VAV0 = WORK2(1)
        VOL0 = WORK2(2)
!****************************************************** MPI COLLECT DATA


        VBC(1) = -VAV0/VOL0
        if (ISHIFT>0) VBC(1) = VBC(1) + E2SHIFT
        VBC(2) = VBC(1)

        if(MYRANK==0) then
          write (6,fmt=9103) VOL0,VAV0,VBC(1)
          write(6,'(79(1H=),/)')
        end if

9103    format ('  VOL INT.',F16.9,'  VAV INT.',F16.9,'  VMT ZERO',F16.9)
! =====================================================================

! ---------------------------------------------------------------------

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

        do I1 = 1,NAEZ

          if(MYLRANK(LMPIC)== &
          MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

            ICELL = NTCELL(I1)
! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
            do ISPIN = 1,NSPIN
              IPOT = NSPIN* (I1-1) + ISPIN

              do IR = 1,IRCUT(IPAN(I1),I1)
                VONS(IR,1,ISPIN) = VONS(IR,1,ISPIN) + RFPI*VBC(ISPIN)
              end do

              call CONVOL(IRCUT(1,I1),IRC(I1),ICELL, &
                          IMAXSH(LMPOT),ILM,IFUNM(1,ICELL),LMPOT,GSH, &
                          THETAS,ZAT(I1),RFPI, &
                          R(1,I1),VONS(1,1,ISPIN),LMSP(1,ICELL), &
                          irid, nfund, irmd, ngshd)

            end do

! -->   final construction of the potentials (straight mixing)
            MIX = MIXING
            RMSAVQ = 0.0D0
            RMSAVM = 0.0D0

            call MIXSTR(RMSAVQ,RMSAVM,LPOT,LMPOT, &
            I1,NSPIN, &
            ITER,RFPI,FPI, &
            MIX, &
            FCM,IRC,IRMIN,R,DRDI,VONS, &
            VISP,VINS, &
            naez, irmd, irnsd)

            I1BRYD=I1
          end if
        end do

! -->  potential mixing procedures: Broyden or Andersen updating schemes
        if (IMIX>=3) then
          call BRYDBM(VISP,VONS,VINS, &
          LMPOT,R,DRDI,MIX, &
          IRC,IRMIN,NSPIN,I1BRYD, &
          IMIX,ITER, &
          UI2,VI2,WIT,SM1S,FM1S, &
          LMPIC,MYLRANK, &
          LCOMM, &
          itdbryd, irmd, irnsd, nspind, &
          LMPID * SMPID * EMPID)
        endif

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do I1 = 1,NAEZ
          if(MYLRANK(LMPIC)== &
          MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

            call resetPotentials(IRC(I1), IRMD, IRMIN(I1), IRMIND, LMPOTD, &
                                 NSPIN, VINS, VISP, VONS)

! ----------------------------------------------------- output_potential
            write(66,rec=I1) VINS,VISP,ECORE
! ----------------------------------------------------- output_potential
          end if
        end do

        close(72)
! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
! ====== write RMS convergency data - not parallelized, written by
! MYRANK=0   (RMSOUT) =================================================

        call RMSOUT(RMSAVQ,RMSAVM,ITER,E2,EFOLD, &
        SCFSTEPS,VBC,QBOUND,NSPIN,NAEZ, &
        KXC,LPOT,A,B,IRC, &
        VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT, &
        ECORE,LCORE,NCORE,ZAT,ITITLE, &
        LMPIC,MYLRANK, &
        LCOMM,LSIZE, &
        irmd, irnsd, lmpid*smpid*empid)

9020    format ('                old', &
        ' E FERMI ',F12.6,' Delta E_F = ',f12.6)
9030    format ('                new', &
        ' E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)



! Wait here in order to guarantee regular and non-errorneous output
! in RESULTS

        call MPI_BARRIER(LCOMM(LMPIC),IERR)


! -----------------------------------------------------------------
! L-MPI: only process with MYLRANK(LMPIC = 1) = 0 is working here
! -----------------------------------------------------------------
        if(MYLRANK(LMPIC)==0) then

          ! DOS was written to file 'results1' and read out here just
          ! to be written in routine wrldos
          ! also other stuff is read from results1
          call RESULTS(LRECRES2,IELAST,ITER,LMAX,NAEZ,NPOL, &
          NSPIN,KPRE,KTE,LPOT,E1,E2,TK,EFERMI, &
          ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ,LDAU, &
          iemxd)

! --> update energy contour

          ! E2 is used twice as a parameter - dangerous!
          call EMESHT(EZ,DEZ,IELAST,E1,E2,E2,TK, &
          NPOL,NPNT1,NPNT2,NPNT3,IEMXD)

          do IE = 1,IELAST
            WEZ(IE) = -2.D0/PI*DEZ(IE)
          end do

          write(6,'(79(1H=))')

! .. get info on MYACTVRANK of this processor: to be used in
!    subsequent reduce-commands
          MYBCRANK = MYACTVRANK
! ..
        endif
! -----------------------------------------------------------------
! L-MPI: only process with MYLRANK(LMPIC = 1) = 0 is working here
! -----------------------------------------------------------------

! ..
      endif
! -----------------------------------------------------------------
! L-MPI: only processes with LMPIC = 1 are working here
! -----------------------------------------------------------------
      close(66)  ! close 'vpotnew'
! -----------------------------------------------------------------

!      CALL MPI_BARRIER(ACTVCOMM,IERR)

      ! why? all processes except 1 have MYBCRANK = 0, this allreduce
      ! tells all the other processes who is the root
      ! not really necessary
      call MPI_ALLREDUCE(MYBCRANK,BCRANK,1,MPI_INTEGER,MPI_MAX, &
      ACTVCOMM,IERR)

      call MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      ACTVCOMM,IERR)

!      CALL MPI_BARRIER(ACTVCOMM,IERR)

      if(MYLRANK(1)==0) then

        open (67,file='energy_mesh',form='unformatted')
        write (67) IELAST,EZ,WEZ,E1,E2
        write (67) NPOL,TK,NPNT1,NPNT2,NPNT3
        close (67)

        write(6,'(79(1H=))')
        write(6,'(19X,A,I3,A,I10)') '       ITERATION : ', &
        ITER,' SUM of QMR ',NOITER_ALL
        write(6,'(79(1H=),/)')
        call SYSTEM_CLOCK(SYSTEM_F,RATETIME,MAXTIME)
        WALLCLOCK_F=MPI_WTIME()
        write(6,'(79(1H=))')
        write(6,*) 'Wallclock, System and CPU-time compared:'
        write(6,*) 'MPI_WTIME=',(WALLCLOCK_F-WALLCLOCK_I),'sec'
        write(6,*) 'SYSTEM_CLOCK=',(SYSTEM_F-SYSTEM_I)/100,'sec'
        call OUTTIME(MYLRANK(1),'end .................', &
        TIME_I,ITER)
        write(6,'(79(1H=))')
        write(2,'(79(1H=))')
        call OUTTIME(MYLRANK(1),'finished in .........', &
        TIME_S,ITER)
        write(2,'(79(1H=))')
        write(6,'(79(1H=),/)')
      endif

! manual exit possible by creation of file 'STOP' in home directory

      inquire(file='STOP',exist=STOPIT)
      if (STOPIT) goto 200


! ######################################################################
! ######################################################################
    enddo          ! SC ITERATION LOOP ITER=1, SCFSTEPS
! ######################################################################
! ######################################################################


200 continue

    if (MYLRANK(1)==0) close(2)    ! TIME

    if (KFORCE==1) close(54)
! ======================================================================
! ======================================================================

! Free communicators and groups ..
! ..
    if (MYLRANK(LMPIC)>=0) then
      call MPI_COMM_FREE(LCOMM(LMPIC),IERR)
    endif

    call MPI_GROUP_FREE(LGROUP(LMPIC),IERR)

    call MPI_COMM_FREE(ACTVCOMM,IERR)
    call MPI_GROUP_FREE(ACTVGROUP,IERR)
! .. .

  endif ! ACTVGROUP

!=====================================================================
!     processors not fitting in NAEZ*LMPID do nothing ...
! ... and wait here
!=====================================================================


!      WRITE(6,*) 'BARRIER i:',MYRANK
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!      WRITE(6,*) 'BARRIER f:',MYRANK
  call MPI_FINALIZE(IERR)

!-----------------------------------------------------------------------------
! Array DEallocations BEGIN
!-----------------------------------------------------------------------------
  deallocate(EZ, stat = memory_stat)
  deallocate(WEZ, stat = memory_stat)
  deallocate(DEZ, stat = memory_stat)
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

!-----------------------------------------------------------------------------
! Array DEallocations END
!-----------------------------------------------------------------------------

end program MAIN2
