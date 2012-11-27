module main2_arrays_mod
  implicit none

  SAVE

  ! KKRnano-Arrays
  double complex, dimension(:,:,:), allocatable ::  DEN
  double complex, dimension(:,:,:), allocatable ::  DSYMLL
  double precision, dimension(:,:), allocatable :: RBASIS

  double precision, dimension(:,:), allocatable :: SMAT
  double precision, dimension(:,:), allocatable :: RNORM
  double precision, dimension(:,:,:), allocatable :: BZKP
  double precision, dimension(:,:), allocatable :: VOLCUB
  double precision, dimension(:), allocatable :: VOLBZ
  double precision, dimension(:,:), allocatable :: RR

  integer, dimension(:), allocatable :: KMESH
  integer, dimension(:), allocatable :: NOFKS
  integer, dimension(:,:), allocatable :: EZOA

  integer, dimension(:), allocatable :: NUMN0
  integer, dimension(:,:), allocatable :: INDN0

  integer, dimension(:), allocatable :: LLDAU    ! LDA+U

  double complex, dimension(:,:,:), allocatable ::  TMATN
  double complex, dimension(:,:,:), allocatable ::  DTDE
  double complex, dimension(:,:,:), allocatable ::  TREFLL
  double complex, dimension(:,:,:), allocatable ::  DTREFLL
  double complex, dimension(:,:,:,:), allocatable ::  DGREFN
  double complex, dimension(:,:,:,:), allocatable ::  GREFN
  double complex, dimension(:,:,:,:), allocatable ::  GMATN

  !----- Initial guess ---------------------------------------------------
  complex, dimension(:,:,:), allocatable :: PRSC

  !----- Lloyd -----------------------------------------------------------
  double complex, dimension(:,:), allocatable ::  LLY_G0TR
  double complex, dimension(:,:), allocatable ::  LLY_GRDT
  double complex, dimension(:), allocatable ::  TR_ALPH
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  double precision, dimension(:), allocatable :: ECOU
  double precision, dimension(:,:), allocatable :: ESPC
  double precision, dimension(:,:), allocatable :: ESPV
  double precision, dimension(:), allocatable :: EXC

  double precision, dimension(:), allocatable :: ZAT

  ! ----------------------------------------------------------------------
  double precision, dimension(:,:,:), allocatable ::  R2NEF
  double precision, dimension(:,:,:), allocatable ::  RHO2NS
  double precision, dimension(:,:), allocatable ::  CHARGE
  !     ..

  double precision, dimension(:), allocatable :: CMINST    ! charge moment of interstitial
  double precision, dimension(:), allocatable :: CMOM      ! LM moment of total charge
  double precision, dimension(:), allocatable :: CATOM     ! total charge per atom

  !     .. FORCES
  !double precision, dimension(:,:), allocatable :: FLM
  !double precision, dimension(:,:), allocatable :: FLMC
  !     .. MIXING
  double precision, dimension(:), allocatable :: SM1S
  double precision, dimension(:), allocatable :: FM1S
  double precision, dimension(:,:), allocatable :: UI2
  double precision, dimension(:,:), allocatable :: VI2
  double precision, dimension(:), allocatable :: WIT

  double precision, dimension(:,:,:), allocatable :: RCLS
  double precision, dimension(:), allocatable :: RMTREF
  double precision, dimension(:), allocatable :: VREF

  integer, dimension(:,:), allocatable :: ATOM
  integer, dimension(:), allocatable :: CLS
  integer, dimension(:), allocatable :: NACLS
  integer, dimension(:), allocatable :: REFPOT

  !==============================================================================
  ! Dimension Parameters.
  !==============================================================================

  integer :: NAEZ
  integer :: LMAXD
  integer :: NREFD
  integer :: IRID
  integer :: BCPD
  integer :: NACLSD
  integer :: IRMD
  integer :: IEMXD
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
  integer::   LMAXD1
  integer::   MMAXD
  integer::   LMXSPD
  integer::   LMPOTD
  integer::   IRMIND
  integer::   LRECRES2
  integer::   NTIRD    ! for Broyden mixing

  !derived parameters
  integer :: LPOT
  integer :: NGUESSD

  !Parallelisation
  integer :: SMPID
  integer :: EMPID
  integer :: NTHRDS


CONTAINS

  !---------------------------------------------------------------------
  !> Reads array dimension parameters (previously in inc.p and inc.cls)
  !> from file inp0.unf
  !---------------------------------------------------------------------
  subroutine read_dimension_parameters()

    implicit none

    integer :: FILEHANDLE = 67

    open (FILEHANDLE, FILE='inp0.unf', FORM='unformatted')

    read(FILEHANDLE) LMAXD
    read(FILEHANDLE) NSPIND
    read(FILEHANDLE) NAEZ
    read(FILEHANDLE) IRNSD
    read(FILEHANDLE) IRMD
    read(FILEHANDLE) NREFD
    read(FILEHANDLE) NRD
    read(FILEHANDLE) IRID
    read(FILEHANDLE) NFUND
    read(FILEHANDLE) NCELLD
    read(FILEHANDLE) NACLSD
    read(FILEHANDLE) NCLSD
    read(FILEHANDLE) IPAND
    read(FILEHANDLE) NXIJD
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

    ! derived parameter
    LPOT = 2*LMAXD

  end subroutine read_dimension_parameters

  !----------------------------------------------------------------------------
  subroutine allocate_main2_arrays()
    use ErrorMessages_mod
    implicit none

    integer :: memory_stat

    allocate(DEN(0:LMAXD1,IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(DSYMLL(LMMAXD,LMMAXD,48), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RBASIS(3,NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(SMAT(LMXSPD,NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RNORM(IEMXD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(BZKP(3,KPOIBZ,MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VOLCUB(KPOIBZ,MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VOLBZ(MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RR(3,0:NRD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(KMESH(IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NOFKS(MAXMSHD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EZOA(NACLSD,NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NUMN0(NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(INDN0(NAEZ,NACLSD), stat = memory_stat)
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
    allocate(PRSC(NGUESSD*LMMAXD,EKMD,NSPIND-SMPID+1), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLY_G0TR(NCLSD,IEMXD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(LLY_GRDT(IEMXD,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(TR_ALPH(NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ECOU(0:LPOT), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ESPC(0:3,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ESPV(0:LMAXD1,NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(EXC(0:LPOT), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ZAT(NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(R2NEF(IRMD,LMPOTD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RHO2NS(IRMD,LMPOTD,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CHARGE(0:LMAXD1,2), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CMINST(LMPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CMOM(LMPOTD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CATOM(NSPIND), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    !allocate(FLM(-1:1,NAEZ), stat = memory_stat)
    !if(memory_stat /= 0) call fatalMemoryError("main2")
    !allocate(FLMC(-1:1,NAEZ), stat = memory_stat)
    !if(memory_stat /= 0) call fatalMemoryError("main2")
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
    allocate(RCLS(3,NACLSD,NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(RMTREF(NREFD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(VREF(NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(ATOM(NACLSD,NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(CLS(NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(NACLS(NCLSD), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(REFPOT(NAEZ), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("main2")

    !initialise to be safe
    RHO2NS = 0.0d0
    R2NEF = 0.0d0
    !DTIXIJ = dcmplx(0.0d0, 0.0d0)

    ! use garbage values
    GMATN = dcmplx(99999.0d0, 99999.0d0)
    LLY_GRDT = dcmplx(99999.0d0, 99999.0d0)

  end subroutine

  !------------------------------------------------------------------------------
  subroutine deallocate_main2_arrays()
    implicit none

    integer :: memory_stat

    deallocate(DEN, stat = memory_stat)
    deallocate(DSYMLL, stat = memory_stat)
    deallocate(RBASIS, stat = memory_stat)
    deallocate(SMAT, stat = memory_stat)
    deallocate(RNORM, stat = memory_stat)
    deallocate(BZKP, stat = memory_stat)
    deallocate(VOLCUB, stat = memory_stat)
    deallocate(VOLBZ, stat = memory_stat)
    deallocate(RR, stat = memory_stat)
    deallocate(KMESH, stat = memory_stat)
    deallocate(NOFKS, stat = memory_stat)
    deallocate(EZOA, stat = memory_stat)
    deallocate(NUMN0, stat = memory_stat)
    deallocate(INDN0, stat = memory_stat)
    deallocate(LLDAU, stat = memory_stat)
    deallocate(TMATN, stat = memory_stat)
    deallocate(DTDE, stat = memory_stat)
    deallocate(TREFLL, stat = memory_stat)
    deallocate(DTREFLL, stat = memory_stat)
    deallocate(DGREFN, stat = memory_stat)
    deallocate(GREFN, stat = memory_stat)
    deallocate(GMATN, stat = memory_stat)
    deallocate(PRSC, stat = memory_stat)
    deallocate(LLY_G0TR, stat = memory_stat)
    deallocate(LLY_GRDT, stat = memory_stat)
    deallocate(TR_ALPH, stat = memory_stat)
    deallocate(ECOU, stat = memory_stat)
    deallocate(ESPC, stat = memory_stat)
    deallocate(ESPV, stat = memory_stat)
    deallocate(EXC, stat = memory_stat)
    deallocate(ZAT, stat = memory_stat)
    deallocate(R2NEF, stat = memory_stat)
    deallocate(RHO2NS, stat = memory_stat)
    deallocate(CHARGE, stat = memory_stat)
    deallocate(CMINST, stat = memory_stat)
    deallocate(CMOM, stat = memory_stat)
    deallocate(CATOM, stat = memory_stat)
    !deallocate(FLM, stat = memory_stat)
    !deallocate(FLMC, stat = memory_stat)
    deallocate(SM1S, stat = memory_stat)
    deallocate(FM1S, stat = memory_stat)
    deallocate(UI2, stat = memory_stat)
    deallocate(VI2, stat = memory_stat)
    deallocate(WIT, stat = memory_stat)
    deallocate(RCLS, stat = memory_stat)
    deallocate(RMTREF, stat = memory_stat)
    deallocate(VREF, stat = memory_stat)
    deallocate(ATOM, stat = memory_stat)
    deallocate(CLS, stat = memory_stat)
    deallocate(NACLS, stat = memory_stat)
    deallocate(REFPOT, stat = memory_stat)

  end subroutine

end module main2_arrays_mod
