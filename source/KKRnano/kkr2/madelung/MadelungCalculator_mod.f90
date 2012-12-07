module MadelungCalculator_mod
  implicit none

  public :: MadelungCalculator

  private :: MadelungLatticeData
  private :: MadelungHarmonics
  private :: MadelungClebschData

  type MadelungLatticeData
    integer NGMAX
    integer NRMAX
    !INTEGER NSG(ISHLD)
    !integer NSR(ISHLD)
    integer, allocatable, dimension(:) :: NSG
    integer, allocatable, dimension(:) :: NSR
    integer NSHLG
    integer NSHLR
    !double precision GN(3,NMAXD)
    !double precision RM(3,NMAXD)
    double precision, allocatable, dimension(:,:) :: GN
    double precision, allocatable, dimension(:,:) :: RM
    integer NMAXD
    integer ISHLD
  end type

  type MadelungHarmonics
    double precision, dimension(:), allocatable :: WG
    double precision, dimension(:,:,:), allocatable :: YRG
    integer :: LASSLD
  end type

  type MadelungClebschData
    double precision, dimension(:), allocatable :: CLEB
    integer, dimension(:,:), allocatable :: ICLEB
    integer, dimension(:), allocatable :: LOFLM
    integer :: IEND
    integer :: NCLEBD
  end type

  type MadelungCalculator
    private

    double precision ALAT
    integer LMXSPD
    integer LPOT
    integer LMPOTD
    integer LASSLD
    double precision VOLUME0

    double precision, dimension(:,:), allocatable :: DFAC

    type (MadelungLatticeData) :: mlattice
    type (MadelungClebschData) :: clebsch

  end type

  contains

  !----------------------------------------------------------------------------
  !> Creates a MadelungCalculator object.
  !> Don't forget to free resources with destroyMadelungCalculator
  subroutine createMadelungCalculator(madelung_calc, lmax, ALAT, RMAX, GMAX, &
                                      BRAVAIS, NMAXD, ISHLD)
    implicit none
    type (MadelungCalculator), intent(inout) :: madelung_calc
    integer, intent(in) :: lmax
    double precision, intent(in) :: ALAT
    double precision, intent(in) :: RMAX
    double precision, intent(in) :: GMAX
    double precision, intent(in)::  BRAVAIS(3,3)

    integer, intent(in) :: NMAXD
    integer, intent(in) :: ISHLD
    !---------------------------------------------------------

    type (MadelungHarmonics) :: harmonics

    double precision ::  RECBV(3,3)

    integer :: LPOT
    integer :: LMXSPD
    integer :: LMPOTD
    integer :: LASSLD
    integer :: memory_stat

    integer, parameter :: no_output = 1

    LPOT = 2*lmax
    LASSLD = 4*lmax
    LMXSPD = (2*LPOT+1)**2
    LMPOTD = (LPOT+1)**2

    madelung_calc%ALAT = ALAT
    madelung_calc%LPOT = LPOT
    madelung_calc%LMXSPD = LMXSPD
    madelung_calc%LMPOTD = LMPOTD
    madelung_calc%LASSLD = LASSLD

    allocate(madelung_calc%DFAC(0:LPOT,0:LPOT), stat = memory_stat)
    if(memory_stat /= 0) then
      write(*,*) "MadelungCalculator: Allocation error, DFAC"
      stop
    end if

    call LATTIX99(ALAT,BRAVAIS,RECBV,madelung_calc%VOLUME0, .false.)

    call createMadelungHarmonics(harmonics, lmax)

    call createMadelungLatticeData(madelung_calc%mlattice, NMAXD, ISHLD)

    call createMadelungClebschData(madelung_calc%clebsch, LMXSPD, LMPOTD)

    ! now wrap nasty function call
    call MADELUNG3D(LPOT,harmonics%YRG,harmonics%WG, &
     ALAT, &
     RMAX,GMAX,BRAVAIS,RECBV, &
     LMXSPD,LASSLD,LPOT,LMPOTD, &
     NMAXD,ISHLD, &
     LMPOTD,madelung_calc%clebsch%CLEB,madelung_calc%clebsch%ICLEB,madelung_calc%clebsch%IEND, &
     madelung_calc%clebsch%NCLEBD,madelung_calc%clebsch%LOFLM, &
     madelung_calc%DFAC, &
     madelung_calc%mlattice%NGMAX,madelung_calc%mlattice%NRMAX, &
     madelung_calc%mlattice%NSG, madelung_calc%mlattice%NSR, &
     madelung_calc%mlattice%NSHLG, madelung_calc%mlattice%NSHLR, &
     madelung_calc%mlattice%GN, madelung_calc%mlattice%RM, &
     no_output)

    call destroyMadelungHarmonics(harmonics)

  end subroutine

  !----------------------------------------------------------------------------
  !> Destroys a MadelungCalculator object.
  subroutine destroyMadelungCalculator(madelung_calc)
    implicit none
    type (MadelungCalculator), intent(inout) :: madelung_calc
    !--------------------------------------------------------

    call destroyMadelungLatticeData(madelung_calc%mlattice)

    call destroyMadelungClebschData(madelung_calc%clebsch)

    deallocate(madelung_calc%DFAC)

  end subroutine

  !----------------------------------------------------------------------------
  !> Calculates Lattice sums needed for Madelung potential calculation.
  !> Needs a properly constructed MadelungCalculator object.
  !> STRMAT wrapper
  subroutine calculateMadelungLatticeSum(madelung_calc, naez, atom_index, rbasis, smat)
    implicit none
    type (MadelungCalculator), intent(inout) :: madelung_calc
    integer, intent(in) :: naez
    integer, intent(in) :: atom_index
    double precision, dimension(3,NAEZ), intent(in) :: rbasis
    double precision, dimension(madelung_calc%LMXSPD,naez), intent(out) :: SMAT   !SMAT(LMXSPD,NAEZ)
    !---------------------------------------------------------


    call STRMAT(madelung_calc%ALAT,madelung_calc%LPOT,NAEZ,madelung_calc%mlattice%NGMAX, &
    madelung_calc%mlattice%NRMAX,madelung_calc%mlattice%NSG,madelung_calc%mlattice%NSR, &
    madelung_calc%mlattice%NSHLG,madelung_calc%mlattice%NSHLR, &
    madelung_calc%mlattice%GN,madelung_calc%mlattice%RM, RBASIS,SMAT, &
    madelung_calc%VOLUME0, &
    madelung_calc%LASSLD,madelung_calc%LMXSPD,naez,atom_index)

  end subroutine

  !----------------------------------------------------------------------------
  !> Add Madelung potential to VONS.
  !> Needs SMAT (Lattice sums from calculateMadelungLatticeSum)
  !> principal input: CMOM, CMINST, SMAT, VONS --> VONS (changed)
  !> Wrapper for VMADELBLK
  subroutine addMadelungPotential_com(madelung_calc, CMOM, CMINST, NSPIN, &
                                      NAEZ, VONS, ZAT, R, IRCUT, IPAN, VMAD, &
                                      SMAT, rank, communicator, comm_size, irmd, ipand)

    implicit none

    type (MadelungCalculator), intent(in) :: madelung_calc
    double precision, intent(inout) :: CMOM(madelung_calc%LMPOTD)
    double precision, intent(inout) :: CMINST(madelung_calc%LMPOTD)
    integer, intent(in) :: nspin
    integer, intent(in) :: naez
    double precision, intent(inout) :: VONS(IRMD,madelung_calc%LMPOTD,2)

    double precision, intent(in) ::  R(IRMD)
    double precision, intent(in) ::  ZAT(naez)
    integer IRCUT(0:IPAND)
    integer IPAN
    double precision VMAD
    double precision, intent(in) ::  SMAT(madelung_calc%LMXSPD,naez)

    integer, intent(in) :: rank
    integer, intent(in) :: communicator
    integer, intent(in) :: comm_size
    integer, intent(in) :: irmd
    integer, intent(in) :: ipand


    call VMADELBLK_new_com(CMOM,CMINST,madelung_calc%LPOT,NSPIN,NAEZ, &
    VONS,ZAT,R,IRCUT,IPAN,VMAD, &
    madelung_calc%LMPOTD, SMAT, &
    madelung_calc%clebsch%CLEB,madelung_calc%clebsch%ICLEB,madelung_calc%clebsch%IEND, &
    madelung_calc%LMXSPD,madelung_calc%clebsch%NCLEBD,madelung_calc%clebsch%LOFLM, &
    madelung_calc%DFAC, &
    rank, communicator, comm_size, irmd, ipand)

  !      SUBROUTINE VMADELBLK_new_com(CMOM,CMINST,LPOT,NSPIN,NAEZ,
  !     &                     VONS,ZAT,R,
  !     &                     IRCUT,IPAN,VMAD,
  !     &                     LMPOT,SMAT,CLEB,ICLEB,IEND,
  !     &                     LMXSPD,NCLEBD,LOFLM,DFAC,
  !     >                     MYLRANK,
  !     >                     communicator,comm_size,
  !     &                     irmd, ipand)

  end subroutine



!============= Helper routines =============================================

  !----------------------------------------------------------------------------
  subroutine createMadelungHarmonics(harmonics, lmax)
    implicit none
    type (MadelungHarmonics), intent(inout) :: harmonics
    integer, intent(in) :: lmax
    !--------------------------------------------------------------------------

    !integer :: memory_stat
    integer :: LASSLD

    LASSLD = 4*lmax
    harmonics%LASSLD = LASSLD

    allocate(harmonics%WG(LASSLD))
    !if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(harmonics%YRG(LASSLD,0:LASSLD,0:LASSLD))
    !if(memory_stat /= 0) call fatalMemoryError("main2")

    call GAUNT2(harmonics%WG,harmonics%YRG,LMAX)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyMadelungHarmonics(harmonics)
    implicit none
    type (MadelungHarmonics), intent(inout) :: harmonics

    deallocate(harmonics%WG)
    deallocate(harmonics%YRG)
  end subroutine

  !----------------------------------------------------------------------------
  subroutine createMadelungLatticeData(mlattice, NMAXD, ISHLD)
    implicit none
    type (MadelungLatticeData), intent(inout) :: mlattice
    integer, intent(in) :: NMAXD
    integer, intent(in) :: ISHLD
    !--------------------------------------------------------------------------

    mlattice%NMAXD = NMAXD
    mlattice%ISHLD = ISHLD

    allocate(mlattice%GN(3,NMAXD))
    !if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(mlattice%RM(3,NMAXD))
    !if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(mlattice%NSG(ISHLD))
    !if(memory_stat /= 0) call fatalMemoryError("main2")
    allocate(mlattice%NSR(ISHLD))
    !if(memory_stat /= 0) call fatalMemoryError("main2")

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyMadelungLatticeData(mlattice)
    implicit none
    type (MadelungLatticeData), intent(inout) :: mlattice
    !--------------------------------------------------------------------------

    deallocate(mlattice%GN)

    deallocate(mlattice%RM)

    deallocate(mlattice%NSG)

    deallocate(mlattice%NSR)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine createMadelungClebschData(clebsch, LMXSPD, LMPOTD)
    implicit none
    type (MadelungClebschData), intent(inout) :: clebsch
    integer, intent(in) :: LMXSPD
    integer, intent(in) :: LMPOTD
    !--------------------------------------------------------------------------

    allocate(clebsch%CLEB(LMXSPD*LMPOTD))
    allocate(clebsch%LOFLM(LMXSPD))
    allocate(clebsch%ICLEB(LMXSPD*LMPOTD,3))
  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyMadelungClebschData(clebsch)
    implicit none
    type (MadelungClebschData), intent(inout) :: clebsch
    !--------------------------------------------------------------------------

    deallocate(clebsch%CLEB)
    deallocate(clebsch%LOFLM)
    deallocate(clebsch%ICLEB)

  end subroutine

end module MadelungCalculator_mod
