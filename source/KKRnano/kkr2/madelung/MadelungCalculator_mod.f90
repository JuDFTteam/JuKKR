!> Instructions:
!> Create Madelung calculator
!> Create Madelung lattice sum(s) and configure with Madelung calculator
!> use...
!> destroy Madelung lattice sum(s)
!> destroy Madelung calculator

module MadelungCalculator_mod
  implicit none
  private

  public :: MadelungCalculator, MadelungClebschData, MadelungLatticeSum, create, destroy
  public :: calculateMadelungLatticeSum, calc_dfac
  public :: createMadelungLatticeSum, destroyMadelungLatticeSum ! deprecated
  public :: createMadelungCalculator, destroyMadelungCalculator ! deprecated

  !----------------------------------------------------------------------------
  type MadelungLatticeData
    integer NGMAX
    integer NRMAX
    integer, allocatable :: NSG(:)
    integer, allocatable :: NSR(:)
    integer NSHLG
    integer NSHLR
    double precision, allocatable :: GN(:,:)
    double precision, allocatable :: RM(:,:)
    integer NMAXD
    integer ISHLD
  endtype

  !----------------------------------------------------------------------------
  type MadelungHarmonics
    double precision, allocatable :: WG(:)
    double precision, allocatable :: YRG(:,:,:)
    integer :: LASSLD
  endtype

  !----------------------------------------------------------------------------
  type MadelungClebschData
    double precision, allocatable :: CLEB(:)
    integer, allocatable :: ICLEB(:,:)
    integer, allocatable :: LOFLM(:)
    integer :: IEND
    integer :: NCLEBD
  endtype

  !----------------------------------------------------------------------------
  type MadelungCalculator
    !private

    double precision :: ALAT
    integer :: LMXSPD
    integer :: LPOT
    integer :: LMPOTD
    integer :: LASSLD
    double precision :: VOLUME0

    double precision, allocatable :: DFAC(:,:)

    type(MadelungLatticeData) :: mlattice
    type(MadelungClebschData) :: clebsch

  endtype

  !----------------------------------------------------------------------------
  type MadelungLatticeSum
    double precision, allocatable :: SMAT(:,:)
    type(MadelungCalculator), pointer :: madelung_calc
    integer :: num_atoms
  endtype

  
  interface create
    module procedure createMadelungCalculator, createMadelungClebschData, createMadelungLatticeSum
  endinterface
  
  interface destroy
    module procedure destroyMadelungCalculator, destroyMadelungClebschData, destroyMadelungLatticeSum
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  !> Creates a MadelungCalculator object.
  !> Don't forget to free resources with destroyMadelungCalculator
  subroutine createMadelungCalculator(madelung_calc, lmax, ALAT, RMAX, GMAX, BRAVAIS, NMAXD, ISHLD)
    use Lattice_mod, only: lattix99
    type(MadelungCalculator), intent(inout) :: madelung_calc
    integer, intent(in) :: lmax
    double precision, intent(in) :: ALAT
    double precision, intent(in) :: RMAX
    double precision, intent(in) :: GMAX
    double precision, intent(in)::  BRAVAIS(3,3)

    integer, intent(in) :: NMAXD
    integer, intent(in) :: ISHLD
    !---------------------------------------------------------

    type(MadelungHarmonics) :: harmonics

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
    endif

    call lattix99(alat, bravais, recbv, madelung_calc%volume0, .false.)

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

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Destroys a MadelungCalculator object.
  subroutine destroyMadelungCalculator(madelung_calc)
    type(MadelungCalculator), intent(inout) :: madelung_calc
    !--------------------------------------------------------

    call destroyMadelungLatticeData(madelung_calc%mlattice)

    call destroyMadelungClebschData(madelung_calc%clebsch)

    deallocate(madelung_calc%DFAC)

  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Creates data storage for a Madelung lattice sum.
  !>
  !> Has to be configured with a properly setup MadelungCalculator.
  !> The MadelungCalculator can be reused for several MadelungLatticeSums
  !> as long as the geometry and lmax does not change
  subroutine createMadelungLatticeSum(madelung_sum, madelung_calc, num_atoms)
    type(MadelungLatticeSum), intent(inout) :: madelung_sum
    type(MadelungCalculator), target, intent(in)    :: madelung_calc
    integer, intent(in) :: num_atoms

    !---- local
    type(MadelungCalculator), pointer    :: mad_ptr

    mad_ptr => madelung_calc
    madelung_sum%madelung_calc => mad_ptr

    madelung_sum%num_atoms = num_atoms

    allocate(madelung_sum%SMAT(madelung_calc%LMXSPD, num_atoms))

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Destroys and frees storage of Madelung lattice sum.
  subroutine destroyMadelungLatticeSum(madelung_sum)
    type(MadelungLatticeSum), intent(inout) :: madelung_sum

    deallocate(madelung_sum%SMAT)

  endsubroutine ! destroy  

  !----------------------------------------------------------------------------
  !> Calculates Lattice sum for atom 'atom_index' for positions 'rbasis'
  !>  needed for Madelung potential calculation.
  !>
  !> Needs a properly constructed MadelungLatticeSum object.
  !> STRMAT wrapper
  subroutine calculateMadelungLatticeSum(madelung_sum, atom_index, rbasis)
    type(MadelungLatticeSum), intent(inout) :: madelung_sum
    integer, intent(in) :: atom_index
    double precision, intent(in) :: rbasis(3,madelung_sum%num_atoms)
    !---------------------------------------------------------
    type(MadelungCalculator), pointer :: madelung_calc
    integer :: naez

    madelung_calc => madelung_sum%madelung_calc
    naez = madelung_sum%num_atoms

    call STRMAT(madelung_calc%ALAT,madelung_calc%LPOT,NAEZ,madelung_calc%mlattice%NGMAX, &
    madelung_calc%mlattice%NRMAX,madelung_calc%mlattice%NSG,madelung_calc%mlattice%NSR, &
    madelung_calc%mlattice%NSHLG,madelung_calc%mlattice%NSHLR, &
    madelung_calc%mlattice%GN,madelung_calc%mlattice%RM, RBASIS, madelung_sum%SMAT, &
    madelung_calc%VOLUME0, &
    madelung_calc%LASSLD,madelung_calc%LMXSPD,naez,atom_index)

  endsubroutine ! calc

  !----------------------------------------------------------------------------
  !
  ! --> calculate:                               (2*(l+l')-1)!!
  !                 dfac(l,l') = (4pi)**2 *  ----------------------
  !                                          (2*l+1)!! * (2*l'+1)!!
  subroutine calc_dfac(dfac, lpot)
    double precision, intent(inout) :: dfac(0:lpot, 0:lpot)
    integer, intent(in) :: lpot

    double precision :: pi, fpi
    integer :: L1, L2
    pi = 4.0d0*atan(1.0d0)
    fpi = 4.0d0*pi

    dfac(0,0) = fpi*fpi
    do l1 = 1,lpot
       dfac(l1,0) = dfac(l1-1,0)*dble(2*l1-1)/dble(2*l1+1)
       dfac(0,l1) = dfac(l1,0)
       do l2 = 1,l1
          dfac(l1,l2) = dfac(l1,l2-1)*dble(2*(l1+l2)-1)/dble(2*l2+1)
          dfac(l2,l1) = dfac(l1,l2)
       enddo
    enddo
  endsubroutine ! calc

  !----------------------------------------------------------------------------
  subroutine initMadelungClebschData(clebsch, lmax)
    type(MadelungClebschData), intent(inout) :: clebsch
    integer, intent(in) :: lmax

    integer lpot, lassld, lmxspd, lmpotd, nclebd, L, M, I
    type(MadelungHarmonics) :: harmonics

    LPOT = 2*lmax
    LASSLD = 4*lmax
    LMXSPD = (2*LPOT+1)**2
    LMPOTD = (LPOT+1)**2
    nclebd = lmxspd*lmpotd

    i = 1
    do l = 0,2*lpot
       do m = -l,l
          clebsch%loflm(i) = l
          i = i + 1
       enddo
    enddo

    call createMadelungHarmonics(harmonics, lmax)

    call madelgaunt(lpot,harmonics%yrg,harmonics%wg,clebsch%cleb,clebsch%icleb,clebsch%iend,lassld,nclebd)
    clebsch%nclebd = nclebd

    call destroyMadelungHarmonics(harmonics)
  endsubroutine ! init

!============= Helper routines =============================================

  !----------------------------------------------------------------------------
  subroutine createMadelungHarmonics(harmonics, lmax)
    use Harmonics_mod, only: Gaunt2 ! initialization of wg and yrg
    type(MadelungHarmonics), intent(inout) :: harmonics
    integer, intent(in) :: lmax
    !--------------------------------------------------------------------------

    !integer :: memory_stat
    integer :: LASSLD

    LASSLD = 4*lmax
    harmonics%LASSLD = LASSLD
    allocate(harmonics%WG(LASSLD))
    allocate(harmonics%YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(harmonics%WG, harmonics%YRG, LMAX)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyMadelungHarmonics(harmonics)
    type(MadelungHarmonics), intent(inout) :: harmonics

    deallocate(harmonics%WG)
    deallocate(harmonics%YRG)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  subroutine createMadelungLatticeData(mlattice, NMAXD, ISHLD)
    type(MadelungLatticeData), intent(inout) :: mlattice
    integer, intent(in) :: NMAXD
    integer, intent(in) :: ISHLD
    !--------------------------------------------------------------------------

    mlattice%NMAXD = NMAXD
    mlattice%ISHLD = ISHLD

    allocate(mlattice%GN(3,NMAXD))
    allocate(mlattice%RM(3,NMAXD))
    allocate(mlattice%NSG(ISHLD))
    allocate(mlattice%NSR(ISHLD))

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyMadelungLatticeData(mlattice)
    type(MadelungLatticeData), intent(inout) :: mlattice
    !--------------------------------------------------------------------------

    deallocate(mlattice%GN)
    deallocate(mlattice%RM)
    deallocate(mlattice%NSG)
    deallocate(mlattice%NSR)

  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  subroutine createMadelungClebschData(clebsch, LMXSPD, LMPOTD)
    type(MadelungClebschData), intent(inout) :: clebsch
    integer, intent(in) :: LMXSPD
    integer, intent(in) :: LMPOTD
    !--------------------------------------------------------------------------

    allocate(clebsch%CLEB(LMXSPD*LMPOTD))
    allocate(clebsch%LOFLM(LMXSPD))
    allocate(clebsch%ICLEB(LMXSPD*LMPOTD,3))
  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyMadelungClebschData(clebsch)
    type(MadelungClebschData), intent(inout) :: clebsch
    !--------------------------------------------------------------------------

    deallocate(clebsch%CLEB)
    deallocate(clebsch%LOFLM)
    deallocate(clebsch%ICLEB)

  endsubroutine ! destroy

endmodule MadelungCalculator_mod
