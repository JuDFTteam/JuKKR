!> Instructions:
!> Create Madelung calculator
!> Create Madelung lattice sum(s) and configure with Madelung calculator
!> use...
!> destroy Madelung lattice sum(s)
!> destroy Madelung calculator

module MadelungCalculator_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private

  public :: MadelungCalculator, MadelungClebschData, MadelungLatticeSum, create, destroy
  public :: calculateMadelungLatticeSum, calc_dfac
  public :: createMadelungLatticeSum, destroyMadelungLatticeSum ! deprecated
  public :: createMadelungCalculator, destroyMadelungCalculator ! deprecated

  !----------------------------------------------------------------------------
  type MadelungLatticeData
    integer :: ngmax, nrmax
    integer :: nshlg, nshlr
    integer, allocatable :: nsg(:), nsr(:)
    double precision, allocatable :: gn(:,:), rm(:,:)
    integer :: nmaxd, ishld
  endtype

  !----------------------------------------------------------------------------
  type MadelungHarmonics
    double precision, allocatable :: wg(:)
    double precision, allocatable :: yrg(:,:,:)
    integer :: lassld
  endtype

  !----------------------------------------------------------------------------
  type MadelungClebschData
    integer :: iend
    double precision, allocatable :: cleb(:)
    integer, allocatable :: icleb(:,:)
    integer, allocatable :: loflm(:)
    integer :: nclebd
  endtype

  !----------------------------------------------------------------------------
  type MadelungCalculator
    double precision :: alat
    integer :: lpot
    double precision :: volume0
    double precision, allocatable :: dfac(:,:)
    type(MadelungLatticeData) :: mlattice
    type(MadelungClebschData) :: clebsch
    integer :: lmxspd
    integer :: lmpotd
    integer :: lassld
  endtype

  !----------------------------------------------------------------------------
  type MadelungLatticeSum
    integer :: num_atoms
    double precision, allocatable :: smat(:,:)
    type(MadelungCalculator), pointer :: madelung_calc
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
  subroutine createMadelungCalculator(madelung_calc, lmax, alat, rmax, gmax, bravais, nmaxd, ishld)
    use lattice_mod, only: lattix99
    type(MadelungCalculator), intent(inout) :: madelung_calc
    integer, intent(in) :: lmax
    double precision, intent(in) :: alat, rmax, gmax, bravais(3,3)
    integer, intent(in) :: nmaxd, ishld

    type(MadelungHarmonics) :: harmonics
    double precision :: recbv(3,3)
    integer :: lpot, lmxspd, lmpotd, lassld
    integer :: memory_stat
    integer, parameter :: no_output = 1

    lpot = 2*lmax
    lassld = 4*lmax
    lmxspd = (2*lpot+1)**2
    lmpotd = (lpot+1)**2

    madelung_calc%alat = alat
    madelung_calc%lpot = lpot
    madelung_calc%lmxspd = lmxspd
    madelung_calc%lmpotd = lmpotd
    madelung_calc%lassld = lassld

    allocate(madelung_calc%dfac(0:lpot,0:lpot), stat=memory_stat)
    if (memory_stat /= 0) die_here("Allocation error, DFAC")

    call lattix99(alat, bravais, recbv, madelung_calc%volume0, .false.)

    call createMadelungHarmonics(harmonics, lmax)

    call createMadelungLatticeData(madelung_calc%mlattice, NMAXD, ISHLD)

    call createMadelungClebschData(madelung_calc%clebsch, LMXSPD, LMPOTD)

    ! now wrap nasty function call
    call madelung3d(lpot, harmonics%yrg, harmonics%wg, alat, rmax, gmax, bravais, recbv, &
     lmxspd, lassld, lpot, lmpotd, nmaxd, ishld, lmpotd, &
     madelung_calc%clebsch%cleb,   madelung_calc%clebsch%icleb, madelung_calc%clebsch%iend, &
     madelung_calc%clebsch%nclebd, madelung_calc%clebsch%loflm, madelung_calc%dfac, &
     madelung_calc%mlattice%ngmax, madelung_calc%mlattice%nrmax, &
     madelung_calc%mlattice%nsg,   madelung_calc%mlattice%nsr, &
     madelung_calc%mlattice%nshlg, madelung_calc%mlattice%nshlr, &
     madelung_calc%mlattice%gn,    madelung_calc%mlattice%rm, no_output)

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

    type(MadelungCalculator), pointer :: mad_ptr

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

    type(MadelungCalculator), pointer :: madelung_calc
    
    madelung_calc => madelung_sum%madelung_calc

    call STRMAT(madelung_calc%ALAT,madelung_calc%LPOT,madelung_sum%num_atoms,madelung_calc%mlattice%NGMAX, &
    madelung_calc%mlattice%NRMAX,madelung_calc%mlattice%NSG,madelung_calc%mlattice%NSR, &
    madelung_calc%mlattice%NSHLG,madelung_calc%mlattice%NSHLR, &
    madelung_calc%mlattice%GN,madelung_calc%mlattice%RM, RBASIS, madelung_sum%SMAT, &
    madelung_calc%VOLUME0, &
    madelung_calc%LASSLD,madelung_calc%LMXSPD,madelung_sum%num_atoms,atom_index)

  endsubroutine ! calc

  !----------------------------------------------------------------------------
  !
  ! --> calculate:                               (2*(l+l')-1)!!
  !                 dfac(l,l') = (4pi)**2 *  ----------------------
  !                                          (2*l+1)!! * (2*l'+1)!!
  subroutine calc_dfac(dfac, lpot)
    use Constants_mod, only: pi
    double precision, intent(inout) :: dfac(0:lpot, 0:lpot)
    integer, intent(in) :: lpot

    double precision :: fpi
    integer :: l1, l2
    
    fpi = 4.d0*pi

    dfac(0,0) = fpi*fpi
    do l1 = 1, lpot
      dfac(l1,0) = dfac(l1-1,0)*dble(2*l1-1)/dble(2*l1+1)
      dfac(0,l1) = dfac(l1,0) ! symmetric
      do l2 = 1, l1
        dfac(l1,l2) = dfac(l1,l2-1)*dble(2*(l1+l2)-1)/dble(2*l2+1)
        dfac(l2,l1) = dfac(l1,l2) ! symmetric
      enddo ! l2
    enddo ! l1
    
  endsubroutine ! calc

  !----------------------------------------------------------------------------
  subroutine initMadelungClebschData(clebsch, lmax)
    type(MadelungClebschData), intent(inout) :: clebsch
    integer, intent(in) :: lmax

    integer lpot, lassld, lmxspd, lmpotd, nclebd, l, m, i
    type(MadelungHarmonics) :: harmonics

    lpot = 2*lmax
    lassld = 4*lmax
    lmxspd = (2*lpot+1)**2
    lmpotd = (lpot+1)**2
    nclebd = lmxspd*lmpotd

    i = 1
    do l = 0, 2*lpot
      do m = -l, l
        clebsch%loflm(i) = l
        i = i + 1
      enddo ! m
    enddo ! l

    call createMadelungHarmonics(harmonics, lmax)

    call madelgaunt(lpot, harmonics%yrg, harmonics%wg, clebsch%cleb, clebsch%icleb, clebsch%iend, lassld, nclebd)
    clebsch%nclebd = nclebd

    call destroyMadelungHarmonics(harmonics)
  endsubroutine ! init

!============= Helper routines =============================================

  !----------------------------------------------------------------------------
  subroutine createMadelungHarmonics(harmonics, lmax)
    use Harmonics_mod, only: Gaunt2 ! initialization of wg and yrg
    type(MadelungHarmonics), intent(inout) :: harmonics
    integer, intent(in) :: lmax

    integer :: lassld!, memory_stat

    lassld = 4*lmax
    harmonics%lassld = lassld
    allocate(harmonics%wg(lassld))
    allocate(harmonics%yrg(lassld,0:lassld,0:lassld))

    call gaunt2(harmonics%wg, harmonics%yrg, lmax)

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

  
  
  !>    @param print_info  0=print some info  1=print nothing
  subroutine madelung3d(lpot, yrg, wg, alat, rmax, gmax, bravais, recbv, &
      lmxspd,lassld,lpotd,lmpotd, nmaxd, ishld, &
      lmpot,cleb,icleb,iend, nclebd,loflm,dfac, &
      ngmax,nrmax,nsg,nsr,nshlg,nshlr,gn,rm, print_info)
    use Constants_mod, only: pi
    ! **********************************************************************
    ! *                                                                    *
    ! * This subroutine calculates the Madelung potential coefficients     *
    ! *                                                                    *
    ! **********************************************************************
    integer, intent(in) :: lpot, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, print_info
    double precision, intent(in) :: alat, rmax, gmax, bravais(3,3), recbv(3,3)
    double precision, intent(in) :: yrg(lassld,0:lassld,0:lassld), wg(lassld)
    integer, allocatable, intent(out) :: nsg(:), nsr(:) ! nsg(ishld), nsr(ishld)
    double precision, allocatable, intent(out) :: gn(:,:), rm(:,:) ! gn(3,nmaxd), rm(3,nmaxd)
    integer, intent(out) :: ngmax, nrmax, nshlg, nshlr
    double precision, intent(out) :: cleb(lmxspd*lmpotd) ! Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
    double precision, intent(out) :: dfac(0:lpotd,0:lpotd)
    integer, intent(out) :: icleb(lmxspd*lmpotd,3), iend
    integer, intent(out) :: loflm(:)
    integer, intent(inout) :: nclebd, lmpot ! has intent in?
    
    double precision :: fpi
    integer :: i, iprint, l, m, l1, l2

    iprint = 0
    nclebd = lmxspd*lmpotd
    fpi = 4.d0*pi
    lmpot = (lpot+1)**2

    ! --> determine the l-value for given lm
    assert( size(loflm, 1) >= (2*lpot+1)**2 )
    i = 0
    do l = 0, 2*lpot
      do m = -l, l
        i = i + 1
        assert(i == l*l + l + m + 1)
        loflm(i) = l
      enddo ! m
    enddo ! l
    
    ! --> calculate:                            (2*(l+l')-1)!!
    !                 dfac(l,l') = 4pi**2 *  ----------------------
    !                                       (2*l+1)!! * (2*l'+1)!!

    dfac(0,0) = fpi*fpi
    do l1 = 1, lpot
      dfac(l1,0) = dfac(l1-1,0)*dble(2*l1-1)/dble(2*l1+1)
      dfac(0,l1) = dfac(l1,0) ! symmetric
      do l2 = 1, l1
        dfac(l1,l2) = dfac(l1,l2-1)*dble(2*(l1+l2)-1)/dble(2*l2+1)
        dfac(l2,l1) = dfac(l1,l2) ! symmetric
      enddo ! l2
    enddo ! l1

    ! --> calculate the gaunt coefficients
    call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)

    call lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, nsr, rmax, gmax, gn, rm, iprint, nmaxd, ishld, print_info)
    
  endsubroutine madelung3d

  
  ! **********************************************************************
  ! *                                                                    *
  ! *  generate lattice vectors of direct and reciprocal space from      *
  ! *  basic translation vectors br                                      *
  ! *                                                                    *
  ! *  alat            : lattice constant                                *
  ! *  br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors                *
  ! *                    *** in a.u. ****                                *
  ! *  rmax            : maximum radius in real space        (input)     *
  ! *  gmax            : maximum radius in reciprocal space  (input)     *
  ! *  ngmax           : Number of reciprocal lattice vectors            *
  ! *  gn(3,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
  ! *  nrmax           : Number of real lattice vectors                  *
  ! *  rm(3,nmaxd)     : x,y,z  of real space vectors                    *
  ! *  nshlg           : shells in reciprocal space                      *
  ! *  nshlr           : shells in real space                            *
  ! *  nsg,nsr         : integer arrays, number of atoms in each shell   *
  ! *                                                                    *
  ! *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
  ! *  one it is used only locally (GNR/RMR)       v.popescu May 2004    *
  ! *                                                                    *
  ! **********************************************************************

subroutine lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, nsr, rmax, gmax, gn, rm, iprint, nmaxd, ishld, print_info)
  use Constants_mod, only: pi
  implicit none
  integer, intent(in) :: iprint, print_info
  integer, intent(in) :: nmaxd, ishld ! old dimensioning
  double precision, intent(in) :: alat !< lattice constant
  
  integer, intent(out) :: ngmax
  double precision, intent(in) :: gmax
  double precision, intent(in) :: recbv(3,3) !< reciprocal basis vectors
  double precision, allocatable, intent(out) :: gn(:,:) ! (3,nmaxd)
  integer, allocatable, intent(out) :: nsg(:) ! (ishld)
  integer, intent(out) :: nshlg
  
  integer, intent(out) :: nrmax
  double precision, intent(in) :: rmax
  double precision, intent(in) :: bravais(3,3) !< bravais matrix == real space basis vectors
  double precision, allocatable, intent(out) :: rm(:,:) ! (3,nmaxd)
  integer, allocatable, intent(out) :: nsr(:) ! (ishld)
  integer, intent(out) :: nshlr

  integer :: i, n, l, k, ist, numr(3), numg(3)
  double precision :: bg(3,3), absg2(3), absgm
  double precision :: br(3,3), absr2(3), absrm
  double precision, allocatable :: gnr(:), rmr(:) ! gnr(nmaxd), rmr(nmaxd)

  if (print_info == 0) write(6,'(5X,2A,/)') '< LATTICE3D > : ','generating direct/reciprocal lattice vectors'

  ! --> basic trans. vectors and basis vectors
  br(1:3,1:3) = bravais(1:3,1:3)*alat

  ! --> generate primitive vectors BG of reciprocal space
  bg(1:3,1:3) = recbv(1:3,1:3)*(2.d0*pi/alat)

  ! --> estimate no. of lattice vectors, todo: compare Fleur routine boxdim
  do i = 1, 3
    absr2(i) = sum(br(1:3,i)**2)
    absg2(i) = sum(bg(1:3,i)**2)
  enddo ! i
  
  absrm = (2.d0*pi)/sqrt(max(absr2(1), absr2(2), absr2(3)))
  absgm = (2.d0*pi)/sqrt(max(absg2(1), absg2(2), absg2(3)))
  numr(1:3) = 2*idint(rmax/absgm) + 3 ! Warning: cross term here: direct depends on reciprocal
  numg(1:3) = 2*idint(gmax/absrm) + 3 ! Warning: cross term here: reciprocal depends on direct

! ! write(0,'(a)') trim(__FILE__+"alat"+alat)
! ! write(0,'(a)') trim(__FILE__+"bravais"+reshape(bravais, [9]))
! ! write(0,'(a)') trim(__FILE__+"br"+reshape(br, [9]))
! ! write(0,'(a)') trim(__FILE__+"recbv"+reshape(recbv, [9]))
! ! write(0,'(a)') trim(__FILE__+"bg"+reshape(bg, [9]))
! ! write(0,*) trim(__FILE__+"loop"+numg+" in g-space from gmax ="+gmax)
! ! write(0,*) trim(__FILE__+"loop"+numr+" in r-space from rmax ="+rmax)

  ! **********************************************************************
  !                 generate lattice vectors of real space
  ! **********************************************************************
  ist = count_vectors_in_sphere(numr, br, rmax, & ! inputs
    vec=rm, rad=rmr, nvecs=nrmax, nshells=nshlr, nsh=nsr, & ! outputs
    tol_origin=1.d-6, tol_newshell=1.d-6, space='r') ! config
  
  if (nshlr <= 1) die_here("cut-off radius RMAX too small!")
  
  ! **********************************************************************
  !                 generate lattice vectors of reciprocal space
  ! **********************************************************************
  ist = count_vectors_in_sphere(numg, bg, gmax, & ! inputs
    vec=gn, rad=gnr, nvecs=ngmax, nshells=nshlg, nsh=nsg, & ! outputs
    tol_origin=1.d-6, tol_newshell=1.d-7, space='g') ! config
  
  if (nshlg <= 1) die_here("cut-off radius GMAX too small!")

  if (print_info == 0) then
     write(6, fmt=99002)
     write(6, fmt=99003) 'Direct  lattice', nrmax, nshlr, rmr(nrmax)
     write(6, fmt=99003) 'Recipr. lattice', ngmax, nshlg, gnr(ngmax)
     write(6, fmt=99004)
  endif

  if (iprint < 3) return

  k = 0
  write(6, fmt=99005) 'real-space'
  do l = 1, nshlr
    write(6, fmt=99006) l, nsr(l), rmr(k+1), rm(1:3,k+1)
    do n = 2, nsr(l)
      write(6, fmt=99007) rm(1:3,k+n)
    enddo ! n
    if (l /= nshlr) write(6, fmt=99008)
    k = k + nsr(l)
  enddo ! l
  write(6, fmt=99009)
  k = 0
  write(6, fmt=99005) 'reciprocal'
  do l = 1, nshlg
    write(6, fmt=99006) l, nsg(l), gnr(k+1), gn(1:3,k+1)
    do n = 2, nsg(l)
      write(6, fmt=99007) gn(1:3,k+n)
    enddo ! n
    if (l /= nshlg) write(6, fmt=99008)
    k = k + nsg(l)
  enddo ! l
  write(6, fmt=99009)

99002 format (10x,'               vectors  shells  max. R ',/,10x,'               ------------------------------')
99003 format (10x,a,i7,2x,i6,2x,f9.5)
99004 format (10x,'               ------------------------------',/)
99005 format (10x,55('+'),/,18x,'generated ',a,' lattice vectors',/,10x,55('+'),/,10x,'shell Nvec    radius          x         y         z',/,10x,55('-'))
99006 format (10x,i5,i5,f12.6,2x,3f10.5)
99007 format (34x,3f10.5)
99008 format (13x,52('-'))
99009 format (10x,55('+'),/)
endsubroutine ! lattice3d
  
  
  integer function count_vectors_in_sphere(num, bm, dmax, vec, rad, nvecs, nshells, nsh, tol_origin, tol_newshell, space) result(ist)
    integer, intent(in) :: num(3)
    double precision, intent(in) :: bm(3,3) !< bravais matrix or reciprocal space basis
    double precision, intent(in) :: dmax !< cutoff radius of the sphere
    double precision, allocatable, intent(out), optional :: vec(:,:), rad(:) ! check where this data is needed --> maybe better in one array (0:3,:)
    integer, intent(out) :: nvecs, nshells ! number
    integer, allocatable, intent(out), optional :: nsh(:)
    double precision, intent(in) :: tol_origin, tol_newshell
    character, intent(in) :: space ! for debug: 'r' or 'g'
    
    integer :: numh(3), i, i1, i2, i3, ivmin, iminl(1), i01, ish, nshl
    double precision :: dmax2, v2, vmin, very_large, da, db, vx(3), vxy(3), vxyz(3)
    double precision, allocatable :: cv(:,:), d2(:)
    integer, allocatable :: nvis(:) ! tmp for nsh
  
    numh = num/2 + 1
    dmax2 = dmax**2 ! radius^2
    
    ! **********************************************************************
    !                 generate lattice vectors of real or reciprocal space
    ! **********************************************************************
   
    do i01 = 0, 1 ! loop runs twice: iteration #0: count & allocate, iteration #1: store
    
      i = 0 ! init
      do i1 = 1 - numh(1), num(1) - numh(1)
        vx(1:3) = i1*bm(1:3,1)
        do i2 = 1 - numh(2), num(2) - numh(2)
          vxy(1:3) = i2*bm(1:3,2) + vx(1:3) 
          do i3 = 1 - numh(3), num(3) - numh(3)
            vxyz(1:3) = i3*bm(1:3,3) + vxy(1:3) 
            
            v2 = sum(vxyz(1:3)**2)

            if (v2 <= dmax2) then
              i = i + 1 ! count up
              if (i01 == 1) then
                cv(1:3,i) = vxyz(1:3) ! store vector components
                d2(i) = v2 ! and length^2 in the second iteration
              endif
            endif
            
          enddo  ! i3
        enddo ! i2
      enddo ! i1
      
      if (i01 == 0) then ! after first iterations (count iteration)
        nvecs = i
        allocate(cv(1:3,nvecs), d2(nvecs), nvis(nvecs), stat=ist)
      endif

    enddo ! i01    
    ! ======================================================================
!!! write(0,*) trim(__FILE__+"loop"+num+" in"+space-"-space finds"+nvecs+"vectors")
    !
    ! --> sort vectors in order of increasing absolute value
    !
    
    if (present(vec)) then; deallocate(vec, stat=ist); allocate(vec(1:3,nvecs), stat=ist); endif
    if (present(rad)) then; deallocate(rad, stat=ist); allocate(rad(nvecs), stat=ist); endif
    
    ! warning: this method of computing the shell structure
    ! scales N^2 with N=nvecs, better would be sorting with N*log(N)
    
    da = tol_origin
    
    ish = 0
    nshl = -1 ! start first shell with one less
    very_large = dmax**2 + 9.d9
    
    do i = 1, nvecs
    
      iminl = minloc(d2) ! find the location of the smallest element in d2
      ivmin = iminl(1) ! pass index
      vmin = sqrt(d2(ivmin))

      nshl = nshl + 1 ! increase the number of points in this shell
      
      if (present(vec)) vec(1:3,i) = cv(1:3,ivmin) ! store vector
      if (present(rad)) rad(i) = vmin ! store radius
      
      db = vmin
      ! ----------------------------------------------------------------------
      if (db > da + tol_newshell) then
        ! create a new shell
        ish = ish + 1 ! create a shell index of the current shell
        nvis(ish) = nshl ! store the number of points in the last shell
!!! write(0,*) trim(__FILE__+"create shell"+ish+"with"+nshl+"points in"+space-"-space")
        
        nshl = 0 ! init number of points for the new shell
        da = db
      endif
      ! ----------------------------------------------------------------------
      
      d2(ivmin) = very_large ! take it out of the loop
    enddo ! i
    deallocate(cv, d2, stat=ist) ! free local arrays

    ! close current shell
    ish = ish + 1
    nshl = nshl + 1
    nvis(ish) = nshl
!!! write(0,*) trim(__FILE__+"create shell"+ish+"with"+nshl+"points in"+space-"-space (last)")

    nshells = ish ! export the max number of shells
    
    if (present(nsh)) then
      deallocate(nsh, stat=ist)
      allocate(nsh(nshells), stat=ist)
      nsh = nvis(1:nshells)
    endif ! present nsh(:)

    deallocate(nvis, stat=ist)
    
  endfunction ! count_vectors_in_sphere
  
  
  subroutine madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)
    integer, intent(in) :: lpot, lassld, nclebd
    integer, intent(out) :: iend
    double precision, intent(in) :: yrg(lassld,0:lassld,0:lassld), wg(lassld)
    double precision, intent(out) :: cleb(:) ! (nclebd) ! Attention: Dimension NCLEBD appears sometimes as NCLEB1, an empirical factor that has to be optimized
    integer, intent(out) :: icleb(:,:) ! (nclebd,3)

    double precision :: clecg, factor, s
    integer :: i, j, l1, l2, l3, m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s, mclebd
    !
    ! --> set up of the gaunt coefficients with an index field
    !     recognize that they are needed here only for l3=l1+l2
    !
    if (2*lpot > lassld) then
      write (6,*) 'Dim ERROR in MADELGAUNT -- 2*LPOT > LASSLD', 2*lpot, lassld
      stop
    endif
    
    mclebd = min(size(cleb, 1), size(icleb, 1))

    i = 0
    do l1 = 0, lpot
      do l2 = 0, lpot
        l3 = l1 + l2
        do m1 = -l1, l1
          do m2 = -l2, l2
            do m3 = -l3, l3
              m1s = sign(1,m1)
              m2s = sign(1,m2)
              m3s = sign(1,m3)
              if (m1s*m2s*m3s >= 0) then
                m1a = abs(m1)
                m2a = abs(m2)
                m3a = abs(m3)
                factor = 0.d0
                if (m1a+m2a == m3a) factor = factor + 0.125d0*(3*m3s+sign(1,-m3))
                if (m1a-m2a == m3a) factor = factor + 0.25d0*m1s
                if (m2a-m1a == m3a) factor = factor + 0.25d0*m2s
                if (factor /= 0.d0) then
                  if (m1s*m2s /= 1 .or. m2s*m3s /= 1 .or. m1s*m3s /= 1) factor = -factor
                  !
                  s = 0.d0
                  do j = 1,lassld
                    s = s + wg(j)*yrg(j,l1,m1a)*yrg(j,l2,m2a)*yrg(j,l3,m3a)
                  enddo ! j
                  !
                  clecg = s*factor
                  if (abs(clecg) > 1.d-10) then
                    i = i + 1
                    if (i <= mclebd) then
                      cleb(i) = clecg
                      icleb(i,1) = l1*l1 + l1 + m1 + 1
                      icleb(i,2) = l2*l2 + l2 + m2 + 1
                      icleb(i,3) = l3*l3 + l3 + m3 + 1
                    endif ! respect array bounds
                  endif ! abs(clecg) > 1.d-10
                endif ! factor /= 0.d0
              endif ! m1s*m2s*m3s >= 0
            enddo ! m3
          enddo ! m2
        enddo ! m1
      enddo ! l2
    enddo ! l1
    iend = i
    
    if (iend > mclebd) then
      write (6,fmt='(3I10)') i, nclebd, mclebd
      stop ' Dim stop in MADELGAUNT '
    endif
    
  endsubroutine madelgaunt
  
endmodule MadelungCalculator_mod
