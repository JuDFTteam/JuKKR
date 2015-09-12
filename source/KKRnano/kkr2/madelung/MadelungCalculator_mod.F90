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
  public :: calculateMadelungLatticeSum, createDfac ! deprecated
  public :: createMadelungLatticeSum, destroyMadelungLatticeSum ! deprecated
  public :: createMadelungCalculator, destroyMadelungCalculator ! deprecated

  !----------------------------------------------------------------------------
  type MadelungLatticeData
    integer :: ngmax, nrmax
    integer :: nshlg, nshlr
    integer, allocatable :: nsg(:), nsr(:)
    double precision, allocatable :: gn(:,:), rm(:,:)
  endtype

  !----------------------------------------------------------------------------
  type MadelungHarmonics
    double precision, allocatable :: wg(:)
    double precision, allocatable :: yrg(:,:,:)
    integer :: lassld = -1
  endtype

  !----------------------------------------------------------------------------
  type MadelungClebschData
    integer :: iend
    double precision, allocatable :: cleb(:)
    integer, allocatable :: icleb(:,:) !< dimension(nclebd,3)
    integer, allocatable :: loflm(:)
    integer :: nclebd
  endtype

  !----------------------------------------------------------------------------
  type MadelungCalculator
    integer :: lpot
    double precision :: alat
    double precision :: volume0
    double precision, allocatable :: dfac(:,:)
    type(MadelungLatticeData) :: lattice
    type(MadelungClebschData) :: clebsch
    integer :: lmxspd ! == (2*lpot+1)**2
    integer :: lmpotd ! == (lpot+1)**2
    integer :: lassld ! == 2*lpot
  endtype

  !----------------------------------------------------------------------------
  type MadelungLatticeSum
    integer :: num_atoms
    double precision, allocatable :: smat(:,:)
    type(MadelungCalculator), pointer :: madelung_calc
  endtype

  
  interface create
    module procedure createMadelungCalculator, createMadelungClebschData, createMadelungLatticeSum, createDfac
  endinterface
  
  interface destroy
    module procedure destroyMadelungCalculator, destroyMadelungClebschData, destroyMadelungLatticeSum, &
                     destroyMadelungHarmonics, destroyMadelungLatticeData
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  !> Creates a MadelungCalculator object.
  !> Don't forget to free resources with destroyMadelungCalculator
  subroutine createMadelungCalculator(self, lmax, alat, rmax, gmax, bravais)
    use lattice_mod, only: lattix99
    type(MadelungCalculator), intent(inout) :: self
    integer, intent(in) :: lmax
    double precision, intent(in) :: alat, rmax, gmax, bravais(3,3)

    double precision :: recbv(3,3)
    integer :: lpot, ist
    
    self%alat = alat
    
    lpot = 2*lmax
    
    self%lpot = lpot
    self%lmxspd = (2*lpot+1)**2
    self%lmpotd = (lpot+1)**2
    self%lassld = 4*lmax ! == 2*lpot
    
    call createDfac(self%dfac, lpot)

    call lattix99(alat, bravais, recbv, self%volume0, .false.)

    call createMadelungClebschData(self%clebsch, lmax)

#define ml self%lattice
    call lattice3d(alat, bravais, recbv, ml%ngmax, ml%nrmax, ml%nshlg, ml%nshlr, ml%nsg, ml%nsr, rmax, gmax, ml%gn, ml%rm, iprint=0, print_info=1) ! 1:no_output
#undef  ml

  endsubroutine ! create

  
  
!   !>    @param print_info  0=print some info  1=print nothing
!   subroutine madelung3d(lpot, yrg, wg, alat, rmax, gmax, bravais, recbv, &
!       lmxspd,lassld,lpotd,lmpotd, &
!       lmpot,cleb,icleb,iend, nclebd,loflm, &
!       ngmax,nrmax,nsg,nsr,nshlg,nshlr,gn,rm, print_info)
!     use Constants_mod, only: pi
!     ! **********************************************************************
!     ! *                                                                    *
!     ! * This subroutine calculates the Madelung potential coefficients     *
!     ! *                                                                    *
!     ! **********************************************************************
!     integer, intent(in) :: lpot, lmxspd, lassld, lpotd, lmpotd, print_info
!     double precision, intent(in) :: alat, rmax, gmax, bravais(3,3), recbv(3,3)
!     double precision, intent(in) :: yrg(lassld,0:lassld,0:lassld), wg(lassld)
!     integer, allocatable, intent(out) :: nsg(:), nsr(:) ! nsg(ishld), nsr(ishld)
!     double precision, allocatable, intent(out) :: gn(:,:), rm(:,:) ! gn(3,nmaxd), rm(3,nmaxd)
!     integer, intent(out) :: ngmax, nrmax, nshlg, nshlr
!     double precision, intent(out) :: cleb(lmxspd*lmpotd) ! Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
!     integer, intent(out) :: icleb(lmxspd*lmpotd,3), iend
!     integer, intent(out) :: loflm(:)
!     integer, intent(inout) :: nclebd, lmpot ! has intent in?
!     
!     integer :: i, iprint, l, m, l1, l2, icleb_dummy(1,3)
!     double precision :: cleb_dummy(1)
! 
!     iprint = 0
!     nclebd = lmxspd*lmpotd
!     lmpot = (lpot+1)**2
! 
! !     ! --> determine the l-value for given lm
! !     assert( size(loflm, 1) >= (2*lpot+1)**2 )
! !     i = 0
! !     do l = 0, 2*lpot
! !       do m = -l, l
! !         i = i + 1
! !         assert(i == l*l + l + m + 1)
! !         loflm(i) = l
! !       enddo ! m
! !     enddo ! l
! !     
! ! !   call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)
! !     iend = madelgaunt(lpot, yrg, wg, cleb, icleb, lassld) ! calculate the gaunt coefficients
! !     allocate(
! !     iend = madelgaunt(lpot, yrg, wg, cleb, icleb, lassld) ! calculate the gaunt coefficients
!     
!     call lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, nsr, rmax, gmax, gn, rm, iprint, print_info)
!     
!   endsubroutine madelung3d
  
  
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

    allocate(madelung_sum%smat(madelung_calc%lmxspd,num_atoms))
  endsubroutine ! create

  
  !----------------------------------------------------------------------------
  !> Destroys a MadelungCalculator object.
  subroutine destroyMadelungCalculator(madelung_calc)
    type(MadelungCalculator), intent(inout) :: madelung_calc

    call destroyMadelungLatticeData(madelung_calc%lattice)
    call destroyMadelungClebschData(madelung_calc%clebsch)
    deallocate(madelung_calc%dfac)
  endsubroutine ! destroy
  
  !----------------------------------------------------------------------------
  !> Destroys and frees storage of Madelung lattice sum.
  subroutine destroyMadelungLatticeSum(madelung_sum)
    type(MadelungLatticeSum), intent(inout) :: madelung_sum

    deallocate(madelung_sum%smat)
  endsubroutine ! destroy  

  !----------------------------------------------------------------------------
  !> Calculates Lattice sum for atom 'atom_index' for positions 'rbasis'
  !>  needed for Madelung potential calculation.
  !>
  !> Needs a properly constructed MadelungLatticeSum object.
  !> STRMAT wrapper
  subroutine calculateMadelungLatticeSum(self, atom_index, rbasis)
    type(MadelungLatticeSum), intent(inout) :: self
    integer, intent(in) :: atom_index
    double precision, intent(in) :: rbasis(3,self%num_atoms)

    type(MadelungCalculator), pointer :: madelung_calc
    
#define mc self%madelung_calc
    call STRMAT(mc%alat, mc%lpot, self%num_atoms, mc%lattice%ngmax, &
      mc%lattice%nrmax, mc%lattice%nsg, mc%lattice%nsr, &
      mc%lattice%nshlg, mc%lattice%nshlr, &
      mc%lattice%gn, mc%lattice%rm, rbasis, self%smat, &
      mc%volume0, &
      mc%lassld, mc%lmxspd, self%num_atoms, atom_index)
#undef  mc

  endsubroutine ! calc

  !----------------------------------------------------------------------------
  !
  ! --> calculate:                               (2*(l+l')-1)!!
  !                 dfac(l,l') = (4pi)**2 *  ----------------------
  !                                          (2*l+1)!! * (2*l'+1)!!
  subroutine createDfac(dfac, lpot)
    use Constants_mod, only: pi
    double precision, allocatable, intent(out) :: dfac(:,:) !(0:lpot,0:lpot)
    integer, intent(in) :: lpot

    integer :: l1, l2, ist

    deallocate(dfac, stat=ist)
    allocate(dfac(0:lpot,0:lpot), stat=ist); if (ist /= 0) die_here("Allocation of dfac failed")
    
    dfac(0,0) = (4.d0*pi)**2
    do l1 = 1, lpot
      dfac(l1,0) = dfac(l1-1,0)*dble(2*l1-1)/dble(2*l1+1)
      dfac(0,l1) = dfac(l1,0) ! symmetric
      do l2 = 1, l1
        dfac(l1,l2) = dfac(l1,l2-1)*dble(2*(l1+l2)-1)/dble(2*l2+1)
        dfac(l2,l1) = dfac(l1,l2) ! symmetric
      enddo ! l2
    enddo ! l1
    
  endsubroutine ! create


  !----------------------------------------------------------------------------
  subroutine createMadelungHarmonics(harmonics, lmax)
    use Harmonics_mod, only: Gaunt2 ! initialization of wg and yrg
    type(MadelungHarmonics), intent(inout) :: harmonics
    integer, intent(in) :: lmax

    integer :: lassld, ist

    lassld = 4*lmax
    if (harmonics%lassld == lassld) return ! already done

    deallocate(harmonics%wg, harmonics%yrg, stat=ist)
    allocate(harmonics%wg(lassld), harmonics%yrg(lassld,0:lassld,0:lassld), stat=ist)
    if (ist /= 0) die_here("Allocation of YRG with lmax="-lassld+"failed!")
    
    call gaunt2(harmonics%wg, harmonics%yrg, lmax)
    
    harmonics%lassld = lassld
  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyMadelungHarmonics(harmonics)
    type(MadelungHarmonics), intent(inout) :: harmonics

    deallocate(harmonics%wg, harmonics%yrg)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  subroutine destroyMadelungLatticeData(lattice)
    type(MadelungLatticeData), intent(inout) :: lattice

    deallocate(lattice%gn, lattice%nsg)
    deallocate(lattice%rm, lattice%nsr)
  endsubroutine ! destroy
 
  !----------------------------------------------------------------------------
  subroutine createMadelungClebschData(self, lmax)
    type(MadelungClebschData), intent(inout) :: self
    integer, intent(in) :: lmax

    integer :: icleb_dummy(1,3), iend
    double precision :: cleb_dummy(1)
    integer :: l, m, i, ist, lpot

    type(MadelungHarmonics) :: hmx

    deallocate(self%cleb, self%icleb, self%loflm, stat=ist)
    lpot = 2*lmax

    call createMadelungHarmonics(hmx, lmax)
    
    iend = madelgaunt(lpot, hmx%yrg, hmx%wg, hmx%lassld, cleb_dummy, icleb_dummy) ! only count the gaunt coefficients
    allocate(self%cleb(iend), self%icleb(iend,3))
    self%nclebd = iend
    self%iend = madelgaunt(lpot, hmx%yrg, hmx%wg, hmx%lassld, self%cleb, self%icleb) ! calculate and store the gaunt coefficients
    if (iend /= self%iend) die_here("first and second call did not return the same number of non-zero Gaunt coefficients!")
    
    call destroyMadelungHarmonics(hmx)
    
    allocate(self%loflm((2*lpot+1)**2))
    i = 0
    do l = 0, 2*lpot
      do m = -l, l
        i = i + 1
        assert(i == l*l + l + m + 1)
        self%loflm(i) = l
      enddo ! m
    enddo ! l
    
  endsubroutine ! create
  
  !----------------------------------------------------------------------------
  subroutine destroyMadelungClebschData(self)
    type(MadelungClebschData), intent(inout) :: self

    deallocate(self%cleb, self%loflm, self%icleb)
  endsubroutine ! destroy

  

  
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

subroutine lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, nsr, rmax, gmax, gn, rm, iprint, print_info)
  use Constants_mod, only: pi
  implicit none
  integer, intent(in) :: iprint, print_info
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

  
  character(len=*), parameter :: F97="(34x,3f10.5)", F98="(13x,52('-'))", F99="(10x,55('+'),/)", &
  F92="(25x,'vectors  shells  max. R ',/,25x,30('-'))", F93="(10x,a,i7,2x,i6,2x,f9.5)", F94="(25x,30('-'),/)", F96="(10x,i5,i5,f12.6,2x,3f10.5)", &
  F95="(10x,55('+'),/,18x,'generated ',a,' lattice vectors',/,10x,55('+'),/,10x,'shell Nvec    radius          x         y         z',/,10x,55('-'))"
  
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
     write(6, fmt=F92)
     write(6, fmt=F93) 'Direct  lattice', nrmax, nshlr, rmr(nrmax)
     write(6, fmt=F93) 'Recipr. lattice', ngmax, nshlg, gnr(ngmax)
     write(6, fmt=F94)
  endif

  if (iprint < 3) return

  k = 0
  write(6, fmt=F95) 'real-space'
  do l = 1, nshlr
    write(6, fmt=F96) l, nsr(l), rmr(k+1), rm(1:3,k+1)
    do n = 2, nsr(l)
      write(6, fmt=F97) rm(1:3,k+n)
    enddo ! n
    if (l /= nshlr) write(6, fmt=F98)
    k = k + nsr(l)
  enddo ! l
  write(6, fmt=F99)
  k = 0
  write(6, fmt=F95) 'reciprocal'
  do l = 1, nshlg
    write(6, fmt=F96) l, nsg(l), gnr(k+1), gn(1:3,k+1)
    do n = 2, nsg(l)
      write(6, fmt=F97) gn(1:3,k+n)
    enddo ! n
    if (l /= nshlg) write(6, fmt=F98)
    k = k + nsg(l)
  enddo ! l
  write(6, fmt=F99)

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
  
  
  integer function madelgaunt(lpot, yrg, wg, lassld, cleb, icleb) result(iend)
    integer, intent(in) :: lpot, lassld
    double precision, intent(in) :: yrg(lassld,0:lassld,0:lassld), wg(lassld)
    double precision, intent(out) :: cleb(:) ! (nclebd) ! Attention: Dimension NCLEBD appears sometimes as NCLEB1, an empirical factor that has to be optimized
    integer, intent(out) :: icleb(:,:) ! (nclebd,3)

    double precision :: clecg, factor, s
    integer :: i, j, l1, l2, l3, m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s, mclebd
    !
    ! --> set up of the gaunt coefficients with an index field
    !     recognize that they are needed here only for l3=l1+l2
    !
    if (2*lpot > lassld) die_here("madelgaunt 2*lpot ="+(2*lpot)+">"+lassld)
    
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
    
    if (iend > mclebd .and. mclebd > 1) & ! if the arrays are dimensioned 1, this is intentional for counting
      warn(6, "coefficients in MadelungGaunt have been truncated: arrays store only"+mclebd+"of"+iend) 
    
  endfunction madelgaunt
  
  
  ! **********************************************************************
  ! *                                                                    *
  !>*  calculation of lattice sums for l .le. 2*lpot :                   *
  !>*                                                                    *
  !>*                   ylm( q(i) - q(j) - rv )                          *
  !>*        sum      ===========================                        *
  !>*                 | q(i) - q(j) - rv |**(l+1)                        *
  !>*                                                                    *
  !>*         - summed over all lattice vectors rv  -                    *
  !>*                                                                    *
  !>*  ylm       : real spherical harmic to given l,m                    *
  !>*  q(i),q(j) : basis vectors of the unit cell                        *
  !>*                                                                    *
  !>*  in the case of i = j, rv = 0 is omitted.                          *
  !>*                                                                    *
  !>*  the ewald method is used to perform the lattice summations        *
  !>*  the splitting parameter lamda is set equal sqrt(pi)/alat          *
  !>*  (alat is the lattice constant) .                                  *
  !>*                                                                    *
  !>*  if the contribution of the last shell of the direct and the       *
  !>*  reciprocal lattice is greater than 1.0e-8 a message is written    *
  !>*                                                                    *
  !>*                                    b.drittler may 1989             *
  !>*                                                                    *
  !>*  Dimension of arrays gv,rv changed from (4,*) to (3,*), the 4th    *
  !>*  one not being used (see also lattice3d)     v.popescu May 2004    *
  ! *                                                                    *
  ! **********************************************************************

! OpenMP parallelised, needs threadsafe erfcex, gamfc and ymy E.R.

subroutine strmat(alat,lpot,naez,ngmax,nrmax,nsg,nsr,nshlg,nshlr, &
     gv,rv,qi0,smat,vol,lassld,lmxspd,naezd,i1) ! todo: remove naezd
  use Harmonics_mod, only: ymy
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  use Constants_mod, only: pi
  implicit none
  ! Parameters
  double complex, parameter :: CI=(0.d0,1.d0)
  double precision, parameter :: BOUND=1.d-8

  ! Arguments
  double precision, intent(in) :: alat
  double precision, intent(in) :: vol
  integer, intent(in) :: lpot
  integer, intent(in) :: naez
  integer, intent(in) :: ngmax
  integer, intent(in) :: nrmax
  integer, intent(in) :: nshlg
  integer, intent(in) :: nshlr
  integer, intent(in) :: lassld
  integer, intent(in) :: lmxspd
  integer, intent(in) :: naezd
  double precision, intent(in) :: gv(3,*)
  double precision, intent(in) :: qi0(3,*)
  double precision, intent(in) :: rv(3,*)
  double precision, intent(inout) :: smat(lmxspd,*)
  integer, intent(in) :: nsg(*)
  integer, intent(in) :: nsr(*)
  integer, intent(in) :: i1

 !local variables of strmat
 double complex :: bfac
 double precision :: alpha, beta
 double precision :: dq1
 double precision :: dq2
 double precision :: dq3
 double precision :: dqdotg
 double precision :: expbsq
 double precision :: fpi
 double precision :: g1
 double precision :: g2
 double precision :: g3
 double precision :: ga
 double precision :: lamda
 double precision :: r
 double precision :: r1
 double precision :: r2
 double precision :: r3
 double precision :: rfac
 double precision :: s
 integer :: i
 integer :: i2
 integer :: it
 integer :: l
 integer :: lm
 integer :: lmx
 integer :: lmxsp
 integer :: m
 integer :: nge
 integer :: ngs
 integer :: nre
 integer :: nrs
 integer :: nstart

 double complex :: stest(lmxspd)
 double precision :: g(0:lassld), ylm(lmxspd)

  lmx = 2*lpot
  lmxsp = (lmx+1)**2
  fpi = 4.0d0*pi

  ! --> choose proper splitting parameter

  lamda = sqrt(pi)/alat

  ! **********************************************************************
  !$omp parallel do private(I2,DQ1,DQ2,DQ3,STEST,LM,NSTART,IT, &
  !$omp                     NRS,NGS,NRE,NGE,I,R1,R2,R3, &
  !$omp                     YLM,R,ALPHA,G,RFAC,L,M, &
  !$omp                     G1,G2,G3,GA,BETA,EXPBSQ,DQDOTG,BFAC,S)
  do i2 = 1, naez
     !======================================================================
     dq1 = (qi0(1,i1) - qi0(1,i2)) * alat
     dq2 = (qi0(2,i1) - qi0(2,i2)) * alat
     dq3 = (qi0(3,i1) - qi0(3,i2)) * alat

     stest(1) = -sqrt(fpi)/vol/(4d0*lamda*lamda)
     do lm = 2,lmxsp
        stest(lm) = 0.0d0
     enddo

     ! --> exclude the origine and add correction if i1 == i2

     if ( i1 == i2 ) then
        stest(1) = stest(1) - lamda/pi
        nstart = 2
     else
        nstart = 1
     endif
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     ! --> loop first over n-1 shells of real and reciprocal lattice - then
     !     add the contribution of the last shells to see convergence

     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     do it = 1, 2
        if ( it == 1 ) then
           nrs = nstart
           ngs = 2
           nre = nrmax - nsr(nshlr)
           nge = ngmax - nsg(nshlg)
        else
           nrs = nre + 1
           ngs = nge + 1
           nre = nrmax
           nge = ngmax
        endif

        ! --> sum over real lattice

        ! ---------------------------------------------------------------------
        do i = nrs, nre
           r1 = dq1 - rv(1,i)
           r2 = dq2 - rv(2,i)
           r3 = dq3 - rv(3,i)

           call ymy(r1, r2, r3, r, ylm, lmx)
           alpha = lamda*r
           g(0:lmx) = gamfc(lmx, alpha, r)

           do l = 0, lmx
              rfac = g(l)/sqrt(pi)
              do m = -l, l
                 lm = l*(l+1) + m + 1
                 stest(lm) = stest(lm) + ylm(lm)*rfac
              enddo ! m
           enddo ! l
        enddo ! i
        ! ---------------------------------------------------------------------

        ! --> sum over reciprocal lattice

        ! ---------------------------------------------------------------------
        do i = ngs, nge
           g1 = gv(1,i)
           g2 = gv(2,i)
           g3 = gv(3,i)

           call ymy(g1,g2,g3,ga,ylm,lmx)
           beta = ga/lamda
           expbsq = exp(beta*beta*0.25d0)
           dqdotg = dq1*g1 + dq2*g2 + dq3*g3

           bfac = fpi*exp(CI*dqdotg)/(ga*ga*expbsq*vol)

           do l = 0,lmx
              do m = -l,l
                 lm = l*(l+1) + m + 1
                 stest(lm) = stest(lm) + ylm(lm)*bfac
              enddo ! m
              bfac = bfac*ga/dble(2*l+1)*(-CI)
           enddo ! l
        enddo ! i
        ! ---------------------------------------------------------------------
        if ( it == 1 ) then
           do lm = 1, lmxsp
              if (abs(dimag(stest(lm))) > BOUND) die_here("Imaginary contribution to REAL lattice sum")
              smat(lm,i2) = dble(stest(lm))
              stest(lm) = 0.d0
           enddo ! lm
        else

           ! --> test convergence

           do lm = 1, lmxsp
              s = dble(stest(lm))
              smat(lm,i2) = smat(lm,i2) + s
              !IF (2 < 1 .AND. ABS(S) > BOUND ) WRITE (6,FMT=99001) I1,I2, &
              !LM,ABS(S)
           enddo ! lm
        endif
        ! ---------------------------------------------------------------------
     enddo ! it
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  enddo ! I2 ! loop over all atoms
  !$omp endparallel do
  ! **********************************************************************

! 99001 format (5x,'WARNING : Convergence of SMAT(',i2,',',i2,') ', &
!        ' for LMXSP =',i3,' is ',1p,d8.2,' > 1D-8',/,15x, &
!        'You should use more lattice vectors (RMAX/GMAX)')
endsubroutine strmat
  
  
  
  function gamfc(lmax, alpha, r) result(glh)
    !----------------------------------------------------------------------
    !      calculation of convergence function
    !
    !       glh = i(alpha,l)/r**(l+1)*sqrt(pi)
    !
    !      with
    !            alpha = r times the splitting paramter lamda
    !      and
    !            i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *
    !
    !                                sum ( 2**i * x**(2i-1) / (2i-1)!! )
    !                              1..i..l
    !
    !
    ! Note: gamfc( alpha -> 0, ... ) => glh = sqrt(pi) for all l  E.R.
    !-----------------------------------------------------------------------
    double precision :: glh(0:lmax) ! result
    integer, intent(in) :: lmax
    double precision, intent(in) :: alpha, r

    double precision, external :: erfcex
    double precision :: arg, facl, fex
    integer :: l

    arg = alpha*alpha
    facl = 2.0d0*alpha

    glh(0) = erfcex(alpha)
    !---> recursion
    do l = 1, lmax
      glh(l) = glh(l-1) + facl
      facl = facl*arg/(l + 0.5d0)
    enddo ! l

    ! if arg is to big then cannot calculate 1/exp(arg) !!

    fex = exp(-arg)

    do l = 0, lmax
      fex = fex/r
      glh(l) = glh(l)*fex
    enddo ! l

  endfunction gamfc
  
  
endmodule MadelungCalculator_mod
