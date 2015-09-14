
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
  
  public :: testdimlat ! public for then --check functionality

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
!   type(MadelungCalculator), pointer :: madelung_calc
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
    integer :: lpot!, ist
    
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

  !----------------------------------------------------------------------------
  !> Creates data storage for a Madelung lattice sum.
  !>
  !> Has to be configured with a properly setup MadelungCalculator.
  !> The MadelungCalculator can be reused for several MadelungLatticeSums
  !> as long as the geometry and lmax does not change
  subroutine createMadelungLatticeSum(madelung_sum, lmxspd, num_atoms)!, madelung_calc
    type(MadelungLatticeSum), intent(inout) :: madelung_sum
    integer, intent(in) :: lmxspd, num_atoms

    madelung_sum%num_atoms = num_atoms
    allocate(madelung_sum%smat(lmxspd,num_atoms))
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
  subroutine calculateMadelungLatticeSum(self, calc, atom_index, rbasis)
    type(MadelungLatticeSum), intent(inout) :: self
    type(MadelungCalculator), intent(in) :: calc ! by introducing this, I try to get rid of the pointers in all the MadelungLatticeSum objects that all point to the same calculator
    integer, intent(in) :: atom_index
    double precision, intent(in) :: rbasis(3,self%num_atoms)

    call strmat(calc%alat, lmax=2*calc%lpot, naez=self%num_atoms, &
      ngmax=calc%lattice%ngmax, nrmax=calc%lattice%nrmax, &
      nlshellg=calc%lattice%nsg(calc%lattice%nshlg), nlshellr=calc%lattice%nsr(calc%lattice%nshlr), &
      gv=calc%lattice%gn, rv=calc%lattice%rm, qv=rbasis, vol=calc%volume0, i1=atom_index, &
      smat=self%smat) ! result

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
 

  ! call TESTDIMLAT(params%ALAT,arrays%BRAVAIS,RECBV,params%RMAX,params%GMAX, dims%NMAXD, dims%ISHLD) ! former call syntax
  subroutine testdimlat(alat, bravais, recbv, rmax, gmax)!, NMAXD, ISHLD)
    double precision, intent(in) :: alat !< lattice constant
    
    double precision, intent(in) :: gmax, rmax
    double precision, intent(in) :: bravais(3,3), recbv(3,3) !< bravais matrix and reciprocal basis vectors
    integer :: NMAXD, ISHLD!!!, intent(out)
    
    integer :: ngmax, nrmax, nshlg, nshlr
    double precision, allocatable :: vec(:,:) ! (3,nmaxd)
    integer, allocatable :: nsh(:)
    
    call lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsh, nsh, rmax, gmax, vec, vec, iprint=0, print_info=0)
    
    NMAXD = max(ngmax, nrmax)
    write(*,'(a,i13,9a)') ' NMAXD  =',NMAXD,'  ! minimum determined in ',__FILE__
    ISHLD = max(nshlg, nshlr)
    write(*,'(a,i13,9a)') ' ISHLD  =',ISHLD,'  ! minimum determined in ',__FILE__
    
  endsubroutine ! testdimlat
      
  
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
    
    integer :: i, n, l, k, ist, numr(3), numg(3), naiver, naiveg
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
    naiver = ceiling(rmax/absgm) ! Warning: cross term here: direct depends on reciprocal
    naiveg = ceiling(gmax/absrm) ! Warning: cross term here: reciprocal depends on direct
! write(0,*) trim(__FILE__+"loop"+(2*naiveg+1)+"in g-space from gmax ="+gmax)
! write(0,*) trim(__FILE__+"loop"+(2*naiver+1)+"in r-space from rmax ="+rmax)
    
    numr(1:3) = boxdim(bmat=br, dcut=rmax, naive=naiver, space='r')
    numg(1:3) = boxdim(bmat=bg, dcut=gmax, naive=naiveg, space='g')
! write(0,*) trim(__FILE__+"loop only"+(2*numg+1)+"in g-space from gmax ="+gmax)
! write(0,*) trim(__FILE__+"loop only"+(2*numr+1)+"in r-space from rmax ="+rmax)
    
  ! ! write(0,'(a)') trim(__FILE__+"alat"+alat)
  ! ! write(0,'(a)') trim(__FILE__+"bravais"+reshape(bravais, [9]))
  ! ! write(0,'(a)') trim(__FILE__+"br"+reshape(br, [9]))
  ! ! write(0,'(a)') trim(__FILE__+"recbv"+reshape(recbv, [9]))
  ! ! write(0,'(a)') trim(__FILE__+"bg"+reshape(bg, [9]))

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
  
  
! #define ORIGINAL_N_SQUARE_ALGORITHM    
  
  integer function count_vectors_in_sphere(numh, bm, dmax, vec, rad, nvecs, nshells, nsh, tol_origin, tol_newshell, space) result(ist)
#ifndef ORIGINAL_N_SQUARE_ALGORITHM
    use Sorting_mod, only: permutation_of
#endif    
    integer, intent(in) :: numh(3)
    double precision, intent(in) :: bm(3,3) !< bravais matrix or reciprocal space basis
    double precision, intent(in) :: dmax !< cutoff radius of the sphere
    double precision, allocatable, intent(out), optional :: vec(:,:), rad(:) ! check where this data is needed --> maybe better in one array (0:3,:)
    integer, intent(out) :: nvecs, nshells ! number
    integer, allocatable, intent(out), optional :: nsh(:)
    double precision, intent(in) :: tol_origin, tol_newshell
    character, intent(in) :: space ! for debug: 'r' or 'g'
    
    integer :: i, i1, i2, i3, ivmin, i01, ish!, num(3) 
    double precision :: dmax2, v2, da, db, vx(3), vxy(3), vxyz(3)
    double precision, allocatable :: cv(:,:), d2(:)
    integer, allocatable :: perm(:), nvis(:) ! tmp for nsh, nvis=number_of_vectors_in_shell
#ifdef  ORIGINAL_N_SQUARE_ALGORITHM
    double precision :: very_large
    integer :: iminl(1)
#endif

!     numh = num/2 + 1
!  ==> 1 - numh(1) : num(1) - numh(1) is equivalent to -num/2 : num - num/2, however, num was chosen odd, so its a symmetric
    
    dmax2 = dmax**2 ! radius^2
    
    ! **********************************************************************
    !                 generate lattice vectors of real or reciprocal space
    ! **********************************************************************

    do i01 = 0, 1 ! loop runs twice: iteration #0: count & allocate, iteration #1: store
    
      i = 0 ! init
      do i1 = -numh(1), numh(1) ! 1 - numh(1), num(1) - numh(1)
        vx(1:3) = i1*bm(1:3,1)
        do i2 = -numh(2), numh(2) ! 1 - numh(2), num(2) - numh(2)
          vxy(1:3) = i2*bm(1:3,2) + vx(1:3) 
          do i3 = -numh(3), numh(3) ! 1 - numh(3), num(3) - numh(3)
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
      
      if (i01 == 0) then ! after first iteration (the "count" iteration) and before the second iteration (the "store" iteration)
        nvecs = i
        allocate(cv(1:3,nvecs), d2(nvecs), nvis(nvecs), stat=ist)
      endif

    enddo ! i01    
    ! ======================================================================
!DBG  write(0,*) trim(__FILE__+"loop"+(2*numh+1)+" in"+space-"-space finds"+nvecs+"vectors")

    if (present(vec)) then; deallocate(vec, stat=ist); allocate(vec(1:3,nvecs), stat=ist); endif
    if (present(rad)) then; deallocate(rad, stat=ist); allocate(rad(nvecs), stat=ist); endif


    ! sort vectors in order of increasing absolute value
    ! warning: the original method of computing the shell structure scales N^2 with N=nvecs, better would be sorting with N*log(N)

#ifdef  ORIGINAL_N_SQUARE_ALGORITHM    
    very_large = dmax**2 + 9.d9
#else
    allocate(perm(1:nvecs), stat=ist)
    perm = permutation_of(d2(1:nvecs)) ! sort such that d2(perm(:)) is smallest first. This scales N*log(N)
#endif

    da = 0.d0 ! = tol_origin ! init ! todo: revise this (the first vector is [0.,0.,0.] always, so we can init with 0.d0)
                                    ! todo discuss why the origin tolerance was introduced althoug there is tol_newshell and (if redundant) remove it
    ish = 1 ! open the first shell
    nvis(ish) = 0 ! init number of vectors for the first shell
    
    do i = 1, nvecs
#ifdef  ORIGINAL_N_SQUARE_ALGORITHM    
      iminl = minloc(d2) ! find the location of the smallest element in d2 --> here, this algorithm scales Order(nvec^2)
      ivmin = iminl(1) ! pass index
#else
      ivmin = perm(i) ! use the index order found by quicksort
#endif

      db = sqrt(d2(ivmin))
      ! ----------------------------------------------------------------------
      if (db > da + tol_newshell) then
        
!DBG  write(0,*) trim(__FILE__+"create shell"+ish+"with"+nvis(ish)+"points in"+space-"-space") ! show the last shell
        
        ish = ish + 1 ! create a shell index of the current shell
        nvis(ish) = 0 ! init number of vectors for the new shell
        
        da = db
      endif ! does not lie in the open shell
      ! ----------------------------------------------------------------------
      
      nvis(ish) = nvis(ish) + 1 ! increase the number of points in this shell

      if (present(vec)) vec(1:3,i) = cv(1:3,ivmin) ! store vector
      if (present(rad)) rad(i) = db ! store radius
      
#ifdef  ORIGINAL_N_SQUARE_ALGORITHM    
      d2(ivmin) = very_large ! take it out of the loop
#endif      
    enddo ! i
    nshells = ish ! export the max number of shells

!DBG  write(0,*) trim(__FILE__+"create shell"+ish+"with"+nvis(ish)+"points in"+space-"-space (last)") ! show the last shell

    deallocate(cv, d2, perm, stat=ist) ! free local arrays
    
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
  !>*                   ylm(q(i) - q(j) - rv)                          *
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

subroutine strmat(alat, lmax, naez, ngmax, nrmax, nlshellg, nlshellr, gv, rv, qv, vol, i1, smat)
  use Harmonics_mod, only: ymy
  use Constants_mod, only: pi

    double precision, intent(in) :: alat
    double precision, intent(in) :: vol
    integer, intent(in) :: lmax !< usually 2*lpot
    integer, intent(in) :: naez !< number of all atoms 
    integer, intent(in) :: ngmax !< number of reciprocal space vectors
    integer, intent(in) :: nrmax !< number of real space vectors
    integer, intent(in) :: nlshellg !< =nsg(nshlg) number of vectors in the largest shell
    integer, intent(in) :: nlshellr !< =nsr(nshlr) number of vectors in the largest shell
    ! (todo for gv and rv: introduce a forth component that has the vector absolute, so the evaluation of ymy can be made faster)
    double precision, intent(in) :: gv(3,*) !< list of vectors in the reciprocal space 
    double precision, intent(in) :: rv(3,*) !< list of vectors in the real space 
    double precision, intent(in) :: qv(3,*) !< list of atomic positions (usuall called rbasis)
    integer, intent(in) :: i1 !< atom index of the source atom
    double precision, intent(out) :: smat(:,:) ! ((lmax+1)^2,naez)

    double complex, parameter :: CI=(0.d0,1.d0)
    double precision, parameter :: BOUND=1.d-8
    double precision :: alpha, lamda, kappa!, beta
    double precision :: dq(3), ga, vr(3), ra, ga2
    double precision :: dqdotg!, expbsq
    double precision :: fpi, rfac, sqrtPiInv, sqrtPi
    integer :: i, i2, i01
    integer :: l, m, lm
    integer :: nge, ngs, nre, nrs, nstart
    double complex :: bfac, stest((lmax+1)**2)
    double precision :: g(0:lmax), ylm((lmax+1)**2)

    fpi = 4.d0*pi
    sqrtPi = sqrt(pi)
    sqrtPiInv = 1./sqrtPi

    assert(size(smat, 1) == (lmax+1)**2) ! check leading dimension

    lamda = sqrt(pi)/alat ! choose proper splitting parameter
    kappa = -0.25d0/(lamda*lamda)

    ! **********************************************************************
    !$omp parallel do private(i2,dq,stest,lm,nstart,i01,nrs,ngs,nre,nge,i,vr,ylm,ra,alpha,g,rfac,l,m,vg,ga,ga2,dqdotg,bfac)
    do i2 = 1, naez
      !======================================================================
      dq(1:3) = (qv(1:3,i1) - qv(1:3,i2)) * alat

      stest(:) = 0.d0
      stest(1) = -sqrtPi/(vol*2.d0*lamda*lamda)
      
      ! --> exclude the origin and add correction if i1 == i2
      if (i1 == i2) then
        stest(1) = stest(1) - lamda/pi
        nstart = 2
      else
        nstart = 1
      endif
      
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! --> loop first over n-1 shells of real and reciprocal lattice - then
      !     add the contribution of the last shells to see convergence
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i01 = 0, 1
      
        if (i01 == 0) then
          nrs = nstart
          ngs = 2
          nre = nrmax - nlshellr ! nsr(nshlr)
          nge = ngmax - nlshellg ! nsg(nshlg)
        else
          nrs = nre + 1
          ngs = nge + 1
          nre = nrmax
          nge = ngmax
        endif

        ! --> sum over real lattice
        do i = nrs, nre
          vr(1:3) = dq(1:3) - rv(1:3,i)

          call ymy(vr(1), vr(2), vr(3), ra, ylm, lmax)
          alpha = lamda*ra
          g(0:lmax) = gamfc(lmax, alpha, ra)

          do l = 0, lmax
            rfac = g(l)*sqrtPiInv!/sqrt(pi)
            do m = -l, l
              lm = l*l + l + m + 1
              stest(lm) = stest(lm) + ylm(lm)*rfac
            enddo ! m
          enddo ! l
         
        enddo ! i

        ! --> sum over reciprocal lattice
        do i = ngs, nge

          call ymy(gv(1,i), gv(2,i), gv(3,i), ga, ylm, lmax)
          ga2 = ga*ga
          dqdotg = dq(1)*gv(1,i) + dq(2)*gv(2,i) + dq(3)*gv(3,i)

          bfac = fpi*exp(dcmplx(kappa*ga2, dqdotg))/(ga2*vol)
          do l = 0, lmax
            do m = -l, l
              lm = l*l + l + m + 1
              stest(lm) = stest(lm) + ylm(lm)*bfac
            enddo ! m
            bfac = (-CI)*bfac*ga/dble(2*l+1)
          enddo ! l
          
        enddo ! i

        if (i01 == 0) then
          if (any(abs(dimag(stest(:))) > BOUND)) die_here("Imaginary contribution to REAL lattice sum")
          smat(:,i2) = dble(stest(:)) ! store only real part
          stest(:) = 0.d0
        else
          smat(:,i2) = smat(:,i2) + dble(stest(:)) ! add only real part
!           ! --> test convergence
!           do lm = 1, (lmax+1)**2
!             !IF (2 < 1 .AND. ABS(dble(stest(lm))) > BOUND) WRITE (6,FMT=99001) I1,I2,LM,ABS(dble(stest(lm)))
! 99001       format (5x,'WARNING : Convergence of SMAT(',i2,',',i2,') for LMXSP =',i3,' is ',1p,d8.2,' > 1D-8',/,15x,'You should use more lattice vectors (RMAX/GMAX)')
!           enddo ! lm
        endif
        
      enddo ! i01
    enddo ! i2 ! loop over all atoms
    !$omp end parallel do
    
  endsubroutine strmat
  
  
  
  function boxdim(bmat, dcut, naive, space) result(ngv)
    integer          :: ngv(3) ! result
    double precision, intent(in) :: bmat(3,3), dcut
    integer , intent(in) :: naive
    character, intent(in) :: space

    integer :: i, j, k, i1, i2, i3
    double precision :: quad(3,3), rr(3,3), determ, wrap(3,3), scl(3,3), wGGw(3,3)

    quad = matmul(transpose(bmat), bmat) ! construct the metric such that [i1,i2,i3]*quad*[i1,i2,i3] = K2

    ! now: K^2 = i_1^2 Q_{11} + 2 i_1 i_2 Q_{12} + i_2^2 Q_{22}
    ! and in 3D:              + 2 i_2 i_3 Q_{23} + i_3^2 Q_{33} + 2 i_3 i_1 Q_{31}
    ! here, we have used that Q_{12} + Q_{21} = 2Q_{12} the metric is symmetric.
    ! for the k-th component i_k: derive the equation above w.r.t. the other two variables
    ! this is the perpendicular gradient. We demand it to vanish where the ellipsoid is at its maximum w.r.t. the coordinate i_k
    ! example: for i_1: derive equation above w.r.t. i_2 and i_3
    ! solve these equations finding i_2(i_1) and i_3(i_1)
    ! plug in i_2(i_1) and i_3(i_1) into the equation above and solve for i_1
       
    determ = invert3x3_r(quad, inverse=rr) ! the inverse solves the linear equations such that the perpendicular gradient is zero
    if (determ < 1d-12) then
      warn(6, "boxdim failed in"+space-"space")
      ngv(1:3) = naive
      return
    endif

    scl = 0. ; do i = 1, 3 ; scl(i,i) = 1./rr(i,i) ; enddo ! create a diagonal matrix
    
    wrap = matmul(rr, scl) ! wrap is unitless and has 1.0 in all three diagonal elements, however, the off-diagonal elements are important
    wGGw = matmul(wrap, matmul(quad, wrap))

    do i = 1, 3 ; ngv(i) = ceiling(dcut/sqrt(wGGw(i,i))) ; enddo ! i

  endfunction ! boxdim
  
  
  !! determinant of a 3 by 3 matrix
  double precision function determinant3x3_r(a) result(det)
    double precision, intent(in) :: a(3,3)
    det = a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
        + a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
        + a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
  endfunction ! determinant3x3

  !! inverts a 3 by 3 matrix, if the determinant is non-zero
  double precision function invert3x3_r(a, inverse) result(det) !! result: determinant of matrix a
    double precision, intent(in)                :: a(3,3) !! matrix
    double precision, intent(out)               :: inverse(3,3) !! inverse
    double precision :: invdet
    det = determinant3x3_r(a)
    if(det == 0.) then
#ifdef DEBUG
      if(o>0) write(o,'(9A)') sym,' invert3x3: determinant = 0.'
#endif
      inverse = 0.
      return
    endif ! det == 0
    invdet = 1./det
    inverse(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*invdet
    inverse(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))*invdet
    inverse(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*invdet

    inverse(1,2) = (a(3,2)*a(1,3)-a(3,3)*a(1,2))*invdet
    inverse(2,2) = (a(3,3)*a(1,1)-a(3,1)*a(1,3))*invdet
    inverse(3,2) = (a(3,1)*a(1,2)-a(3,2)*a(1,1))*invdet

    inverse(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*invdet
    inverse(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))*invdet
    inverse(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*invdet
  endfunction ! invert3x3
  
  
  
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
    !                                sum (2**i * x**(2i-1) / (2i-1)!!)
    !                              1..i..l
    !
    !
    ! Note: gamfc(alpha -> 0, ...) => glh = sqrt(pi) for all l  E.R.
    !-----------------------------------------------------------------------
    double precision :: glh(0:lmax) ! result
    integer, intent(in) :: lmax
    double precision, intent(in) :: alpha, r

    double precision, external :: erfcex
    double precision :: arg, facl, fex
    integer :: l

    arg = alpha*alpha
    facl = 2.d0*alpha

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
