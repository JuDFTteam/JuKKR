#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! JUST FOR TESTING purposes
! replace by proper implementation
module TEST_lcutoff_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private

  public :: initLcutoffNew
  
  integer, parameter, public :: N_radii = 4

  integer,                protected, public :: lm_low(N_radii)
  double precision,       protected, public :: cutoff_radius(N_radii)
  integer, allocatable,   protected, public :: lmarray(:)
  integer,                protected, public :: cutoffmode
  logical,                protected, public :: DEBUG_dump_matrix = .false.
  integer,                protected, public :: num_truncated(0:N_radii)

  contains

  !----------------------------------------------------------------------------
  subroutine initLcutoffNew(trunc_zone, atom_ids, arrays)
    use Main2Arrays_mod, only: Main2Arrays
    use TruncationZone_mod, only: TruncationZone, createTruncationZone

    type(TruncationZone), intent(inout) :: trunc_zone
    type(Main2Arrays), intent(in) :: arrays
    integer, intent(in) :: atom_ids(:) ! list of global atom IDs

    integer :: lmmaxd, atomindex, ilocal, ii, ind, ios, num_local_atoms, num, ist, nrad, irad
    integer, allocatable :: lmarray_temp(:), lmarray_full(:)

    lmmaxd = arrays%lmmaxd

    allocate(lmarray_full(size(arrays%rbasis, 2)))
    allocate(lmarray_temp(size(arrays%rbasis, 2)))

    nrad = 0
    cutoff_radius(:) = 9.d9 ! effectively infinity
    lm_low(:) = lmmaxd
    cutoffmode = 4 ! 3:iterative solver, 4:full solver
    
    open(91, file='lcutoff', form='formatted', action='read', status='old', iostat=ios)
    if (ios == 0) then ! opening was successful
      nrad = 1
      read(91,*) cutoff_radius(nrad)
      read(91,*) lm_low(nrad)
      read(91,*) cutoffmode

      lm_low(nrad+1:) = -1
      cutoff_radius(nrad+1:) = 0.d0
      do while (cutoffmode > 4)
        cutoffmode = cutoffmode - 2
        ! cutoff-mode 5: iterative solver with 2 cutoffs,
        !             6: direct solver with 2 cutoffs
        nrad = nrad+1
        read(91,*) cutoff_radius(nrad)
        read(91,*) lm_low(nrad)
      enddo ! while
      close(91)
    else
      write(6,*) 'No file "lcutoff" found, use defaults.' ! todo: convert to warning
    endif

    lmarray_full = 0 ! init
    num_local_atoms = size(atom_ids)
    
    if (num_local_atoms > 1 .and. nrad > 0) &
      warn(6, "cannot handle more than one local atom correctly with truncation, but found"+num_local_atoms)

    do ilocal = 1, num_local_atoms
      atomindex = atom_ids(ilocal) ! global atom index

      lmarray_temp(:) = lmmaxd ! init with the maximum

      do irad = 1, nrad
      
        ! truncation zones
        call calcCutoffarray(lmarray_temp, arrays%rbasis, arrays%rbasis(:,atomindex), arrays%bravais, cutoff_radius(irad), lm_low(irad))

      enddo ! irad
        
      lmarray_full = max(lmarray_full, lmarray_temp) ! reduction: merge truncation zones of local atoms
    enddo ! ilocal

    ! TODO: a bit confusing, is never deallocated
    call createTruncationZone(trunc_zone, lmarray_full)

    num_truncated(:) = 0
    num_truncated(0) = count(lmarray_full == lmmaxd)
    do irad = 1, nrad
      num_truncated(irad) = count(lmarray_full == lm_low(irad))
    enddo ! irad

    num = count(lmarray_full > 0)
    allocate(lmarray(num)) ! never deallocated - who cares

    ind = 0
    do ii = 1, size(lmarray_full)
      if (lmarray_full(ii) > 0) then
        ind = ind + 1
        lmarray(ind) = lmarray_full(ii)
        CHECKASSERT(trunc_zone%index_map(ii) == ind) ! TODO
      endif
    enddo ! ii

    deallocate(lmarray_temp, lmarray_full, stat=ist)
  endsubroutine ! initLcutoffNew

  
#if 0
!!!! never used !!

  !----------------------------------------------------------------------------
  !> Set all entries in gllh to zero which are out of range due to truncation
  subroutine cropGLLH(gllh, lmmaxd, naclsd, naezd, lmarray, numn0, indn0)
    double complex, intent(inout) :: gllh(lmmaxd,lmmaxd*naclsd,naezd)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: naclsd
    integer, intent(in) :: naezd
    integer, intent(in) :: lmarray(naezd)
    integer, intent(in) :: numn0(naezd)
    integer, intent(in) :: indn0(naezd,naclsd)

    integer ii, jj, lmmax1, lmmax2, lm1, lm2, clustersitelm
    double complex, parameter :: ZERO = (0.d0, 0.d0)

    do ii = 1, naezd
      do jj = 1, numn0(ii)

        lmmax1 = lmarray(ii)
        lmmax2 = lmarray(indn0(ii,jj))

        do lm2 = 1, lmmaxd
          do lm1 = 1, lmmaxd

            clustersitelm = lm2 + lmmaxd*(jj-1)

            if (lm1 > lmmax1 .or. lm2 > lmmax2)  gllh(lm1,clustersitelm,ii) = ZERO

          enddo ! lm1
        enddo ! lm2
        
      enddo ! jj
    enddo ! ii
    
  endsubroutine ! crop


  !------------------------------------------------------------------------------
  !> Generates matrix (\Delta T G_ref - 1) BUT PERFORMS L-CUTOFF BY TRUNCATING
  !> T-MATRIX
  !> on input: gllh contains G_ref, on output: gllh contains coefficient matrix
  subroutine generateCoeffMatrixCROPPED(gllh, numn0, indn0, tmatll, naez, lmmaxd, naclsd, lmarray)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: naclsd
    integer, intent(in) :: naez
    double complex, intent(inout) :: gllh(lmmaxd,lmmaxd*naclsd,naez)
    integer, intent(in) :: indn0(naez,naclsd)
    integer, intent(in) :: numn0(naez)
    double complex, intent(in) :: tmatll(lmmaxd,lmmaxd,naez)
    integer, intent(in) :: lmarray(:)

    double complex, parameter :: CONE = (1.d0,0.d0), ZERO = (0.d0,0.d0)
    double complex :: tgh(lmmaxd) ! temporary
    integer :: il1b, il2b
    integer :: lm1, lm2, lm3
    integer :: lmmax1, lmmax2, lmmax3
    integer :: site_index
    integer :: site_lm_index
    integer :: cluster_site_index
    integer :: cluster_site_lm_index

    ! -------------- Calculation of (Delta_t * G_ref - 1) ---------------
    !
    ! NUMN0(site_index) is the number of atoms in the reference cluster
    ! of atom/site 'site_index' (inequivalent atoms only!)
    ! INDN0 stores the index of the atom in the basis corresponding to
    ! the reference cluster atom (inequivalent atoms only!)
    ! -------------------------------------------------------------------

    !$omp parallel do private(site_index, site_lm_index, cluster_site_index, cluster_site_lm_index, il1b, il2b, lm1, lm2, lm3, tgh)
    do site_index = 1, naez
      il1b = lmmaxd*(site_index-1) ! offset
      do cluster_site_index = 1, numn0(site_index)

        lmmax2 = lmarray(indn0(site_index,cluster_site_index))
        lmmax1 = lmarray(site_index)
        lmmax3 = lmmax1
        il2b = lmmaxd*(indn0(site_index,cluster_site_index)-1) ! offset

        do lm2 = 1, lmmaxd
          cluster_site_lm_index = lm2 + lmmaxd*(cluster_site_index-1)

          tgh(:) = ZERO
          do lm3 = 1, lmmax3
            do lm1 = 1, lmmax1
              tgh(lm1) = tgh(lm1) + tmatll(lm1,lm3,site_index)*gllh(lm3,cluster_site_lm_index,site_index)
            enddo ! lm1
          enddo ! lm3

          do lm1 = 1, lmmaxd
            site_lm_index = il1b + lm1
            gllh(lm1,cluster_site_lm_index,site_index) = tgh(lm1)

            if (site_lm_index == il2b + lm2) then
              ! substract 1 only at the 'diagonal'
              gllh(lm1,cluster_site_lm_index,site_index) = gllh(lm1,cluster_site_lm_index,site_index) - CONE
            endif

            if (lm1 > lmmax1 .or. lm2 > lmmax2) gllh(lm1,cluster_site_lm_index,site_index) = ZERO

          enddo ! lm1
        enddo ! lm2

      enddo ! cluster_site_index
    enddo ! site_index
    !$omp endparallel do
  endsubroutine ! generateCoeffMatrixCROPPED

#endif  
  
  !----------------------------------------------------------------------------
  !> Modifies the array 'cutoffarray' on positions that correspond to sites that
  !> are further away from 'center' than 'dist_cut'. The value 'lm_low' is written
  !> at the modified positions.
  !>
  !> If merging = .true.: change the lm value only when it is larger than the
  !> original value -> this merges the truncation zones
  subroutine calcCutoffarray(cutoffarray, rbasis, center, bravais, dist_cut, lm_low)
    integer, intent(inout) :: cutoffarray(:)
    double precision, intent(in) :: rbasis(:,:) ! assumed(1:3,*)
    double precision, intent(in) :: center(3)
    double precision, intent(in) :: bravais(3,3)
    double precision, intent(in) :: dist_cut
    integer, intent(in) :: lm_low

    integer :: ii
    
    do ii = 1, size(rbasis, 2)
      if (distance2_pbc(rbasis(1:3,ii), center(1:3), bravais) > dist_cut*dist_cut) cutoffarray(ii) = lm_low
    enddo ! ii

  endsubroutine ! calc
  
  !> Calculate distance of two points taking into account the
  !> periodic boundary conditions
  !> Assume that points are in same unit cell
  !> Is this a valid assumption???
  double precision function distance2_pbc(point1, point2, bravais) result(dist_sq)
    double precision, intent(in) :: point1(3), point2(3)
    double precision, intent(in) :: bravais(3,3)

    double precision :: vec(3), vt(3)
    integer :: nx, ny, nz

    vec(1:3) = point2(1:3) - point1(1:3)

    dist_sq = vec(1)**2 + vec(2)**2 + vec(3)**2

    ! brute force distance checking
    do nx = -1, 1
      do ny = -1, 1
        do nz = -1, 1
          vt(:) = vec(:) + nx*bravais(:,1) + ny*bravais(:,2) + nz*bravais(:,3)
          dist_sq = min(dist_sq, vt(1)*vt(1) + vt(2)*vt(2) + vt(3)*vt(3))
        enddo ! nz
      enddo ! ny
    enddo ! nx

  endfunction ! distance squared

  
!!! never used 
!   subroutine getLMarray(lmarray, rbasis, center, bravais, dist_cut, lm_high, lm_low)
!     integer, intent(out) :: lmarray(:)
!     double precision, intent(in) :: rbasis(:,:) ! assumed(1:3,*)
!     double precision, intent(in) :: center(3)
!     double precision, intent(in) :: bravais(3,3)
!     double precision, intent(in) :: dist_cut
!     integer, intent(in):: lm_high, lm_low
! 
!     integer :: ii
!     double precision :: dist2
! 
!     do ii = 1, size(rbasis, 2)
!       dist2 = distance2_pbc(rbasis(1:3,ii), center, bravais)
!       if (dist2 > dist_cut*dist_cut) then
!         lmarray(ii) = lm_low
!       else
!         lmarray(ii) = lm_high
!       endif
!     enddo ! ii
! 
!   endsubroutine ! get
  
endmodule ! TEST_lcutoff_mod
