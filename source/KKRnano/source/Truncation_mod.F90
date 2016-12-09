
! JUST FOR TESTING purposes
! replace by proper implementation
module Truncation_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
#include "DebugHelpers/logging_macros.h"
  use Logging_mod, only:    ! import no name, just mention it for the module dependency 
  USE_LOGGING_MOD
  implicit none
  private

  public :: initTruncation

#ifndef ell_int_t
#define ell_int_t integer(kind=1)
#endif
  
  ! union of the truncation zones over all local atoms
  ell_int_t, allocatable, protected, public :: lmax_a_array(:,:) ! ell_max for each atom-atom pair
  ell_int_t, allocatable, protected, public :: lmax_array(:) ! ell_max for each atom in the truncation region

  integer, allocatable, protected, public :: num_nonzero(:) ! dim(num_local_atoms)
  integer, allocatable, protected, public :: num_elements(:) ! dim(num_local_atoms)
  integer, protected, public :: num_truncated(-1:9)

  contains
  
#define useStatistics

  !----------------------------------------------------------------------------
  subroutine initTruncation(self, atom_ids, lmax, bravais, rbasis, lcutoff_radii, cutoff_radius, comm)
    use TruncationZone_mod, only: TruncationZone, create

#ifdef  useStatistics
    use Statistics_mod, only: SimpleStats, init, add, allreduce, eval
    type(SimpleStats) :: stats(1)
#endif

    type(TruncationZone), intent(inout) :: self
    integer, intent(in) :: lmax
    double precision, intent(in) :: bravais(3,3)
    double precision, intent(in) :: rbasis(:,:)
    integer, intent(in) :: atom_ids(:) !< list of global atom IDs
    double precision, intent(in) :: lcutoff_radii(0:) !< 
    double precision, intent(in) :: cutoff_radius !< 
    integer, intent(in) :: comm
    
    double precision, parameter :: R_active = 1.e-6 ! truncation radii below this are inactive
    integer :: naez, atomindex, ila, num_local_atoms, ist, nradii, ell
    ell_int_t, allocatable :: lmax_atom(:,:), lmax_full(:)
    ell_int_t :: l_lim(9), ellmax
    double precision :: r2lim(9)

    ellmax = lmax ! convert to ell_int_t

    l_lim = -1
    nradii = 0 ! 0:no truncation
    if (cutoff_radius > R_active) then ! a single cutoff radius is active
      ! this is equivalent to lcutoff_radii=[0 ... 0 cutoff_radius 0 ... 0] at position lmax 
      nradii = 1
      l_lim(nradii) = lmax ! convert to ell_int_t
      r2lim(nradii) = cutoff_radius**2
      if (any(lcutoff_radii > R_active)) &
        warn(6, "ell-dependent truncation (lcutoff_radii) is deactivated by cutoff_radius for ell="-lmax) 
    else
      ! now we check if the lcutoff_radii option has been specified in the input file
      do ell = min(ubound(lcutoff_radii, 1), lmax, 8), 0, -1
        if (lcutoff_radii(ell) > R_active) then
          nradii = nradii + 1 
          l_lim(nradii) = ell ! convert to ell_int_t
          r2lim(nradii) = lcutoff_radii(ell)**2
        endif
      enddo ! ell
    endif 
    l_lim = min(l_lim, ellmax) ! never higher than this

    num_local_atoms = size(atom_ids)
    naez = size(rbasis, 2)
    allocate(lmax_full(naez), lmax_atom(naez,num_local_atoms), stat=ist)
    if (ist /= 0) die_here("allocation of masks failed, requested"+(naez*.5**20*(num_local_atoms + 1))+"MiByte")
    
#ifdef  useStatistics    
    call init(stats, name=["ntrunc"])
#endif

    lmax_full = -1 ! init as truncated

    do ila = 1, num_local_atoms
      atomindex = atom_ids(ila) ! atom index into rbasis

      if (nradii > 0) then
        ! compute truncation zones
        call calcCutoffarray(lmax_atom(:,ila), l_lim, r2lim, rbasis, rbasis(:,atomindex), bravais, nradii)
        lmax_atom(atomindex,ila) = ellmax ! diagonal element is always full
      else
        lmax_atom(:,ila) = ellmax ! init as full interaction, no truncation
      endif

!! #define DEBUG_me      
#ifdef  DEBUG_me
      num_truncated(:) = 0
      do ell = -1, lmax
        num_truncated(ell) = count(lmax_atom(:,ila) == ell)
      enddo ! ell
      write(*,'(a,2(i0,a),9(" ",i0))') "atom #",atomindex,':  ',num_truncated(-1),' outside and inside s,p,d,f,... :',num_truncated(0:)
#endif

#ifdef  useStatistics    
      call add(stats(1), dble(count(lmax_atom(:,ila) >= 0)))
#endif

      lmax_full(:) = max(lmax_full(:), lmax_atom(:,ila)) ! reduction: merge truncation zones of local atoms
    enddo ! ila
    
#ifdef  useStatistics    
    ist = allreduce(stats, comm)
    WRITELOG(0,*) "truncation stats: ",trim(eval(stats(1)))
#endif

    num_truncated(:) = 0
    do ell = -1, lmax
      num_truncated(ell) = count(lmax_full == ell)
    enddo ! ell

    call create(self, mask=lmax_full)

    allocate(lmax_array(self%naez_trc)) ! lmax_array is never deallocated - who cares
    lmax_array(:) = lmax_full(self%global_atom_id(:)) ! compression to cluster atoms only

    ! also store the information for each local atom in module vars
    allocate(num_nonzero(num_local_atoms), num_elements(num_local_atoms))
    allocate(lmax_a_array(self%naez_trc,num_local_atoms)) ! todo: deallocate somewhen
    do ila = 1, num_local_atoms
      lmax_a_array(:,ila) = lmax_atom(self%global_atom_id(:),ila) ! compression to cluster atoms only
      num_nonzero(ila) = count(lmax_a_array(:,ila) >= 0)
      num_elements(ila) = sum((lmax_a_array(:,ila) + 1)**2)
    enddo ! ila

    if( .not. all(lmax_array == maxval(lmax_a_array, dim=2)) ) die_here('masks are inconsistent')

    deallocate(lmax_atom, lmax_full, stat=ist) ! ignore status
  endsubroutine ! initTruncation

  
  !----------------------------------------------------------------------------
  !> Modifies the array 'lmax_atom' on positions that correspond to sites that
  !> are further away from 'center' than 'dist_cut'. The value 'lm_low' is written
  !> at the modified positions.
  !>
  !> If merging = .true.: change the lm value only when it is larger than the
  !> original value -> this merges the truncation zones
  subroutine calcCutoffarray(lmax_atom, l_lim, radii2, rbasis, center, bravais, nradii)
    ell_int_t, intent(out) :: lmax_atom(:)
    ell_int_t, intent(in)  :: l_lim(nradii)
    double precision, intent(in) :: radii2(nradii)
    double precision, intent(in) :: rbasis(:,:) !> dim(1:3,*)
    double precision, intent(in) :: center(3)
    double precision, intent(in) :: bravais(3,3)
    integer, intent(in) :: nradii

    integer :: ii, il
    double precision :: d2

    !! warning, this operation is N^2 in the total number of atoms !!

!$OMP PARALLEL DO PRIVATE(ii, d2, il) schedule(static, 1)
    do ii = 1, size(rbasis, 2) ! loop over all atoms
      lmax_atom(ii) = -1 ! init as truncated
      d2 = distance2_pbc(rbasis(1:3,ii), center(1:3), bravais)

      !! we initialize lmax_atom with -1 and assume that there is at least one positive radius:
      do il = 1, nradii
        if (d2 <= radii2(il)) then
!!DEBUG     if(l_lim(il) > lmax_atom(ii)) write(*,'(2(a,i0),9(a,f0.3))') 'extend target #',ii,' to l=',l_lim(il),' rad(l)=',sqrt(radii2(il)),' distance=',sqrt(d2),' alat'
          lmax_atom(ii) = max(lmax_atom(ii), l_lim(il))
        endif
      enddo ! l

    enddo ! ii
!$OMP END PARALLEL DO

  endsubroutine ! calc
  
  !> Calculate distance of two points taking into account the periodic boundary conditions
  !> Assume that points are in same unit cell
  double precision function distance2_pbc(point1, point2, bravais) result(dist_sq)
    double precision, intent(in) :: point1(3), point2(3)
    double precision, intent(in) :: bravais(3,3)

    double precision :: vec(3), vt(3)
    double precision :: vtx(3), vtxy(3)
    integer :: nx, ny, nz

    vec(1:3) = point2(1:3) - point1(1:3)

    dist_sq = vec(1)**2 + vec(2)**2 + vec(3)**2 ! init with the distance without any periodic shift
    
    ! brute force distance checking in a periodic arrangement of cells
    do nx = -1, 1
      vtx(1:3) = vec(1:3) + nx*bravais(1:3,1)
      do ny = -1, 1
        vtxy(1:3) = vtx(1:3) + ny*bravais(1:3,2)
        do nz = -1, 1
!         vt(:) = vec(:) + nx*bravais(:,1) + ny*bravais(:,2) + nz*bravais(:,3) ! direct formula
          vt(1:3) = vtxy(1:3) + nz*bravais(1:3,3)
          dist_sq = min(dist_sq, vt(1)*vt(1) + vt(2)*vt(2) + vt(3)*vt(3))
        enddo ! nz
      enddo ! ny
    enddo ! nx

  endfunction ! distance squared

endmodule ! Truncation_mod
