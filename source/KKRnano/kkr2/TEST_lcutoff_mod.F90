#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! JUST FOR TESTING purposes
! replace by proper implementation
module TEST_lcutoff_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private

  public :: initLcutoffNew

#ifndef ell_int_t
#define ell_int_t integer(kind=1)
#endif
  
  ! union of the truncation zones over all local atoms
  ell_int_t, allocatable, protected, public :: lmax_array(:) ! ell_max
  ell_int_t, allocatable, protected, public :: lmax_a_array(:,:) ! ell_max for each atom-atom pair
  
  logical, protected, public :: DEBUG_dump_matrix = .false.
  integer, protected, public :: num_truncated(-1:8)
  integer, protected, public :: cutoffmode = 0

  contains

  !----------------------------------------------------------------------------
  subroutine initLcutoffNew(trunc_zone, atom_ids, arrays, lcutoff_radii, cutoff_radius, solver_type)
    use Main2Arrays_mod, only: Main2Arrays
    use TruncationZone_mod, only: TruncationZone, create

    type(TruncationZone), intent(inout) :: trunc_zone
    type(Main2Arrays), intent(in) :: arrays
    integer, intent(in) :: atom_ids(:) !< list of global atom IDs
    double precision, intent(in) :: lcutoff_radii(0:) !< 
    double precision, intent(in) :: cutoff_radius !< 
    integer, intent(in) :: solver_type !<
    
    double precision, parameter :: R_active = 1.e-6 ! truncation radii below this are inactive
    integer :: naez, lmax, atomindex, ila, num_local_atoms, ist, nradii, l
    ell_int_t, allocatable :: lmax_atom(:,:), lmax_full(:)
    ell_int_t :: l_lim(9), ellmax
    double precision :: r2lim(9)   
 
    cutoffmode = solver_type

    lmax = 0; do while ((lmax + 1)**2 < arrays%lmmaxd); lmax = lmax + 1; enddo ! find back global lmax
    ellmax = lmax  

    l_lim = -1
    nradii = 0 ! 0:no truncation
    if (cutoff_radius > R_active) then ! a single cutoff radius is active
      ! this is equivalent to lcutoff_radii=[0 ... 0 cutoff_radius 0 ... 0] at position lmax 
      nradii = 1
      l_lim(nradii) = lmax
      r2lim(nradii) = cutoff_radius**2
      if (any(lcutoff_radii > R_active)) &
        warn(6, "l-dependent truncation (lcutoff_radii) is deactivated by cutoff_radius for l="-lmax) 
    else
      ! now we check if the lcutoff_radii option has been specified in the input file
      do l = min(ubound(lcutoff_radii, 1), lmax, 8), 0, -1
        if (lcutoff_radii(l) > R_active) then
          nradii = nradii + 1 
          l_lim(nradii) = l
          r2lim(nradii) = lcutoff_radii(l)**2
        endif
      enddo ! l
    endif 

    num_local_atoms = size(atom_ids)
    naez = size(arrays%rbasis, 2)
    allocate(lmax_full(naez), lmax_atom(naez,num_local_atoms))
    
    lmax_full = -1 ! init as truncated

    do ila = 1, num_local_atoms
      atomindex = atom_ids(ila) ! atom index into rbasis

      if (nradii > 0) then
        lmax_atom(:,ila) = -1 ! init as truncated

        ! compute truncation zones
        call calcCutoffarray(lmax_atom(:,ila), arrays%rbasis, arrays%rbasis(:,atomindex), arrays%bravais, nradii, l_lim, r2lim)
        lmax_atom(:,ila) = min(lmax_atom(:,ila), ellmax)
        lmax_atom(atomindex,ila) = lmax ! diagonal element is full always
      else
        lmax_atom(:,ila) = lmax ! init as full interaction, no truncation
      endif
 
!! #define DEBUG_me      
#ifdef  DEBUG_me
      num_truncated(:) = 0
      do l = -1, lmax
        num_truncated(l) = count(lmax_atom(:,ila) == l)
      enddo ! l
      write(*,'(a,2(i0,a),9(" ",i0))') "atom #",atomindex,':  ',num_truncated(-1),' outside and inside s,p,d,f,... :',num_truncated(0:)
#endif

      lmax_full(:) = max(lmax_full(:), lmax_atom(:,ila)) ! reduction: merge truncation zones of local atoms
    enddo ! ila

    num_truncated(:) = 0
    do l = -1, lmax
      num_truncated(l) = count(lmax_full == l)
    enddo ! l

    if (num_local_atoms > 1 .and. num_truncated(lmax) /= naez) &
      warn(6, "cannot handle more than one local atom correctly with truncation, but found"+num_local_atoms)

    ! TODO: a bit confusing, is never deallocated
    call create(trunc_zone, mask=lmax_full, masks=lmax_atom)

    allocate(lmax_array(trunc_zone%naez_trc)) ! lmax_array is never deallocated - who cares
    lmax_array(:) = lmax_full(trunc_zone%trunc2atom_index(:)) ! compression to cluster atoms only

    ! also store the information for each local atom in module vars
    allocate(lmax_a_array(trunc_zone%naez_trc,num_local_atoms)) ! todo: deallocate somewhen
    do ila = 1, num_local_atoms
      lmax_a_array(:,ila) = lmax_atom(trunc_zone%trunc2atom_index(:),ila) ! compression to cluster atoms only
    enddo ! ila

    deallocate(lmax_atom, lmax_full, stat=ist)
  endsubroutine ! initLcutoffNew

  
  !----------------------------------------------------------------------------
  !> Modifies the array 'lmax_atom' on positions that correspond to sites that
  !> are further away from 'center' than 'dist_cut'. The value 'lm_low' is written
  !> at the modified positions.
  !>
  !> If merging = .true.: change the lm value only when it is larger than the
  !> original value -> this merges the truncation zones
  subroutine calcCutoffarray(lmax_atom, rbasis, center, bravais, nradii, l_lim, radii2)
    ell_int_t, intent(inout) :: lmax_atom(:)
    double precision, intent(in) :: rbasis(:,:) ! assumed(1:3,*)
    double precision, intent(in) :: center(3)
    double precision, intent(in) :: bravais(3,3)
    integer, intent(in) :: nradii
    ell_int_t, intent(in) :: l_lim(nradii)
    double precision, intent(in) :: radii2(nradii)

    integer :: ii, il
    double precision :: d2
    
    !! warning, this operation is N^2 in the total number of atoms !!

    do ii = 1, size(rbasis, 2) ! loop over all atoms
      d2 = distance2_pbc(rbasis(1:3,ii), center(1:3), bravais)

      !! we initialize lmax_atom with -1 and assume that there is at least one positive radius:
      do il = 1, nradii
        if (d2 <= radii2(il)) then
!!DEBUG     if(l_lim(il) > lmax_atom(ii)) write(*,'(2(a,i0),9(a,f0.3))') 'extend target #',ii,' to l=',l_lim(il),' rad(l)=',sqrt(radii2(il)),' distance=',sqrt(d2),' alat'
          lmax_atom(ii) = max(lmax_atom(ii), l_lim(il))
        endif
      enddo ! l

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

    ! brute force distance checking in a periodic arrangement of cells
    do nx = -1, 1
      do ny = -1, 1
        do nz = -1, 1
          vt(:) = vec(:) + nx*bravais(:,1) + ny*bravais(:,2) + nz*bravais(:,3)
          dist_sq = min(dist_sq, vt(1)*vt(1) + vt(2)*vt(2) + vt(3)*vt(3))
        enddo ! nz
      enddo ! ny
    enddo ! nx

  endfunction ! distance squared

endmodule ! TEST_lcutoff_mod
