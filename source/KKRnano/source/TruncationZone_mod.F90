module TruncationZone_mod
!-------------------------------------------------------------------------------
!> Summary: Set up modified index arrays for real space truncation.
!> Author: Elias Rabel, Alexander R Thiess, Paul F Baumeister
!> Category: KKRnano, geometry
!>
!> This is a convenient workaround to make real space truncation possible.
!> Note: O(N**2) scaling in storage and setup time
!-------------------------------------------------------------------------------
#define CHECKALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#include "macros.h"
  implicit none
  private
  public :: TruncationZone, create, destroy

!!! indices into rbasis(1:3,:)
#define global_atom_id_t integer(kind=4)

!!! indices inside the truncation zone
#define trunc_index_t integer(kind=2)

  type TruncationZone
    integer :: naez_trc = 0 !> number of atoms in truncation zone
    integer :: naez_all = 0 !> number of all atoms in the system
    trunc_index_t, allocatable :: trunc_atom_idx(:) !> map atom index to truncation zone atom index (-1 == not in truncation zone == OUTSIDE)
    global_atom_id_t, allocatable :: global_atom_id(:) !> map truncation zone atom index to atom index ("inverse" of trunc_atom_idx)
  endtype

  interface create
    module procedure createTruncationZone
  endinterface

  interface destroy
    module procedure destroyTruncationZone
  endinterface

  integer, parameter, private :: OUTSIDE = -1 ! must be -1 in C-language, could be 0 in Fortran

#ifndef NDEBUG
  global_atom_id_t, allocatable, protected, public :: global_atom_id(:) ! see above. This array is accessible from everywhere for debug purposes
#endif

!!! ell quantum numbers
#ifndef ell_int_t
#define ell_int_t integer(kind=1)
#endif

  contains

  subroutine createTruncationZone(self, mask)
    type(TruncationZone), intent(inout) :: self
    ell_int_t, intent(in) :: mask(:)

    integer :: memory_stat
    global_atom_id_t :: gid
    trunc_index_t  :: idx
    
    self%naez_all = size(mask)
    if (self%naez_all > huge(gid)) stop 'integer-kind not sufficient for global_atom_id in TruncationZone_mod.F90!' ! value range exceeded
    
    self%naez_trc = count(mask >= 0)
    if (self%naez_trc > huge(idx)) stop 'integer-kind not sufficient for trunc_atom_idx in TruncationZone_mod.F90!' ! value range exceeded 

    ! setup index map from global indices to a process local view
    ALLOCATECHECK(self%trunc_atom_idx(OUTSIDE:self%naez_all))
    ALLOCATECHECK(self%global_atom_id(self%naez_trc)) ! setup "inverse" of trunc_atom_idx

    self%trunc_atom_idx(OUTSIDE:0) = OUTSIDE

    idx = 0
    do gid = 1, self%naez_all
      if (mask(gid) >= 0) then
        idx = idx + 1
        self%trunc_atom_idx(gid) = idx ! enumerate
        self%global_atom_id(idx) = gid ! quasi inverse ...
        ! inverse means that all(self%trunc_atom_idx(self%global_atom_id(:)) == [1, 2, 3, ..., naez_trc])
      else
        self%trunc_atom_idx(gid) = OUTSIDE ! atom not in truncation cluster
      endif
    enddo ! gid

#ifndef NDEBUG
    allocate(global_atom_id(self%naez_trc))
    global_atom_id(:) = self%global_atom_id(:) ! make a copy that we can use for DEBUG purposes by use TruncationZone_mod, only: global_atom_id
#endif

  endsubroutine ! create
  
  elemental subroutine destroyTruncationZone(self)
    type(TruncationZone), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%trunc_atom_idx, self%global_atom_id, stat=ist)
    self%naez_all = 0
    self%naez_trc = 0
  endsubroutine ! destroy

endmodule ! TruncationZone_mod
