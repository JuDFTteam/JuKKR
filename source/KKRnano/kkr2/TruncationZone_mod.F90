#define CHECKALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

!------------------------------------------------------------------------------
!> Set up modified index arrays for real space truncation.
!> This is a convenient workaround to make real space truncation possible.
!> Note: O(N**2) scaling in storage and setup time
module TruncationZone_mod
#include "macros.h"
  implicit none
  private
  public :: TruncationZone, create, destroy

  type TruncationZone
    integer :: naez_trc = 0 !> number of atoms in truncation zone
    integer :: naez_all = 0 !> number of all atoms in the system
    integer(kind=2), allocatable :: local_atom_idx(:) !> map atom index to truncation zone atom index (-1 == not in truncation zone == OUTSIDE)
    integer(kind=4), allocatable :: global_atom_id(:) !> map truncation zone atom index to atom index ("inverse" of local_atom_idx)
  endtype

  interface create
    module procedure createTruncationZone
  endinterface

  interface destroy
    module procedure destroyTruncationZone
  endinterface

  integer, parameter, private :: OUTSIDE = -1 ! must be -1 in C-language, could be 0 in Fortran

  contains
  
#ifndef ell_int_t
#define ell_int_t integer(kind=1)
#endif

  subroutine createTruncationZone(self, mask, masks)
    type(TruncationZone), intent(inout) :: self
    ell_int_t, intent(in) :: mask(:)
    ell_int_t, intent(in), optional :: masks(:,:) !> one mask per local atom, dim(naez_all,num_local_atoms)

    integer :: gid, idx, memory_stat
    integer(kind=8), parameter :: Two = 2
    
    self%naez_all = size(mask)
    if (self%naez_all >= Two**31) stop 'integer(kind=4) not sufficient for global_atom_id in TruncationZone_mod.F90!' ! value range exceeded
    
    self%naez_trc = count(mask >= 0)
    if (self%naez_trc >= Two**15) stop 'integer(kind=2) not sufficient for local_atom_idx in TruncationZone_mod.F90!' ! value range exceeded 

    ! setup index map from global indices to a process local view
    ALLOCATECHECK(self%local_atom_idx(OUTSIDE:self%naez_all))
    ALLOCATECHECK(self%global_atom_id(self%naez_trc)) ! setup "inverse" of local_atom_idx

    self%local_atom_idx(OUTSIDE:0) = OUTSIDE

    idx = 0
    do gid = 1, self%naez_all
      if (mask(gid) >= 0) then
        idx = idx + 1
        self%local_atom_idx(gid) = idx ! enumerate
        self%global_atom_id(idx) = gid ! quasi inverse
        ! inverse means that all(self%local_atom_idx(self%global_atom_id(:)) == [1, 2, 3, ..., naez_trc])
      else
        self%local_atom_idx(gid) = OUTSIDE ! atom not in truncation cluster
      endif
    enddo ! gid

    ! if (.not. present(masks)) return
    ! this must hold: (mask(:) >= 0) .eqv. any(masks >= 0, dim=2)
    !
    ! ToDo: use masks to determine more

  endsubroutine ! create
  
  elemental subroutine destroyTruncationZone(self)
    type(TruncationZone), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%local_atom_idx, self%global_atom_id, stat=ist)
    self%naez_all = 0
    self%naez_trc = 0
  endsubroutine ! destroy

endmodule ! TruncationZone_mod
