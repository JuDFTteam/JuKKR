#define CHECKALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif


!------------------------------------------------------------------------------
!> Set up modified index arrays for real space truncation.
!> This is a convenient workaround to make real space truncation possible.
!> Note: O(N**2) scaling in storage and setup time
module TruncationZone_mod
  implicit none
  private
  public :: TruncationZone, create, destroy

  type TruncationZone
    integer :: naez_trc = 0 !> number of atoms in truncation zone
    integer :: naez_all = 0 !> number of all atoms
    integer(kind=2), allocatable :: local_atom_idx(:) !> map atom index to truncation zone atom index (-1 == not in truncation zone == OUTSIDE)
    integer(kind=4), allocatable :: global_atom_id(:) !> map truncation zone atom index to atom index ("inverse" of local_atom_idx)
  endtype

  interface create
    module procedure createTruncationZone
  endinterface

  interface destroy
    module procedure destroyTruncationZone
  endinterface

  integer, parameter, private :: OUTSIDE = -1 ! should be -1 in C-language

  contains
  
#ifndef ell_int_t
#define ell_int_t integer(kind=1)
#endif

  subroutine createTruncationZone(self, mask, masks)
    type(TruncationZone), intent(inout) :: self
    ell_int_t, intent(in) :: mask(:)
    ell_int_t, intent(in), optional :: masks(:,:) !> one mask per local atom, dim(naez_all,num_local_atoms)

    integer :: ii, ind, memory_stat
    
    self%naez_all = size(mask)
    self%naez_trc = count(mask >= 0)

    ! setup index map from global indices to a process local view
    ALLOCATECHECK(self%local_atom_idx(OUTSIDE:self%naez_all))
    ALLOCATECHECK(self%global_atom_id(self%naez_trc)) ! setup "inverse" of local_atom_idx

    self%local_atom_idx(OUTSIDE:0) = OUTSIDE

    ind = 0
    do ii = 1, self%naez_all
      if (mask(ii) >= 0) then
        ind = ind + 1
        self%local_atom_idx(ii) = ind
        self%global_atom_id(ind) = ii ! quasi inverse
        ! inverse means that all(self%local_atom_idx(self%global_atom_id(:)) == [1, 2, 3, ..., naez_trc])
      else
        self%local_atom_idx(ii) = OUTSIDE ! atom not in truncation cluster
      endif
    enddo ! ii

    if (ind >= 2**15) stop 'integer(kind=2) not sufficient in TruncationZone_mod.F90!' ! value range exceeded 

    if (.not. present(masks)) return

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
