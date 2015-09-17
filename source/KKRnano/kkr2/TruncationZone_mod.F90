#define CHECKALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if((STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif


!------------------------------------------------------------------------------
!> Set up modified index arrays for real space truncation.
!> This is a convinient workaround to make real space truncation possible.
!> Note: O(N**2) scaling in storage and setup time
module TruncationZone_mod
! use Main2Arrays_mod
  implicit none
  private
  public :: TruncationZone, create, destroy
  public :: createTruncationZone, destroyTruncationZone
  public :: translate, translateInd

  type TruncationZone
    integer :: naez_trc !< number of atoms in truncation zone
    integer :: naclsd
    !> map atom index to truncation zone atom index (-1 = not in trunc. zone)
    integer, allocatable :: index_map(:)
    !> map truncation zone atom index to atom index ("inverse" of index_map)
    integer, allocatable :: trunc2atom_index(:)
  endtype

  interface createTruncationZone
    module procedure createTruncationZoneOld, createTruncationZoneNew
  endinterface

  interface create
    module procedure createTruncationZoneNew
  endinterface

  interface destroy
    module procedure destroyTruncationZone
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  subroutine createTruncationZoneOld(self, mask, arrays)
    use Main2Arrays_mod, only: Main2Arrays
    type(TruncationZone), intent(inout) :: self
    integer, intent(in) :: mask(:)
    type(Main2Arrays), intent(in) :: arrays

    call createTruncationZoneNew(self, mask)
  endsubroutine ! create

  !----------------------------------------------------------------------------
  ! TODO: FIXME
  subroutine createTruncationZoneNew(self, mask)
    type(TruncationZone), intent(inout) :: self
    integer, intent(in) :: mask(:)

    integer :: num_atoms, ii, ind, naez_trc, memory_stat

    num_atoms = size(mask)

    ALLOCATECHECK(self%index_map(num_atoms))

    ind = 0
    do ii = 1, num_atoms
      if (mask(ii) > 0) then
        ind = ind + 1
        self%index_map(ii) = ind
      else
        self%index_map(ii) = -1 ! atom not in trunc. cluster
      endif
    enddo ! ii
    naez_trc = ind

    ! setup "inverse" of index_map
    ALLOCATECHECK(self%trunc2atom_index(naez_trc))
    ind = 0
    do ii = 1, num_atoms
      if (self%index_map(ii) > 0) then
        ind = ind + 1
        self%trunc2atom_index(ind) = ii
      endif
    enddo ! ii

    self%naez_trc = naez_trc

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyTruncationZone(self)
    type(TruncationZone), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%index_map)
    DEALLOCATECHECK(self%trunc2atom_index)

  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Translate atom index 'ind' to new atom index and return value
  elemental integer function translateInd(self, ind)
    type(TruncationZone), intent(in) :: self
    integer, intent(in) :: ind

    translateInd = -1
!   if (ind > ubound(self%index_map, 1)) stop 'TruncationZone: translateInd: ind too large!' ! not allowed in a pure function
    if (ind > 0) translateInd = self%index_map(ind)
  endfunction ! translateInd

  !----------------------------------------------------------------------------
  subroutine filter1d(mask, array, new_array)
    integer, intent(in) :: mask(:)
    integer, intent(in) :: array(:)
    integer, intent(inout) :: new_array(:)

    integer :: ii, ind

    ind = 0
    do ii = 1, size(mask)
      if (mask(ii) > 0) then
        ind = ind + 1
        new_array(ind) = array(ii)
      endif
    enddo ! ii
  endsubroutine

  !----------------------------------------------------------------------------
  subroutine filter2d1(mask, array, new_array)
    integer, intent(in) :: mask(:) 
    integer, intent(in) :: array(:,:)
    integer, intent(inout) :: new_array(:,:)

    integer :: ii, ind

    ind = 0
    do ii = 1, size(mask)
      if (mask(ii) > 0) then
        ind = ind + 1
        new_array(ind,:) = array(ii,:)
      endif
    enddo ! ii
  endsubroutine ! filter2d1

  !----------------------------------------------------------------------------
  subroutine filter2d2(mask, array, new_array)
    integer, intent(in) :: mask(:)
    integer, intent(in) :: array(:,:)
    integer, intent(inout) :: new_array(:,:)

    integer :: ii, ind

    ind = 0
    do ii = 1, size(mask)
      if (mask(ii) > 0) then
        ind = ind + 1
        new_array(:,ind) = array(:,ii)
      endif
    enddo ! ii
  endsubroutine ! filter2d2

  !----------------------------------------------------------------------------
  !> Translate atom indices in 'array' to new atom indices using 'index_map'
  subroutine translate(index_map, array)
    integer, intent(in) :: index_map(:)
    integer, intent(inout) :: array(:)

    integer :: ii

    do ii = 1, size(array)
      if (array(ii) > 0) array(ii) = index_map(array(ii))
    enddo ! ii
  endsubroutine ! translate

endmodule TruncationZone_mod
