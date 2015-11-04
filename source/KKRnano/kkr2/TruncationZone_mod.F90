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
! use Main2Arrays_mod
  implicit none
  private
  public :: TruncationZone, create, destroy
  public :: createTruncationZone!, destroyTruncationZone ! deprecated
! public :: translateInd

  type TruncationZone
    integer :: naez_trc = 0 !< number of atoms in truncation zone
    integer :: naez_all = 0 !< number of all atoms
    integer, allocatable :: index_map(:) !> map atom index to truncation zone atom index (-1 = not in trunc. zone = OUTSIDE)
    integer, allocatable :: trunc2atom_index(:) !> map truncation zone atom index to atom index ("inverse" of index_map)
  endtype

  interface create
    module procedure createTruncationZone
  endinterface

  interface destroy
    module procedure destroyTruncationZone
  endinterface

  integer, parameter, private :: OUTSIDE = -1
  
  contains

  !----------------------------------------------------------------------------
  ! TODO: FIXME
  subroutine createTruncationZone(self, mask)
    type(TruncationZone), intent(inout) :: self
    integer, intent(in) :: mask(:)

    integer :: naez_all, ii, ind, naez_trc, memory_stat

    naez_all = size(mask)

    ALLOCATECHECK(self%index_map(OUTSIDE:naez_all))

    self%index_map(OUTSIDE:0) = OUTSIDE
    
    ind = 0
    do ii = 1, naez_all
      if (mask(ii) > 0) then
        ind = ind + 1
        self%index_map(ii) = ind
      else
        self%index_map(ii) = OUTSIDE ! atom not in truncation cluster
      endif
    enddo ! ii
    naez_trc = ind

    ! setup "inverse" of index_map
    ALLOCATECHECK(self%trunc2atom_index(naez_trc))
    
    ind = 0
    do ii = 1, naez_all
      if (self%index_map(ii) > 0) then
        ind = ind + 1
        self%trunc2atom_index(ind) = ii
      endif
    enddo ! ii

    ! inverse means that all(self%index_map(self%trunc2atom_index(:)) == [1, 2, 3, ..., naez_trc])
    
    self%naez_trc = naez_trc
    self%naez_all = naez_all

  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyTruncationZone(self)
    type(TruncationZone), intent(inout) :: self

    integer :: ist ! ignore status

    deallocate(self%index_map, self%trunc2atom_index, stat=ist)
    self%naez_trc = 0
    
  endsubroutine ! destroy


endmodule ! TruncationZone_mod
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#if 0
!!! not used

  !----------------------------------------------------------------------------
  !> Translate atom index 'ind' to new atom index and return value
  elemental integer function translateInd(self, ind)
    type(TruncationZone), intent(in) :: self
    integer, intent(in) :: ind

!     translateInd = OUTSIDE
! !   if (ind > ubound(self%index_map, 1)) stop 'TruncationZone: translateInd: ind too large!' ! stop not allowed in a pure function
! 
!     if (ind > 0) &
      translateInd = self%index_map(ind)
  endfunction ! translateInd


  !----------------------------------------------------------------------------
  !> Translate atom indices in 'array' to new atom indices using 'index_map'
  subroutine translate(index_map, array)
    integer, intent(in) :: index_map(:)
    integer, intent(inout) :: array(:)

    where(array > 0)
      array = index_map(array)
    endwhere
    
!     integer :: ii
! 
!     do ii = 1, size(array)
!       if (array(ii) > 0) array(ii) = index_map(array(ii))
!     enddo ! ii
  endsubroutine ! translate
  
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
  endsubroutine ! filter1d

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
#endif

