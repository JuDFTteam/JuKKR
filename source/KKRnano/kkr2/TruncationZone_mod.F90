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
  public :: createTruncationZone!, destroyTruncationZone ! deprecated
! public :: translateInd

  type TruncationZone
    integer :: naez_trc = 0 !> number of atoms in truncation zone
    integer :: naez_all = 0 !> number of all atoms
    integer(kind=2), allocatable :: index_map(:) !> map atom index to truncation zone atom index (-1 = not in trunc. zone = OUTSIDE)
    integer(kind=4), allocatable :: trunc2atom_index(:) !> map truncation zone atom index to atom index ("inverse" of index_map)
  endtype

  public :: clear_non_existing_entries
  integer(kind=2), allocatable :: nzero_list(:), izero_list(:,:) !> NEW FEATURE: store
  
  interface create
    module procedure createTruncationZone
  endinterface

  interface destroy
    module procedure destroyTruncationZone
  endinterface

  integer, parameter, private :: OUTSIDE = 0 ! should be -1 in C-language
  
  contains

  !----------------------------------------------------------------------------
  ! TODO: FIXME
  subroutine createTruncationZone(self, mask, masks)
    type(TruncationZone), intent(inout) :: self
    integer(kind=1), intent(in) :: mask(:)
    integer(kind=1), intent(in), optional :: masks(:,:) !> one mask per local atom, dim(naez_all,num_local_atoms)

    integer :: ii, ind, num_local_atoms, il, iz, memory_stat

    self%naez_all = size(mask)
    self%naez_trc = count(mask > 0)

    ! setup index map from global indices to a process local view
    ALLOCATECHECK(self%index_map(OUTSIDE:self%naez_all))
    ALLOCATECHECK(self%trunc2atom_index(self%naez_trc)) ! setup "inverse" of index_map

    self%index_map(OUTSIDE:0) = OUTSIDE
    
    ind = 0
    do ii = 1, self%naez_all
      if (mask(ii) > 0) then
        ind = ind + 1
        self%index_map(ii) = ind
        self%trunc2atom_index(ind) = ii ! inverse means that all(self%index_map(self%trunc2atom_index(:)) == [1, 2, 3, ..., naez_trc])
      else
        self%index_map(ii) = OUTSIDE ! atom not in truncation cluster
      endif
    enddo ! ii
    
    if (ind >= 2**15) stop 'integer(kind=2) not sufficient in TruncationZone_mod.F90!' ! this is needed if

    if (.not. present(masks)) return
    
    ! this must hold: (mask(:) > 0) == any(masks > 0, dim=2)
    
#define SELF    
    num_local_atoms = size(masks, 2) ! number of local atoms
    allocate(SELF nzero_list(num_local_atoms))
    do il = 1, num_local_atoms
      SELF nzero_list(il) = count(masks(self%trunc2atom_index(:),il) <= 0)
    enddo ! il
    
    allocate(SELF izero_list(maxval(SELF nzero_list),num_local_atoms))

    do il = 1, num_local_atoms
      iz = 0
      do ind = 1, self%naez_trc
        ii = self%trunc2atom_index(ind)
        if (masks(ii,il) <= 0) then
          iz = iz + 1
          SELF izero_list(iz,il) = ind
        endif
      enddo ! ind
      if (iz /= SELF nzero_list(il)) stop __FILE__ ! fatal error
    enddo ! il
    
!   write(*,'(A,999(" ",i0))') 'Truncation zone created, Zero-List: ',SELF nzero_list
#undef SELF
  endsubroutine ! create
  
  subroutine clear_non_existing_entries(vec, N)
    double complex, intent(inout) :: vec(1:,1:)
    integer, intent(in), optional :: N !> block size
    integer :: il, iz, ic, BS
    
    if( .not. allocated(nzero_list)) return
    
    BS = (size(vec, 2) - 1)/size(nzero_list) + 1
    if (present(N)) BS = N
!!!!!!! WORKAROUND for truncation and num_local_atoms > 1, works only if all blocks have size N
    ! set elements in b to zero which only exist because some other of the local atoms is non-zero there
    ! shape(b) = [leaddim_b, ncols], ncols = N*num_local_atoms
    do il = 1, size(nzero_list) ! == num_local_atoms
      do iz = 1, nzero_list(il)
        ic = izero_list(iz,il)
        vec(ic*BS-BS+1:BS*ic,il*BS-BS+1:BS*il) = 0.d0
      enddo ! iz
    enddo ! il
!!!!!!! WORKAROUND
  endsubroutine
  
  
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

