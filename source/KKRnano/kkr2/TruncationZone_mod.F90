#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif


!------------------------------------------------------------------------------
!> Set up modified index arrays for real space truncation.
!> This is a convinient workaround to make real space truncation possible.
!> Note: O(N**2) scaling in storage and setup time
module TruncationZone_mod

  use Main2Arrays_mod

  private :: filter1d
  private :: filter2d1
  private :: filter2d2
  private :: translate

  type TruncationZone
    integer :: naez_trc
    integer :: naclsd
    integer, dimension(:), allocatable :: index_map
    integer, dimension(:), allocatable :: numn0_trc
    integer, dimension(:,:), allocatable :: indn0_trc
    integer, dimension(:,:), allocatable :: atom_trc
    integer, dimension(:,:), allocatable :: ezoa_trc
    integer, dimension(:), allocatable :: cls_trc
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createTruncationZone(self, mask, arrays)
    implicit none
    type (TruncationZone), intent(inout) :: self
    integer, dimension(:), intent(in) :: mask
    type (Main2Arrays), intent(in) :: arrays

    !-----
    integer :: num_atoms
    integer :: ii, jj
    integer :: ind
    integer :: ind_cls, mapped_index
    integer :: naez_trc
    integer :: naclsd
    integer :: memory_stat

    num_atoms = size(mask)
    naclsd = arrays%naclsd

    ALLOCATECHECK(self%index_map(num_atoms))

    ind = 0
    do ii = 1, num_atoms
      if (mask(ii) > 0) then
        ind = ind + 1
        self%index_map(ii) = ind
      else
        self%index_map(ii) = -1 ! atom not in trunc. cluster
      end if
    end do
    naez_trc = ind

    ALLOCATECHECK(self%numn0_trc(naez_trc))
    self%numn0_trc = -1
    call filter1d(mask, arrays%numn0, self%numn0_trc)

    ALLOCATECHECK(self%cls_trc(naez_trc))
    self%cls_trc = -1
    call filter1d(mask, arrays%cls, self%cls_trc)

    !!!ALLOCATECHECK(self%indn0_trc(naez_trc, maxval(self%numn0_trc)))
    ! overdimensioned for compatibility reasons
    ALLOCATECHECK(self%indn0_trc(naez_trc, arrays%naclsd))
    self%indn0_trc = -1

    ind = 0
    do ii = 1, num_atoms
      if (mask(ii) > 0) then
        ind = ind + 1
        ind_cls = 0
        do jj = 1, arrays%numn0(ii)
          mapped_index = self%index_map(arrays%indn0(ii, jj))
          if (mapped_index > 0) then
            ind_cls = ind_cls + 1
            self%indn0_trc(ind, ind_cls) = mapped_index
          end if
        end do
        self%numn0_trc(ind) = ind_cls
        ! reference clusters must be contained in truncation zone
        CHECKASSERT( self%numn0_trc(ind) == arrays%numn0(ii) )
      end if
    end do

    !write (*,*) "indn0     :", arrays%indn0
    !write (*,*) "indn0_trc :", self%indn0_trc
    !write (*,*) "lm        :", mask

    ALLOCATECHECK(self%atom_trc(arrays%naclsd, naez_trc))
    self%atom_trc = -1
    call filter2d2(mask, arrays%atom, self%atom_trc)

    ! atom_trc contains atom indices therefore translate
    ! to new atom indices
    do ii = 1, naez_trc
      call translate(self%index_map, self%atom_trc(:,ii))
    end do

    ! ezoa_trc contains indices of lattice vectors
    ! -> no translation necessary
    ALLOCATECHECK(self%ezoa_trc(arrays%naclsd, naez_trc))
    self%ezoa_trc = -1
    call filter2d2(mask, arrays%ezoa, self%ezoa_trc)

    self%naez_trc = naez_trc
    self%naclsd = arrays%naclsd
  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyTruncationZone(self)
    implicit none
    type (TruncationZone), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%index_map)
    DEALLOCATECHECK(self%numn0_trc)
    DEALLOCATECHECK(self%cls_trc)
    DEALLOCATECHECK(self%indn0_trc)
    DEALLOCATECHECK(self%atom_trc)
    DEALLOCATECHECK(self%ezoa_trc)

  end subroutine

  !----------------------------------------------------------------------------
  !> In-place reordering of a matrix array.
  !> Works only if index_map is monotonously increasing.
  !> Used to reorder T-matrix array
  subroutine reorderMatrices(self, mat_array)
    implicit none
    type (TruncationZone), intent(in) :: self
    double complex, dimension(:,:,:), intent(inout) :: mat_array

    !---------
    integer :: num, ii
    integer :: old_ind, ind

    num = size(self%index_map)
    CHECKASSERT(size(mat_array, 3) == num)

    old_ind = -1
    do ii = 1, num
      ind = self%index_map(ii)
      if (ind > 0 .and. ind /= ii) then
        CHECKASSERT( ind > old_ind .and. ind <= ii )
        mat_array(:,:,ind) = mat_array(:,:, ii)
      end if
      old_ind = ind
    end do

  end subroutine

  !----------------------------------------------------------------------------
  subroutine filter1d(mask, array, new_array)
    implicit none
    integer, dimension(:), intent(in) :: mask
    integer, dimension(:), intent(in) :: array
    integer, dimension(:), intent(inout) :: new_array

    integer :: ii, ind

    ind = 0
    do ii = 1, size(mask)
      if (mask(ii) > 0) then
        ind = ind + 1
        new_array(ind) = array(ii)
      end if
    end do
  end subroutine

  !----------------------------------------------------------------------------
  subroutine filter2d1(mask, array, new_array)
    implicit none
    integer, dimension(:), intent(in) :: mask
    integer, dimension(:,:), intent(in) :: array
    integer, dimension(:,:), intent(inout) :: new_array

    integer :: ii, ind

    ind = 0
    do ii = 1, size(mask)
      if (mask(ii) > 0) then
        ind = ind + 1
        new_array(ind, :) = array(ii, :)
      end if
    end do
  end subroutine

  !----------------------------------------------------------------------------
  subroutine filter2d2(mask, array, new_array)
    implicit none
    integer, dimension(:), intent(in) :: mask
    integer, dimension(:,:), intent(in) :: array
    integer, dimension(:,:), intent(inout) :: new_array

    integer :: ii, ind

    ind = 0
    do ii = 1, size(mask)
      if (mask(ii) > 0) then
        ind = ind + 1
        new_array(:, ind) = array(:, ii)
      end if
    end do
  end subroutine

  !----------------------------------------------------------------------------
  !> Translate atom indices in 'array' to new atom indices using 'index_map'
  subroutine translate(index_map, array)
    implicit none
    integer, dimension(:), intent(in) :: index_map
    integer, dimension(:), intent(inout) :: array

    integer :: ii
    do ii = 1, size(array)
      if (array(ii) > 0) then
        CHECKASSERT(array(ii) > 0)
        array(ii) = index_map(array(ii))
        CHECKASSERT(array(ii) > 0)
      end if
    end do
  end subroutine

end module TruncationZone_mod
