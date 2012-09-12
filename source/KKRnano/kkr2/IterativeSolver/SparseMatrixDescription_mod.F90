module SparseMatrixDescription_mod
  implicit none

  !> description of a (square) sparse matrix in VBR format
  !> VBR = variable block row
  type SparseMatrixDescription
    integer, allocatable, dimension(:) :: kvstr
    integer, allocatable, dimension(:) :: ia
    integer, allocatable, dimension(:) :: ja
    integer, allocatable, dimension(:) :: ka
    integer :: blk_nrows = 0
    integer :: max_blockdim = 0
    integer :: max_blocks_per_row = 0
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createSparseMatrixDescription(sparse, blk_nrows, max_num_blocks)
    implicit none
    type (SparseMatrixDescription), intent(inout) :: sparse
    integer, intent(in) :: blk_nrows
    integer, intent(in) :: max_num_blocks

    allocate(sparse%ia(blk_nrows + 1))
    allocate(sparse%kvstr(blk_nrows + 1))
    allocate(sparse%ja(max_num_blocks))
    allocate(sparse%ka(max_num_blocks + 1))

    sparse%blk_nrows = blk_nrows

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroySparseMatrixDescription(sparse)
    implicit none
    type (SparseMatrixDescription), intent(inout) :: sparse

    deallocate(sparse%ia)
    deallocate(sparse%kvstr)
    deallocate(sparse%ja)
    deallocate(sparse%ka)

    sparse%blk_nrows = 0
    sparse%max_blockdim = 0
    sparse%max_blocks_per_row = 0
  end subroutine

end module SparseMatrixDescription_mod
