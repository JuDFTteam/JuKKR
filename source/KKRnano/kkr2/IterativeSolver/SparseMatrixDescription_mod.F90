module SparseMatrixDescription_mod
  implicit none
  private
  public :: SparseMatrixDescription, create, destroy, dump
  public :: createSparseMatrixDescription, destroySparseMatrixDescription ! deprecated
  public :: dumpSparseMatrixDescription, createSparseMatrixDescriptionFromFile
  public :: getNNZ
  
  !> description of a (square) sparse matrix in VBR format
  !> VBR = variable block row.
  !>
  !> For a description of the VBR format see :
  !> Y. Saad, SPARSKIT: a basic tool kit for sparse matrix computations - Version 2 (1994).
  !>
  !
  !> The actual matrix data has to be stored separately
  !> in an 1D array.
  type SparseMatrixDescription
    !> block dimensions of block-rows and block-cols
    integer, allocatable :: kvstr(:)
    !> For each block row, give index in ja (and ka)
    integer, allocatable :: ia(:)
    !> Column-indices of non-zero blocks
    integer, allocatable :: ja(:)
    !> ia indices into ka - gives start indices of non-zero blocks
    !> in matrix data array
    integer, allocatable :: ka(:)
    !> number of block rows
    integer :: blk_nrows = 0
    !> maximal block dimension
    integer :: max_blockdim = 0
    !> Maximum number of non-zero blocks per block-row
    integer :: max_blocks_per_row = 0
  endtype

  interface create
    module procedure createSparseMatrixDescription, createSparseMatrixDescriptionFromFile
  endinterface
  
  interface destroy
    module procedure destroySparseMatrixDescription
  endinterface

  interface dump
    module procedure dumpSparseMatrixDescription
  endinterface

  contains

  !----------------------------------------------------------------------------
  !> Creates SparseMatrixDescription object.
  !
  !> Creates data structure that contains sparsity information of a
  !> square VBR (variable block row) matrix.
  subroutine createSparseMatrixDescription(sparse, blk_nrows, max_num_blocks)
    type(SparseMatrixDescription), intent(inout) :: sparse
    integer, intent(in) :: blk_nrows
    integer, intent(in) :: max_num_blocks

    allocate(sparse%ia(blk_nrows + 1))
    allocate(sparse%kvstr(blk_nrows + 1))
    allocate(sparse%ja(max_num_blocks))
    allocate(sparse%ka(max_num_blocks + 1))

    sparse%ia = 0
    sparse%kvstr = 0
    sparse%ja = 0
    sparse%ka = 0

    sparse%blk_nrows = blk_nrows
    sparse%max_blockdim = 0
    sparse%max_blocks_per_row = 0

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Returns number of non-zero elements (only if properly setup!).
  integer function getNNZ(sparse)
    type(SparseMatrixDescription), intent(in) :: sparse

    getNNZ = sparse%ka(sparse%ia(sparse%blk_nrows + 1)) - 1

  endfunction ! get

  !----------------------------------------------------------------------------
  !> Destroys SparseMatrixDescription object.
  elemental subroutine destroySparseMatrixDescription(sparse)
    type(SparseMatrixDescription), intent(inout) :: sparse

    integer :: ist ! ignore status
    deallocate(sparse%ia, stat=ist)
    deallocate(sparse%kvstr, stat=ist)
    deallocate(sparse%ja, stat=ist)
    deallocate(sparse%ka, stat=ist)

    sparse%blk_nrows = 0
    sparse%max_blockdim = 0
    sparse%max_blocks_per_row = 0
  endsubroutine ! destroy

  !---------------------------------------------------------------------------
  !> Writes SparseMatrixDescription to formatted file - useful for testing.
  subroutine dumpSparseMatrixDescription(sparse, filename)
    type(SparseMatrixDescription), intent(in) :: sparse
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='formatted', action='write')

    write(fu, *) sparse%blk_nrows, size(sparse%ja)
    write(fu, *) sparse%kvstr
    write(fu, *) sparse%ia
    write(fu, *) sparse%ja
    write(fu, *) sparse%ka
    write(fu, *) sparse%max_blockdim
    write(fu, *) sparse%max_blocks_per_row

    close(fu)
  endsubroutine ! dump

  !---------------------------------------------------------------------------
  !> Creates and reads SparseMatrixDescription from formatted file
  !> - useful for testing.
  subroutine createSparseMatrixDescriptionFromFile(sparse, filename)
    type(SparseMatrixDescription), intent(inout) :: sparse
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: blk_nrows, max_num_blocks

    call destroy(sparse)
    
    open(fu, file=filename, form='formatted', action='read', status='old')

    read(fu, *)  blk_nrows, max_num_blocks
    call createSparseMatrixDescription(sparse, blk_nrows, max_num_blocks)
    read(fu, *) sparse%kvstr
    read(fu, *) sparse%ia
    read(fu, *) sparse%ja
    read(fu, *) sparse%ka
    read(fu, *) sparse%max_blockdim
    read(fu, *) sparse%max_blocks_per_row

    close(fu)
  endsubroutine ! create

endmodule ! SparseMatrixDescription_mod
