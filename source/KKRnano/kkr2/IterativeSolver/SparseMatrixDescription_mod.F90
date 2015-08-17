module SparseMatrixDescription_mod
  implicit none
  private
  public :: SparseMatrixDescription, create, destroy
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
    integer, allocatable, dimension(:) :: kvstr
    !> For each block row, give index in ja (and ka)
    integer, allocatable, dimension(:) :: ia
    !> Column-indices of non-zero blocks
    integer, allocatable, dimension(:) :: ja
    !> ia indices into ka - gives start indices of non-zero blocks
    !> in matrix data array
    integer, allocatable, dimension(:) :: ka
    !> number of block rows
    integer :: blk_nrows = 0
    !> maximal block dimension
    integer :: max_blockdim = 0
    !> Maximum number of non-zero blocks per block-row
    integer :: max_blocks_per_row = 0
  end type

  interface create
    module procedure createSparseMatrixDescription, createSparseMatrixDescriptionFromFile
  endinterface
  
  interface destroy
    module procedure destroySparseMatrixDescription
  endinterface
  
  CONTAINS

  !----------------------------------------------------------------------------
  !> Creates SparseMatrixDescription object.
  !
  !> Creates data structure that contains sparsity information of a
  !> square VBR (variable block row) matrix.
  subroutine createSparseMatrixDescription(sparse, blk_nrows, max_num_blocks)
    type (SparseMatrixDescription), intent(inout) :: sparse
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

  end subroutine

  !----------------------------------------------------------------------------
  !> Returns number of non-zero elements (only if properly setup!).
  integer function getNNZ(sparse)
    type (SparseMatrixDescription), intent(in) :: sparse

    getNNZ = sparse%ka(sparse%ia(sparse%blk_nrows + 1)) - 1

  end function

  !----------------------------------------------------------------------------
  !> Destroys SparseMatrixDescription object.
  subroutine destroySparseMatrixDescription(sparse)
    type (SparseMatrixDescription), intent(inout) :: sparse

    deallocate(sparse%ia)
    deallocate(sparse%kvstr)
    deallocate(sparse%ja)
    deallocate(sparse%ka)

    sparse%blk_nrows = 0
    sparse%max_blockdim = 0
    sparse%max_blocks_per_row = 0
  end subroutine

  !---------------------------------------------------------------------------
  !> Writes SparseMatrixDescription to formatted file - useful for testing.
  subroutine dumpSparseMatrixDescription(sparse, filename)
    type (SparseMatrixDescription), intent(inout) :: sparse
    character(len=*), intent(in) :: filename

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='formatted')

    write(FILEHANDLE, *) sparse%blk_nrows, size(sparse%ja)
    write(FILEHANDLE, *) sparse%kvstr
    write(FILEHANDLE, *) sparse%ia
    write(FILEHANDLE, *) sparse%ja
    write(FILEHANDLE, *) sparse%ka
    write(FILEHANDLE, *) sparse%max_blockdim
    write(FILEHANDLE, *) sparse%max_blocks_per_row

    close(FILEHANDLE)
  end subroutine

  !---------------------------------------------------------------------------
  !> Creates and reads SparseMatrixDescription from formatted file
  !> - useful for testing.
  subroutine createSparseMatrixDescriptionFromFile(sparse, filename)
    type (SparseMatrixDescription), intent(inout) :: sparse
    character(len=*), intent(in) :: filename

    integer, parameter :: FILEHANDLE = 97
    integer :: blk_nrows
    integer :: max_num_blocks

    open(FILEHANDLE, file=filename, form='formatted')

    read(FILEHANDLE, *)  blk_nrows, max_num_blocks
    call createSparseMatrixDescription(sparse, blk_nrows, max_num_blocks)
    read(FILEHANDLE, *) sparse%kvstr
    read(FILEHANDLE, *) sparse%ia
    read(FILEHANDLE, *) sparse%ja
    read(FILEHANDLE, *) sparse%ka
    read(FILEHANDLE, *) sparse%max_blockdim
    read(FILEHANDLE, *) sparse%max_blocks_per_row

    close(FILEHANDLE)
  end subroutine

end module SparseMatrixDescription_mod
