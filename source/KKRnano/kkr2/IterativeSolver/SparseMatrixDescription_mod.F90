module SparseMatrixDescription_mod
  implicit none
  private
  public :: SparseMatrixDescription, create, destroy, dump, getNNZ, getNrows
  
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
  subroutine createSparseMatrixDescription(self, blk_nrows, max_num_blocks)
    type(SparseMatrixDescription), intent(inout) :: self
    integer, intent(in) :: blk_nrows
    integer, intent(in) :: max_num_blocks

    allocate(self%ia(blk_nrows + 1))
    allocate(self%kvstr(blk_nrows + 1))
    allocate(self%ja(max_num_blocks))
    allocate(self%ka(max_num_blocks + 1))

    self%ia = 0
    self%kvstr = 0
    self%ja = 0
    self%ka = 0

    self%blk_nrows = blk_nrows
    self%max_blockdim = 0
    self%max_blocks_per_row = 0

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Returns number of non-zero elements (only if properly setup!).
  integer function getNNZ(self)
    type(SparseMatrixDescription), intent(in) :: self

    getNNZ = self%ka(self%ia(self%blk_nrows + 1)) - 1

  endfunction ! get

  !----------------------------------------------------------------------------
  !> Returns number of non-zero elements (only if properly setup!).
  integer function getNrows(self, naez)
    type(SparseMatrixDescription), intent(in) :: self
    integer, intent(in), optional :: naez

    if (present(naez)) then
      getNrows = self%kvstr(naez + 1) - 1 ! never happens, ToDo: remove
    else
      getNrows = self%kvstr(self%blk_nrows + 1) - 1
    endif
    
  endfunction ! get
  
  !----------------------------------------------------------------------------
  !> Destroys SparseMatrixDescription object.
  elemental subroutine destroySparseMatrixDescription(self)
    type(SparseMatrixDescription), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%ia, self%kvstr, self%ja, self%ka, stat=ist)

    self%blk_nrows = 0
    self%max_blockdim = 0
    self%max_blocks_per_row = 0
  endsubroutine ! destroy

  !---------------------------------------------------------------------------
  !> Writes SparseMatrixDescription to formatted file - useful for testing.
  subroutine dumpSparseMatrixDescription(self, filename)
    type(SparseMatrixDescription), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='formatted', action='write')
    write(fu, *) self%blk_nrows, size(self%ja)
    write(fu, *) self%kvstr
    write(fu, *) self%ia
    write(fu, *) self%ja
    write(fu, *) self%ka
    write(fu, *) self%max_blockdim
    write(fu, *) self%max_blocks_per_row
    close(fu)
  endsubroutine ! dump

  !---------------------------------------------------------------------------
  !> Creates and reads SparseMatrixDescription from formatted file
  !> - useful for testing.
  subroutine createSparseMatrixDescriptionFromFile(self, filename)
    type(SparseMatrixDescription), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: blk_nrows, max_num_blocks

    call destroy(self)
    
    open(fu, file=filename, form='formatted', action='read', status='old')
    read(fu, *)  blk_nrows, max_num_blocks
    call createSparseMatrixDescription(self, blk_nrows, max_num_blocks)
    read(fu, *) self%kvstr
    read(fu, *) self%ia
    read(fu, *) self%ja
    read(fu, *) self%ka
    read(fu, *) self%max_blockdim
    read(fu, *) self%max_blocks_per_row
    close(fu)
  endsubroutine ! create

endmodule ! SparseMatrixDescription_mod
