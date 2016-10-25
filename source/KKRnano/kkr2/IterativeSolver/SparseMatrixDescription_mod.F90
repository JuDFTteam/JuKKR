module SparseMatrixDescription_mod
  implicit none
  private
  public :: SparseMatrixDescription, create, destroy, dump
  
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
    !> block dimension
    integer :: BlockDim = 0
    !> number of block rows
    integer :: nRows = 0
    !> maximum column index
    integer :: nCols = 0
    !> number of non-zero blocks
    integer :: nnzb = 0
    !> For each block row, give start index in ja
    integer, allocatable :: ia(:) !> dim(nRows + 1)
    !> Column-indices of non-zero blocks
    integer, allocatable :: ja(:) !> dim(nnzb)
    !>
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
  subroutine createSparseMatrixDescription(self, nRows, nnzb, nCols)
    type(SparseMatrixDescription), intent(inout) :: self
    integer, intent(in) :: nRows !> number of rows
    integer, intent(in) :: nnzb !> number of non-zero blocks
    integer, intent(in), optional :: nCols !> number of columns

    allocate(self%ia(nRows + 1))
    allocate(self%ja(nnzb))

    self%ia(:) = 0
    self%ja(:) = 0

    self%nRows = nRows
    if (present(nCols)) then
      self%nCols = nCols
    else       
      self%nCols = nRows ! so far, this routine is only used for the KKR operator which is square
    endif
    self%nnzb = nnzb
    self%BlockDim = 0

  endsubroutine ! create
  
  !----------------------------------------------------------------------------
  !> Destroys SparseMatrixDescription object.
  elemental subroutine destroySparseMatrixDescription(self)
    type(SparseMatrixDescription), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%ia, self%ja, stat=ist)

    self%nRows = 0
    self%nCols = 0
    self%BlockDim = 0
    self%nnzb = 0
  endsubroutine ! destroy

  !---------------------------------------------------------------------------
  !> Writes SparseMatrixDescription to formatted file - useful for testing.
  subroutine dumpSparseMatrixDescription(self, filename)
    type(SparseMatrixDescription), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='formatted', action='write')
    write(fu, *) self%nRows, self%nCols, self%nnzb
    write(fu, *) self%ia
    write(fu, *) self%ja
    write(fu, *) self%BlockDim
    close(fu)
  endsubroutine ! dump

  !---------------------------------------------------------------------------
  !> Creates and reads SparseMatrixDescription from formatted file
  !> - useful for testing.
  subroutine createSparseMatrixDescriptionFromFile(self, filename)
    type(SparseMatrixDescription), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: nRows, nCols, nnzb

    call destroy(self)
    
    open(fu, file=filename, form='formatted', action='read', status='old')
    read(fu, *)  nRows, nCols, nnzb
    call createSparseMatrixDescription(self, nRows, nnzb, nCols)
    read(fu, *) self%ia
    read(fu, *) self%ja
    read(fu, *) self%BlockDim
    self%nnzb = nnzb
    close(fu)
  endsubroutine ! create

endmodule ! SparseMatrixDescription_mod
