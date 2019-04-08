module SparseMatrixDescription_mod
  implicit none
  private
  public :: SparseMatrixDescription, create, destroy, dump, exists, subset
  
  !> description of a (square) block sparse matrix structure in Block Spare Row (BSR) format
  !
  !> The actual matrix data has to be stored separately in an 3D array.
  type SparseMatrixDescription
    !> number of block rows
    integer :: nRows = 0
    !> maximum column index
    integer :: nCols = 0
    !> number of non-zero blocks
    integer :: nnzb = 0
    !> For each block row, give start index in ColIndex
    integer, allocatable :: RowStart(:) !> dim(nRows + 1)
    !> Column-indices of non-zero blocks
    integer, allocatable :: ColIndex(:) !> dim(nnzb)
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

  integer function exists(bsr, row, col) result(Ind)
    type(SparseMatrixDescription), intent(in) :: bsr
    integer, intent(in) :: row, col

    Ind = BSR_entry_exists(bsr%RowStart, bsr%Colindex, row, col)
    
  endfunction ! exists

  integer function BSR_entry_exists(RowStart, Colindex, row, col) result(Ind)
    integer, intent(in) :: RowStart(:), ColIndex(:), row, col

    ! ToDo: the ColIndex list should be sorted ascendingly, so bisection search will be faster
    do Ind = RowStart(row), RowStart(row + 1) - 1
      if (ColIndex(Ind) == col) return ! Ind
    enddo ! Ind
    Ind = -1 ! not found

  endfunction ! exists

  integer function subset(set, sub, list) result(nfail) ! returns 0 un success
    !! this is needed to subtract or add the values of two operators whereas SUB must be a subset of SET
    type(SparseMatrixDescription), intent(in) :: set, sub
    integer, intent(out)                      :: list(:) !< dim(set%nnzb), list of indices into SET for which SUB exists
    ! local vars
    integer :: iRow, jCol, SubInd, SetInd
    
    nfail = 0 ! init
    if (size(list) < sub%nnzb) then
#ifdef DEBUG
      write(*,*) __FILE__,": BSR subset: ERROR: ",size(list)," = size(list) < sub%nnzb = ",sub%nnzb
#endif
      nfail = -sub%nnzb
      return
    endif ! list too short
    
    do iRow = 1, sub%nRows
      do SubInd = sub%RowStart(iRow), sub%RowStart(iRow + 1) - 1
        jCol = sub%ColIndex(SubInd)
        SetInd = exists(set, iRow, jCol)
        if (SetInd > 0) then
          list(SubInd) = SetInd ! entry exists
        else
          nfail = nfail + 1
#ifdef DEBUG
          write(*,*) __FILE__,": BSR subset: ERROR: entry(iRow, jCol) of subset does not exist in set. iRow, jCol =",iRow, jCol
#endif
        endif
      enddo ! SubInd
    enddo ! iRow
  endfunction ! subset

  !----------------------------------------------------------------------------
  !> Creates SparseMatrixDescription object.
  !
  !> Creates data structure that contains sparsity information of a
  !> square BSR (block sparse row) matrix.
  subroutine createSparseMatrixDescription(self, nRows, nnzb, nCols)
    type(SparseMatrixDescription), intent(inout) :: self
    integer, intent(in) :: nRows !> number of rows
    integer, intent(in) :: nnzb !> number of non-zero blocks
    integer, intent(in), optional :: nCols !> number of columns

    self%nRows = nRows
    allocate(self%RowStart(nRows + 1))
    self%RowStart(:) = 0
    
    self%nnzb = nnzb
    allocate(self%ColIndex(nnzb))
    self%ColIndex(:) = 0

    ! so far, this routine is only used for the KKR operator which is square
    self%nCols = nRows ; if (present(nCols)) self%nCols = nCols

  endsubroutine ! create
  
  !----------------------------------------------------------------------------
  !> Destroys SparseMatrixDescription object.
  elemental subroutine destroySparseMatrixDescription(self)
    type(SparseMatrixDescription), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%RowStart, self%ColIndex, stat=ist)

    self%nRows = 0
    self%nCols = 0
    self%nnzb = 0
  endsubroutine ! destroy

  !---------------------------------------------------------------------------
  !> Writes SparseMatrixDescription to formatted file
  subroutine dumpSparseMatrixDescription(self, filename)
    type(SparseMatrixDescription), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='formatted', action='write')
    write(fu, *) self%nRows, self%nCols, self%nnzb
    write(fu, '(9999(i0," "))') self%RowStart
    write(fu, '(9999(i0," "))') self%ColIndex
    close(fu)
  endsubroutine ! dump

  !---------------------------------------------------------------------------
  !> Creates and reads SparseMatrixDescription from formatted file
  subroutine createSparseMatrixDescriptionFromFile(self, filename)
    type(SparseMatrixDescription), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: nRows, nCols, nnzb

    call destroy(self)
    open(fu, file=filename, form='formatted', action='read', status='old')
    read(fu, *)  nRows, nCols, nnzb
    call createSparseMatrixDescription(self, nRows, nnzb, nCols)
    read(fu, *) self%RowStart
    read(fu, *) self%ColIndex
    close(fu)
  endsubroutine ! create

endmodule ! SparseMatrixDescription_mod
