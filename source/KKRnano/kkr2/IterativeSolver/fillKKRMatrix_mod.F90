!> Build coefficient matrix for solution of Dyson equation.
!>
!> @author Elias Rabel

module fillKKRMatrix_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: getKKRMatrixStructure, buildKKRCoeffMatrix, buildRightHandSide, solveFull, convertToFullMatrix
  public :: dump

  double complex, parameter :: ZERO=(0.d0, 0.d0), CONE=(1.d0, 0.d0)
  
  interface dump
    module procedure dumpDenseMatrix, dumpSparseMatrixData
  endinterface

  contains

  !------------------------------------------------------------------------------
  !> Setup of the sparsity pattern of the KKR-Matrix.
  subroutine getKKRMatrixStructure(lmax_array, numn0, indn0, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    integer(kind=1), intent(in) :: lmax_array(:) !< lmax each row, dim(nRows)
    integer, intent(in) :: numn0(:) !< dim(nRows)
    integer(kind=2), intent(in) :: indn0(:,:) !< dim(maxval(numn0),nRows)
    type(SparseMatrixDescription), intent(inout) :: sparse

    integer :: nnzb, nRows, ncols, ij, irow, icol, jcol, lmmaxd

    nRows = size(lmax_array)
    lmmaxd = (maxval(lmax_array) + 1)**2
    ncols = nRows ! a logical square matrix

    assert( size(numn0) >= nRows )
    assert( size(indn0, 2) >= nRows )
    assert( size(sparse%ia) == 1 + nRows )

    nnzb = sum(numn0(1:nRows)) ! number of non-zero blocks

    assert( size(sparse%ja) >= nnzb )

    sparse%ia(:) = 0
    sparse%ja(:) = 0 ! init block number
    sparse%BlockDim = lmmaxd

    ij = 0
    do irow = 1, nRows
      sparse%ia(irow) = ij + 1 ! start indices into ja
      do icol = 1, numn0(irow)
        assert( ij < nnzb )

        assert( icol <= size(indn0, 1) )
        jcol = indn0(icol,irow)
        assert( 1 <= jcol .and. jcol <= ncols )

        ij = ij + 1
        sparse%ja(ij) = jcol

      enddo ! icol
    enddo ! irow
    sparse%ia(nRows + 1) = ij + 1 ! final, important since the ranges are always [ia(i) ... ia(i+1)-1]
    assert( ij == nnzb ) ! check

  endsubroutine ! getKKRMatrixStructure

  !---------------------------------------------------------------------
  !> Given G_ref - build (1 - TG_ref) -> coefficent matrix.
  !>
  !> @param smat  on input:  sparse block matrix containing Gref
  !>              on output: coefficient matrix (1 - TG_ref)
  !> @param sparse  sparse matrix description
  !> ia    for each row give index of first non-zero block in ja
  !> ja    column index array of non-zero blocks
  subroutine buildKKRCoeffMatrix(smat, tmatLL, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    double complex, intent(inout) :: smat(:,:,:) !< dim(lmmaxa,lmmaxd,nnzB)
    double complex, intent(in) :: tmatLL(:,:,:) !< dim(lmmaxd,lmmaxd,nRows)
    type(SparseMatrixDescription), intent(in) :: sparse

    double complex :: temp(size(smat, 1))
    integer :: lmmaxa, lmmaxd, nRows, nCols, iRow, jCol, lm2, lm3, Aind

    lmmaxa = size(smat, 1) ! may be memory-aligned
    lmmaxd = size(tmatLL, 2)
    assert( size(tmatLL, 1) == lmmaxd )
    assert( lmmaxa >= lmmaxd )
    nRows = sparse%nRows
    assert( nRows == size(tmatLL, 3) )
    nCols = nRows ! a logical square matrix

    do iRow = 1, nRows
      do Aind = sparse%ia(iRow), sparse%ia(iRow + 1) - 1
        jCol = sparse%ja(Aind) ! ja gives the block-column indices of non-zero blocks

#ifndef NDEBUG
        if (1 > iRow .or. iRow > nRows) then
          write (*,*) "buildKKRCoeffMatrix: invalid iRow", iRow
          stop
        endif

        if (1 > jCol .or. jCol > nCols) then
          write (*,*) "buildKKRCoeffMatrix: invalid jCol", jCol
          stop
        endif
#endif
        ! Note: naive truncation - truncate T matrix to square matrix with
        ! dimension lmmax1 * lmmax1
        ! maybe it would be better to calculate T-matrix with all
        ! smaller lmax

        ! something to think about:
        ! should Gref also be truncated to square?

        !    T (square mat.)          lmmax2                   lmmax2
        !  |----|                  |----------|             |----------|
        !  |    |  lmmax1     *    | G_ref    | lmmax3  =   |   T*G    |  lmmax1
        !  |----|                  |----------| (==lmmax1)  |----------|
        !  lmmax3==lmmax1
        
        do lm2 = 1, lmmaxd

          temp(:) = ZERO
          do lm3 = 1, lmmaxd
            temp(:) = temp(:) + tmatLL(:,lm3,iRow) * smat(lm3,lm2,Aind) ! T*G
          enddo ! lm3

          if (iRow == jCol) temp(lm2) = temp(lm2) - CONE ! subtract 1.0 from the diagonal

          smat(1:lmmaxd,lm2,Aind) = temp(1:lmmaxd)
        enddo ! lm2

      enddo ! Aind ! block columns
    enddo ! block rows

  endsubroutine ! buildKKRCoeffMatrix

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, lmmaxd, atom_indices, tmatLL)
    double complex, intent(out) :: mat_B(:,:,:) !> dim(lmmaxa,lmmaxd,nRHSs)
    integer, intent(in) :: lmmaxd ! ToDo: remove from interface and derive lmmaxd from size(tmatLL, 2)
    integer(kind=2), intent(in) :: atom_indices(:) ! truncation zone indices of the local atoms
    double complex, intent(in), optional :: tmatLL(:,:,:) !< dim(lmmaxd,lmmaxd,nRows)

    integer :: iRHS, atom_index, lm2

    mat_B = ZERO
    assert( size(tmatLL, 1) == lmmaxd )
    assert( size(tmatLL, 2) == lmmaxd )

    do iRHS = 1, size(atom_indices)
      atom_index = atom_indices(iRHS)

      do lm2 = 1, lmmaxd
        if (present(tmatLL)) then
          mat_B(:,lm2+lmmaxd*(iRHS - 1),atom_index) = tmatLL(:,lm2,atom_index)
        else
          mat_B(lm2,lm2+lmmaxd*(iRHS - 1),atom_index) = CONE ! set the block to a unity matrix
        endif
      enddo ! lm2

    enddo ! iRHS

  endsubroutine ! buildRightHandSide


  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertToFullMatrix(smat, ia, ja, BlockDim, full_A)
    double complex, intent(in) :: smat(:,:,:) !< dim(BlockDim,BlockDim,nnzb)
    integer, intent(in) :: ia(:) !> dim(nRows + 1)
    integer, intent(in) :: ja(:)
    integer, intent(in) :: BlockDim
    double complex, intent(out) :: full_A(:,:)

    integer :: iRow, jCol, nRows, Aind

    full_A = ZERO

    assert( size(smat, 1) == BlockDim )
    assert( size(smat, 2) == BlockDim )
    
    nRows = size(ia) - 1

    do iRow = 1, nRows
      do Aind = ia(iRow), ia(iRow+1) - 1
        jCol = ja(Aind)

        full_A(BlockDim*(iRow - 1) + 1:BlockDim*iRow,BlockDim*(jCol - 1) + 1:BlockDim*jCol) = smat(:,:,Aind)

      enddo ! Aind
    enddo ! iRow

  endsubroutine ! convertToFullMatrix

  !----------------------------------------------------------------------------
  !> Solution of a system of linear equations with multiple right hand sides,
  !> using standard dense matrix LAPACK routines.
  integer function solveFull(full_A, full_X) result(info)
    double complex, intent(inout) :: full_A(:,:)
    double complex, intent(inout) :: full_X(:,:) ! on entry this contains full_B, on exit the solution

    integer, allocatable :: ipvt(:)
    integer :: ndim, nRHSs
    external :: zgetrf, zgetrs ! LAPACK

    ndim  = size(full_A, 1)
    nRHSs = size(full_X, 2)
    assert( size(full_A, 2) == ndim ) ! must be square
    assert( size(full_X, 1) == ndim ) ! must match the dims of A
    
    allocate(ipvt(ndim))

    ! factorization
    call zgetrf(ndim, ndim, full_A, ndim, ipvt, info) ! LU-factorize
!   if (info /= 0) ! ToDo: warn

    call zgetrs('n', ndim, nRHSs, full_A, ndim, ipvt, full_X, ndim, info) ! solve the system of linear equations
!   if (info /= 0) ! ToDo: warn

    deallocate(ipvt, stat=info)
  endfunction ! solveFull
  
  
  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to a formatted or an unformatted file
  subroutine dumpSparseMatrixData(smat, filename, formatted)
    double complex, intent(in) :: smat(:,:,:)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: formatted

    if (formatted) then
      call dumpSparseMatrixDataFormatted(smat, filename)
    else
      call dumpSparseMatrixDataBinary(smat, filename)
    endif
  endsubroutine ! dump
  
  !----------------------------------------------------------------------------
  !> Write dense matrix to unformatted file
  subroutine dumpDenseMatrix(mat, filename, formatted)
    double complex, intent(in) :: mat(:,:)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: formatted

    if (formatted) then
      call dumpDenseMatrixFormatted(mat, filename)
    else
      call dumpDenseMatrixBinary(mat, filename)
    endif
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to unformatted file
  subroutine dumpSparseMatrixDataBinary(smat, filename)
    double complex, intent(in) :: smat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) smat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read sparse matrix data from unformatted file
  subroutine loadSparseMatrixData(smat, filename)
    double complex, intent(out) :: smat(:,:,:)
    character(len=*), intent(in) :: filename
    
    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) smat
    close(fu)
  endsubroutine ! load

  !----------------------------------------------------------------------------
  !> Write dense matrix to unformatted file
  subroutine dumpDenseMatrixBinary(mat, filename)
    double complex, intent(in) :: mat(:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) mat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read dense matrix data from unformatted file
  !> - useful for testing.
  subroutine loadDenseMatrix(mat, filename)
    double complex, intent(out) :: mat(:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) mat
    close(fu)
  endsubroutine ! load

  !----------------------------------------------------------------------------
  !> Write GLLh to unformatted file
  !> - useful for testing.
  subroutine dumpGLLhBinary(mat, filename)
    double complex, intent(in) :: mat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) mat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read GLLh from unformatted file
  !> - useful for testing.
  subroutine loadGLLh(mat, filename)
    double complex, intent(out) :: mat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) mat
    close(fu)
  endsubroutine ! load

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to formatted file
  !> - useful for testing.
  subroutine dumpSparseMatrixDataFormatted(smat, filename)
    double complex, intent(in) :: smat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: fi, si, bi

    open(fu, file=filename, form='formatted', action='write')
    do bi = 1, size(smat, 3) ! block index
      do si = 1, size(smat, 2) !  slow index
        do fi = 1, size(smat, 1) !  fast index
          write(fu, *) real(smat(fi,si,bi)), aimag(smat(fi,si,bi))
        enddo ! fi
      enddo ! si
    enddo ! bi
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Write dense matrix to formatted file
  !> First line gives matrix dimension: rows cols
  !> - useful for testing.
  subroutine dumpDenseMatrixFormatted(mat, filename)
    double complex, intent(in) :: mat(:,:)
    character(len=*), intent(in) :: filename
    
    integer, parameter :: fu = 97
    integer :: ii, jj

    open(fu, file=filename, form='formatted', action='write')
    write(fu, *) size(mat, 1), size(mat, 2)
    do jj = 1, size(mat, 2)
      do ii = 1, size(mat, 1)
         write(fu, *) real(mat(ii,jj)), aimag(mat(ii,jj))
      enddo ! ii
    enddo ! jj
    close(fu)
  endsubroutine ! dump

endmodule ! fillKKRMatrix_mod
