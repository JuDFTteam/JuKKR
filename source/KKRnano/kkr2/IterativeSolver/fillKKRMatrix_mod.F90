!> Build coefficient matrix for solution of Dyson equation.
!>
!> @author Elias Rabel
#include "macros.h"

#ifndef NDEBUG
#define ASSERT(CONDITION) CHECKASSERT(CONDITION)
#else
#define ASSERT(CONDITION)
#endif

module fillKKRMatrix_mod
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

    integer(kind=1), intent(in) :: lmax_array(:) !< lmax each row, dim(nrows)
    integer, intent(in) :: numn0(:) !< dim(nrows)
    integer(kind=2), intent(in) :: indn0(:,:) !< dim(maxval(numn0),nrows)
    type(SparseMatrixDescription), intent(inout) :: sparse

    integer :: nnzb, nrows, ncols, ij, irow, icol, jcol, start_address, row_block, col_block, lmmaxd, lmax

    nrows = size(lmax_array)
    lmax = maxval(lmax_array)
    lmmaxd = (lmax + 1)**2
    ncols = nrows ! a logical square matrix

    ASSERT( size(numn0) >= nrows )
    ASSERT( size(indn0, 2) >= nrows )
    ASSERT( size(sparse%ia) == 1 + nrows )
    ASSERT( size(sparse%kvstr) == 1 + nrows )

    nnzb = sum(numn0(1:nrows)) ! number of non-zero blocks

    ASSERT( size(sparse%ja) >= nnzb )
    ASSERT( size(sparse%ka) >= 1 + nnzb )

    sparse%ia = 0
    sparse%ja = 0 ! init block number
    sparse%ka = 0 ! init block start address into data array
    sparse%kvstr = 0 ! init block start indices
    sparse%max_blockdim = 0
    sparse%max_blocks_per_row = 0

    start_address = 1
    sparse%kvstr(1) = 1
    ij = 1
    do irow = 1, nrows
      row_block = (lmax_array(irow) + 1)**2
      ASSERT( row_block == lmmaxd )
      sparse%kvstr(irow+1) = sparse%kvstr(irow) + row_block
      sparse%ia(irow) = ij ! start indices into ja
      sparse%max_blockdim = max(sparse%max_blockdim, row_block)
      ASSERT( sparse%kvstr(irow+1) == irow*sparse%max_blockdim + 1 ) ! the needs to hold when we want to take out kvstr
      sparse%max_blocks_per_row = max(sparse%max_blocks_per_row, numn0(irow))
      do icol = 1, numn0(irow)
        ASSERT( ij <= nnzb )
        sparse%ka(ij) = start_address ! start indices into the data array

        ASSERT( icol <= size(indn0, 1) )
        jcol = indn0(icol,irow)
        ASSERT( 1 <= jcol .and. jcol <= ncols )

        sparse%ja(ij) = jcol

        col_block = (lmax_array(jcol) + 1)**2
        ASSERT( col_block == lmmaxd )
        start_address = start_address + row_block*col_block

        ij = ij + 1
      enddo ! icol
    enddo ! irow
    sparse%ka(ij) = start_address
    ASSERT( sparse%ka(ij) == sparse%max_blockdim**2 * (ij - 1) + 1 ) ! the needs to hold when we want to take out ka
    sparse%ia(nrows+1) = ij ! final, important since the ranges are always [ia(i) ... ia(i+1)-1]
    ASSERT( ij == 1 + nnzb ) ! check

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

    double complex, intent(inout) :: smat(:)
    double complex, intent(in) :: tmatLL(:,:,:) !< dim(lmmaxd,lmmaxd,naez_trc)
    type(SparseMatrixDescription), intent(in) :: sparse

    double complex :: temp(size(tmatLL, 1))
    integer :: lmmaxd, naez_trc
    integer :: block_row, block_col, start, lm2, lm3, lmmax1, lmmax2, lmmax3
    integer :: istart_row, istop_row, istart_col, istop_col, ind_ia

    lmmaxd = size(tmatLL, 1)
    ASSERT( size(tmatLL, 2) == lmmaxd )
    naez_trc = size(tmatLL, 3)
    
    start = 0
    do block_row = 1, naez_trc
!     istart_row = sparse%kvstr(block_row)
!     istop_row  = sparse%kvstr(block_row+1)
      istart_row = sparse%max_blockdim*(block_row - 1) + 1
      istop_row  = sparse%max_blockdim*block_row + 1
      do ind_ia = sparse%ia(block_row), sparse%ia(block_row+1)-1

        block_col = sparse%ja(ind_ia) ! ja gives the block-column indices of non-zero blocks

!       istart_col = sparse%kvstr(block_col)
!       istop_col  = sparse%kvstr(block_col+1)
        istart_col = sparse%max_blockdim*(block_col - 1) + 1
        istop_col  = sparse%max_blockdim*block_col + 1

#ifndef NDEBUG
        if (1 > block_row .or. block_row > naez_trc) then
          write (*,*) "buildKKRCoeffMatrix: invalid block_row", block_row
          stop
        endif

        if (1 > block_col .or. block_col > naez_trc) then
          write (*,*) "buildKKRCoeffMatrix: invalid block_col", block_col
          stop
        endif
#endif
        ! Note: naive truncation - truncate T matrix to square matrix with
        ! dimension lmmax1 * lmmax1
        ! maybe it would be better to calculate T-matrix with all
        ! smaller lmax

        lmmax1 = istop_row - istart_row
        lmmax2 = istop_col - istart_col
        lmmax3 = lmmax1

        ! something to think about:
        ! should Gref also be truncated to square?

        !    T (square mat.)          lmmax2                   lmmax2
        !  |----|                  |----------|             |----------|
        !  |    |  lmmax1     *    | G_ref    | lmmax3  =   |   T*G    |  lmmax1
        !  |----|                  |----------| (==lmmax1)  |----------|
        !  lmmax3==lmmax1

        do lm2 = 1, lmmax2

          temp(:) = ZERO
          do lm3 = 1, lmmax3
            temp(:) = temp(:) + tmatLL(:,lm3,block_row) * smat(start+lm3+lmmax3*(lm2-1)) ! T*G
          enddo ! lm3
          
          if (block_row == block_col) temp(lm2) = temp(lm2) - CONE ! subtract 1.0 from the diagonal

          smat(start+lmmax1*(lm2-1)+ 1:lmmax1 +start+lmmax1*(lm2-1)) = temp(1:lmmax1)
        enddo ! lm2

        start = start + lmmax2*lmmax1
      enddo ! ind_ia ! block columns
    enddo ! block rows

  endsubroutine ! buildKKRCoeffMatrix

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, lmmaxd, atom_indices, kvstr, tmatLL)
    double complex, intent(out) :: mat_B(:,:)
    integer, intent(in) :: lmmaxd
    integer(kind=2), intent(in) :: atom_indices(:) ! truncation zone indices of the local atoms
    integer, intent(in) :: kvstr(:)
    double complex, intent(in), optional :: tmatLL(lmmaxd,lmmaxd,*) !< dim(lmmaxd,lmmaxd,naez_trc)

    integer :: start, ii, num_atoms, atom_index, lm2, lmmax1, lmmax2 

    mat_B = ZERO

    num_atoms = size(atom_indices)

    do ii = 1, num_atoms

      atom_index = atom_indices(ii)

!     lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)
      lmmax1 = lmmaxd

#ifndef NDEBUG
      if (lmmax1 /= lmmaxd) then
        write (*,*) "Central atom not treated with highest lmax", atom_index
        stop
      endif
#endif

      !lmmax2 = lmmaxd
      lmmax2 = lmmax1
      ! use naive truncation: lmmax1 = lmmax2
      ! Note: this is irrelevant, since the central atom
      ! should always be treated with the highest lmax

      start = (atom_index - 1)*lmmaxd ! = kvstr(atom_index) - 1
      do lm2 = 1, lmmax2
        if (present(tmatLL)) then
          mat_B(start+ 1:lmmax1 + start,lm2+lmmax2*(ii-1)) = tmatLL(1:lmmax1,lm2,atom_index)
        else
          mat_B(start+lm2,lm2+lmmax2*(ii-1)) = CONE ! set the block to a unity matrix
        endif
      enddo ! lm2

    enddo ! ii

  endsubroutine ! buildRightHandSide

!==============================================================================
! Routines for testing
!==============================================================================

  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertToFullMatrix(smat, ia, ja, ka, kvstr, kvstc, BlockDim, full)
    double complex, intent(in) :: smat(:)
    integer, intent(in) :: ia(:)
    integer, intent(in) :: ja(:)
    integer, intent(in) :: ka(:)
    integer, intent(in) :: kvstr(:), kvstc(:)
    integer, intent(in) :: BlockDim
    double complex, intent(out) :: full(:,:)

    integer :: ibrow, ibcol, irow, icol, nrows, ind, ind_ia

    full = ZERO

    nrows = size(ia) - 1

    do ibrow = 1, nrows
      ind = ka(ia(ibrow))
      do ind_ia = ia(ibrow), ia(ibrow+1) - 1

        ibcol = ja(ind_ia)

!       do   icol = kvstc(ibcol), kvstc(ibcol+1) - 1
!         do irow = kvstr(ibrow), kvstr(ibrow+1) - 1
      do   icol = BlockDim*(ibcol - 1) + 1, BlockDim*ibcol
        do irow = BlockDim*(ibrow - 1) + 1, BlockDim*ibrow

            full(irow,icol) = smat(ind)

            ind = ind + 1
          enddo ! irow
        enddo ! icol

      enddo ! ind_ia
    enddo ! ibrow

  endsubroutine ! convertToFullMatrix

  !----------------------------------------------------------------------------
  !> Solution of a system of linear equations with multiple right hand sides,
  !> using standard dense matrix LAPACK routines.
  subroutine solveFull(full, mat_B, mat_X)
    double complex, intent(inout) :: full(:,:)
    double complex, intent(in)  :: mat_B(:,:)
    double complex, intent(out) :: mat_X(:,:)

    integer, allocatable :: ipvt(:)
    integer :: ndim, nrhs, info
    external :: zgetrf, zgetrs ! LAPACK

    ndim = size(full, 1)
    allocate(ipvt(ndim))

    nrhs = size(mat_B, 2)
    mat_X(:,:) = mat_B(:,:)
    
    call zgetrf(ndim, ndim, full, ndim, ipvt, info) ! LU-factorize
!   if (info /= 0) ! ToDo: warn
    call zgetrs('n', ndim, nrhs, full, ndim, ipvt, mat_x, ndim, info) ! solve the system of linear equations
!   if (info /= 0) ! ToDo: warn

    deallocate(ipvt, stat=info)
  endsubroutine ! solveFull

  
  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to a formatted or an unformatted file
  subroutine dumpSparseMatrixData(smat, filename, formatted)
    double complex, intent(in) :: smat(:)
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
    double complex, intent(in) :: smat(:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) smat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read sparse matrix data from unformatted file
  subroutine loadSparseMatrixData(smat, filename)
    double complex, intent(out) :: smat(:)
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
    double complex, intent(in) :: smat(:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97
    integer :: ii

    open(fu, file=filename, form='formatted', action='write')

    do ii = 1, size(smat)
      write(fu, *) real(smat(ii)), aimag(smat(ii))
    enddo ! ii

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
