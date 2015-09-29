!> Build coefficient matrix for solution of Dyson equation.
!>
!> @author Elias Rabel

#ifndef NDEBUG
#define ASSERT(CONDITION) if (.not. (CONDITION)) then; write(*,*) "Assertion ", #CONDITION, " failed: ",  __FILE__, __LINE__; endif
#else
#define ASSERT(CONDITION)
#endif

module fillKKRMatrix_mod
  implicit none
  private
  public :: getKKRMatrixStructure, buildKKRCoeffMatrix, buildRightHandSide, solveFull, convertToFullMatrix
  public :: dumpDenseMatrix, dumpDenseMatrixFormatted, dumpSparseMatrixData, dumpSparseMatrixDataFormatted

  double complex, parameter :: ZERO=(0.d0, 0.d0), CONE=(1.d0, 0.d0)
  
  contains

!------------------------------------------------------------------------------
!> Setup of the sparsity pattern of the KKR-Matrix.
  subroutine getKKRMatrixStructure(lmmaxd_array, numn0, indn0, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    
    integer, intent(in) :: lmmaxd_array(:)
    integer, intent(in) :: numn0(:)
    integer, intent(in) :: indn0(:,:)
    type(SparseMatrixDescription), intent(inout) :: sparse

    integer :: nnz_blocks, nrows, ii, irow, icol, start_address

    nrows = size(lmmaxd_array)

    ASSERT(size(numn0) == nrows)
    ASSERT(size(indn0, 1) == nrows)
    ASSERT(size(sparse%ia) == nrows + 1)
    ASSERT(size(sparse%kvstr) == nrows + 1)

    sparse%ia = 0
    sparse%kvstr = 0
    sparse%ja = 0
    sparse%ka = 0

    nnz_blocks = 0
    do ii = 1, nrows
      nnz_blocks = nnz_blocks + numn0(ii)
    enddo ! ii

    ASSERT(nnz_blocks <= size(sparse%ja))
    ASSERT(nnz_blocks + 1 <= size(sparse%ka))

    sparse%kvstr(1) = 1
    do ii = 2, nrows + 1
      sparse%kvstr(ii) = sparse%kvstr(ii-1) + lmmaxd_array(ii - 1)
    enddo ! ii

    ii = 1
    do irow = 1, nrows
      sparse%ia(irow) = ii
      do icol = 1, numn0(irow)  ! square matrix
        ASSERT(icol <= size(indn0, 2))
        ASSERT(ii <= nnz_blocks)
        sparse%ja(ii) = indn0(irow, icol)
        ii = ii + 1
      enddo ! icol
    enddo ! irow
    sparse%ia(nrows + 1) = ii

    start_address = 1
    ii = 1
    do irow = 1, nrows
      do icol = 1, numn0(irow)
        sparse%ka(ii) = start_address

        ASSERT( indn0(irow, icol) >= 1 .and. indn0(irow, icol) <= nrows )

        start_address = start_address + lmmaxd_array(irow)*lmmaxd_array(indn0(irow,icol))
        ii = ii + 1
      enddo ! icol
    enddo ! irow
    sparse%ka(ii) = start_address

    sparse%max_blockdim = maxval(lmmaxd_array)
    sparse%max_blocks_per_row = maxval(numn0)

  endsubroutine ! getKKRMatrixStructure

  !---------------------------------------------------------------------
  !> Given G_ref - build (1 - TG_ref) -> coefficent matrix.
  !>
  !> @param smat  on input:  sparse block matrix containing Gref
  !>              on output: coefficient matrix (1 - TG_ref)
  !> @param sparse  sparse matrix description
  !> ia    for each row give index of first non-zero block in ja
  !> ja    column index array of non-zero blocks
  subroutine buildKKRCoeffMatrix(smat, TMATLL, lmmaxd, num_atoms, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    double complex, intent(inout) :: smat(:)
    double complex, intent(in) :: TMATLL(lmmaxd,lmmaxd,num_atoms)
    integer, intent(in) :: lmmaxd, num_atoms
    type(SparseMatrixDescription), intent(inout) :: sparse

    double complex :: temp(lmmaxd)
    integer :: block_row, block_col, start, lm1, lm2, lm3, lmmax1, lmmax2, lmmax3
    integer :: istart_row, istop_row, istart_col, istop_col, ind_ia

    start = 0
    do block_row = 1, num_atoms
      istart_row = sparse%kvstr(block_row)
      istop_row  = sparse%kvstr(block_row+1)
      do ind_ia = sparse%ia(block_row), sparse%ia(block_row+1)-1

        block_col = sparse%ja(ind_ia)  !ja gives the block-column indices of non-zero blocks

        istart_col = sparse%kvstr(block_col)
        istop_col  = sparse%kvstr(block_col+1)

#ifndef NDEBUG
        if (block_row < 1 .or. block_row > num_atoms) then
          write (*,*) "buildKKRCoeffMatrix: invalid block_row", block_row
          STOP
        endif

        if (block_col < 1 .or. block_col > num_atoms) then
          write (*,*) "buildKKRCoeffMatrix: invalid block_col", block_col
          STOP
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
        ! -|    |  lmmax1     *    | G_ref    | lmmax3  =   |  -T*G    |  lmmax1
        !  |----|                  |----------| (=lmmax1)   |----------|
        !  lmmax3=lmmax1

        do lm2 = 1, lmmax2

          temp(:) = ZERO
          do lm1 = 1, lmmax1
            do lm3 = 1, lmmax3
              temp(lm1) = temp(lm1) - TMATLL(lm1,lm3,block_row) * smat(start+(lm2-1)*lmmax3+lm3)  ! -T*G
            enddo ! lm3
          enddo ! lm1

          do lm1 = 1, lmmax1

            if (block_row == block_col .and. lm1 == lm2) temp(lm1) = temp(lm1) + CONE  ! add 1 on the diagonal

            smat(start+(lm2-1)*lmmax1+lm1) = temp(lm1)
            !smat(start + (lm2 - 1) * lmmax1 + lm1) = (block_row * 100.d0 + block_col) * CONE
          enddo ! lm1

        enddo ! lm2

        start = start + lmmax2*lmmax1
      enddo ! ind_ia ! block columns
    enddo ! block rows

  endsubroutine ! buildKKRCoeffMatrix

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, TMATLL, lmmaxd, atom_indices, kvstr)
    double complex, intent(inout) :: mat_B(:,:)
    double complex, intent(in) :: TMATLL(lmmaxd,lmmaxd,*)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: atom_indices(:)
    integer, intent(in) :: kvstr(:)

    integer :: start, ii, num_atoms, atom_index, lm2, lmmax1, lmmax2 

    mat_B = ZERO

    num_atoms = size(atom_indices)

    do ii = 1, num_atoms

      atom_index = atom_indices(ii)

      lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)

#ifndef NDEBUG
          if (lmmax1 /= lmmaxd) then
            write (*,*) "Central atom not treated with highest lmax", atom_index
            STOP
          endif
#endif

      !lmmax2 = lmmaxd
      lmmax2 = lmmax1
      ! use naive truncation: lmmax1 = lmmax2
      ! Note: this is irrelevant, since the central atom
      ! should always be treated with the highest lmax

      start = kvstr(atom_index) - 1
      do lm2 = 1, lmmax2
!         do lm1 = 1, lmmax1
                      ! TODO: WHY DO I NEED A MINUS SIGN HERE? CHECK
!           mat_B(start+lm1,(ii-1)*lmmax2+lm2) = - TMATLL(lm1,lm2,atom_index)
!           if (lm1 == lm2) mat_B(start+lm1,(ii-1)*lmmax2+lm2) = - CONE
!         enddo ! lm1
        mat_B(start+lm2,(ii-1)*lmmax2+lm2) = - CONE
      enddo ! lm2

    enddo ! ii

  endsubroutine ! buildRightHandSide

!==============================================================================
! Routines for testing
!==============================================================================

  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertToFullMatrix(smat, ia, ja, ka, kvstr, kvstc, full)
    double complex, intent(in) :: smat(:)
    integer, intent(in) :: ia(:)
    integer, intent(in) :: ja(:)
    integer, intent(in) :: ka(:)
    integer, intent(in) :: kvstr(:)
    integer, intent(in) :: kvstc(:)
    double complex, intent(out) :: full(:,:)

    integer :: iblockrow, iblockcol
    integer :: irow, icol, nrows, ind, ind_ia
    integer :: istartcol, istartrow, istopcol, istoprow

    full = ZERO

    nrows = size(ia) - 1

    do iblockrow = 1, nrows
      ind = ka(ia(iblockrow))
      do ind_ia = ia(iblockrow), ia(iblockrow+1) - 1

        iblockcol = ja(ind_ia)

        istartcol = kvstc(iblockcol)
        istopcol  = kvstc(iblockcol+1) - 1
        istartrow = kvstr(iblockrow)
        istoprow  = kvstr(iblockrow+1) - 1

        do icol = istartcol, istopcol
          do irow = istartrow, istoprow
            full(irow,icol) = smat(ind)
            ind = ind + 1
          enddo ! irow
        enddo ! icol

      enddo ! ind_ia
    enddo ! iblockrow

  endsubroutine ! convertToFullMatrix

  !----------------------------------------------------------------------------
  !> Solution of a system of linear equations with multiple right hand sides,
  !> using standard dense matrix LAPACK routines.
  subroutine solveFull(full, mat_B, mat_X)
    double complex, intent(inout) :: full(:,:)
    double complex, intent(in)  :: mat_B(:,:)
    double complex, intent(out) :: mat_X(:,:)

    integer, allocatable :: ipvt(:)
    integer :: ndim, num_rhs, info

    ndim = size(full, 1)
    allocate(ipvt(ndim))

    num_rhs = size(mat_B, 2)
    mat_X(:,:) = mat_B(:,:)
    
    CALL ZGETRF(ndim,ndim,full,ndim,IPVT,info)
    CALL ZGETRS('N',ndim,num_rhs,full,ndim,IPVT,mat_X,ndim,info)

    deallocate(ipvt, stat=info)
  endsubroutine ! solveFull


  !------------------------------------------------------------------------------
  !> Convert solution with l-cutoff to format of solution without l-cutoff.
  !>
  !> @param GLLKE1 output: solution in old format
  !> @param mat_X: solution with l-cutoff
  subroutine toOldSolutionFormat(gllke1, mat_x, lmmaxd, kvstr)
    double complex, intent(out) :: gllke1(:,:)
    double complex, intent(in) :: mat_x(:,:)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: kvstr(:)

    integer :: num_atoms, atom_index, start, lm1, lmmax1

    gllke1 = ZERO
    num_atoms = size(kvstr) - 1

    do atom_index = 1, num_atoms

      start = kvstr(atom_index) - 1
      lmmax1 = kvstr(atom_index + 1) - kvstr(atom_index)

      do lm1 = 1, lmmax1
        gllke1((atom_index-1)*lmmaxd+lm1,:) = mat_x(start+lm1,:)
      enddo ! lm1

    enddo ! atom_index

  endsubroutine ! toOldSolutionFormat

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to unformatted file
  !> - useful for testing.
  subroutine dumpSparseMatrixData(smat, filename)
    double complex, intent(in) :: smat(:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) smat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read sparse matrix data from unformatted file
  !> - useful for testing.
  subroutine readSparseMatrixData(smat, filename)
    double complex, intent(out) :: smat(:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) smat
    close(fu)
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Write dense matrix to unformatted file
  !> - useful for testing.
  subroutine dumpDenseMatrix(mat, filename)
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
  subroutine readDenseMatrix(mat, filename)
    double complex, intent(out) :: mat(:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) mat
    close(fu)
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Write GLLH to unformatted file
  !> - useful for testing.
  subroutine dumpGLLH(mat, filename)
    double complex, intent(in) :: mat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) mat
    close(fu)
  endsubroutine ! dump

  !----------------------------------------------------------------------------
  !> Read GLLH from unformatted file
  !> - useful for testing.
  subroutine readGLLH(mat, filename)
    double complex, intent(out) :: mat(:,:,:)
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 97

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) mat
    close(fu)
  endsubroutine ! read

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

endmodule
