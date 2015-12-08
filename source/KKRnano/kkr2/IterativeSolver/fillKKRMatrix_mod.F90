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
  public :: dump

  double complex, parameter :: ZERO=(0.d0, 0.d0), CONE=(1.d0, 0.d0)
  
  interface dump
    module procedure dumpDenseMatrix, dumpSparseMatrixData
  endinterface
  
  contains

  !------------------------------------------------------------------------------
  !> Setup of the sparsity pattern of the KKR-Matrix.
  subroutine getKKRMatrixStructure(lmmaxd_array, numn0, indn0, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    integer, intent(in) :: lmmaxd_array(:) !< block size of each row, dim(nrows)
    integer, intent(in) :: numn0(:) !< dim(nrows)
    integer, intent(in) :: indn0(:,:) !< dim(nrows,maxval(numn0))
    type(SparseMatrixDescription), intent(inout) :: sparse

    integer :: nnzb, nrows, ij, irow, icol, start_address

    nrows = size(lmmaxd_array)

    ASSERT(size(numn0) == nrows)
    ASSERT(size(indn0, 1) == nrows)
    ASSERT(size(sparse%ia) == nrows+1)
    ASSERT(size(sparse%kvstr) == nrows+1)

    sparse%ia = 0
    sparse%ja = 0
    sparse%ka = 0
    sparse%kvstr = 0 ! init block start indices

    nnzb = sum(numn0(1:nrows)) ! number of non-zero blocks

    ASSERT(nnzb <= size(sparse%ja))
    ASSERT(nnzb + 1 <= size(sparse%ka))

    sparse%kvstr(1) = 1
    do ij = 1, nrows
      sparse%kvstr(ij+1) = sparse%kvstr(ij) + lmmaxd_array(ij)
    enddo ! ij

    ij = 1
    do irow = 1, nrows
      sparse%ia(irow) = ij
      do icol = 1, numn0(irow) ! square matrix
      
        ASSERT(icol <= size(indn0, 2))
        ASSERT(ij <= nnzb)
        
        sparse%ja(ij) = indn0(irow,icol)
        ij = ij + 1
      enddo ! icol
    enddo ! irow
    sparse%ia(nrows+1) = ij ! final, important since the ranges are always [ia(i) ... ia(i+1)-1]

    start_address = 1
    ij = 1
    do irow = 1, nrows
      do icol = 1, numn0(irow)
        sparse%ka(ij) = start_address

        ASSERT( 1 <= indn0(irow,icol) .and. indn0(irow,icol) <= nrows )

        start_address = start_address + lmmaxd_array(irow)*lmmaxd_array(indn0(irow,icol))
        ij = ij + 1
      enddo ! icol
    enddo ! irow
    sparse%ka(ij) = start_address

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
  subroutine buildKKRCoeffMatrix(smat, tmatLL, lmmaxd, num_atoms, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    double complex, intent(inout) :: smat(:)
    double complex, intent(in) :: tmatLL(lmmaxd,lmmaxd,num_atoms)
    integer, intent(in) :: lmmaxd, num_atoms
    type(SparseMatrixDescription), intent(in) :: sparse

    double complex :: temp(lmmaxd)
    integer :: block_row, block_col, start, lm1, lm2, lm3, lmmax1, lmmax2, lmmax3
    integer :: istart_row, istop_row, istart_col, istop_col, ind_ia

    start = 0
    do block_row = 1, num_atoms
      istart_row = sparse%kvstr(block_row)
      istop_row  = sparse%kvstr(block_row+1)
      do ind_ia = sparse%ia(block_row), sparse%ia(block_row+1)-1

        block_col = sparse%ja(ind_ia) ! ja gives the block-column indices of non-zero blocks

        istart_col = sparse%kvstr(block_col)
        istop_col  = sparse%kvstr(block_col+1)

#ifndef NDEBUG
        if (1 > block_row .or. block_row > num_atoms) then
          write (*,*) "buildKKRCoeffMatrix: invalid block_row", block_row
          stop
        endif

        if (1 > block_col .or. block_col > num_atoms) then
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
        ! -|    |  lmmax1     *    | G_ref    | lmmax3  =   |  -T*G    |  lmmax1
        !  |----|                  |----------| (==lmmax1)  |----------|
        !  lmmax3==lmmax1

        do lm2 = 1, lmmax2

          temp(:) = ZERO
          do lm3 = 1, lmmax3
            temp(:) = temp(:) - tmatLL(:,lm3,block_row) * smat(start+lm3+lmmax3*(lm2-1)) ! -T*G
          enddo ! lm3
          
          if (block_row == block_col) temp(lm2) = temp(lm2) + CONE ! add 1.0 on the diagonal

          do lm1 = 1, lmmax1

            smat(start+lm1+lmmax1*(lm2-1)) = temp(lm1)
            
          enddo ! lm1
        enddo ! lm2

        start = start + lmmax2*lmmax1
      enddo ! ind_ia ! block columns
    enddo ! block rows

  endsubroutine ! buildKKRCoeffMatrix

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, lmmaxd, atom_indices, kvstr, tmatLL)
    double complex, intent(inout) :: mat_B(:,:)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: atom_indices(:)
    integer, intent(in) :: kvstr(:)
    double complex, intent(in), optional :: tmatLL(lmmaxd,lmmaxd,*)

    integer :: start, ii, num_atoms, atom_index, lm1, lm2, lmmax1, lmmax2 

    mat_B = ZERO

    num_atoms = size(atom_indices)

    do ii = 1, num_atoms

      atom_index = atom_indices(ii)

      lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)

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

      start = kvstr(atom_index) - 1
      do lm2 = 1, lmmax2
        if (present(tmatLL)) then
          do lm1 = 1, lmmax1
            ! TODO: WHY DO I NEED A MINUS SIGN HERE? CHECK
            mat_B(start+lm1,lm2+lmmax2*(ii-1)) = -tmatLL(lm1,lm2,atom_index)
          enddo ! lm1
        else
          mat_B(start+lm2,lm2+lmmax2*(ii-1)) = -CONE
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
  subroutine convertToFullMatrix(smat, ia, ja, ka, kvstr, kvstc, full)
    double complex, intent(in) :: smat(:)
    integer, intent(in) :: ia(:)
    integer, intent(in) :: ja(:)
    integer, intent(in) :: ka(:)
    integer, intent(in) :: kvstr(:)
    integer, intent(in) :: kvstc(:)
    double complex, intent(out) :: full(:,:)

    integer :: ibrow, ibcol, irow, icol, nrows, ind, ind_ia

    full = ZERO

    nrows = size(ia) - 1

    do ibrow = 1, nrows
      ind = ka(ia(ibrow))
      do ind_ia = ia(ibrow), ia(ibrow+1) - 1

        ibcol = ja(ind_ia)

        do icol = kvstc(ibcol), kvstc(ibcol+1) - 1
          do irow = kvstr(ibrow), kvstr(ibrow+1) - 1
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
    use TruncationZone_mod, only: clear_non_existing_entries
    double complex, intent(inout) :: full(:,:)
    double complex, intent(in)  :: mat_B(:,:)
    double complex, intent(out) :: mat_X(:,:)

    integer, allocatable :: ipvt(:)
    integer :: ndim, num_rhs, info

    ndim = size(full, 1)
    allocate(ipvt(ndim))

    num_rhs = size(mat_B, 2)
    mat_X(:,:) = mat_B(:,:)
    
    call zgetrf(ndim,ndim,full,ndim,ipvt,info)
    call zgetrs('n',ndim,num_rhs,full,ndim,ipvt,mat_x,ndim,info)

    deallocate(ipvt, stat=info)
    
    call clear_non_existing_entries(mat_X) ! make up for treating more than one atom with truncation
  endsubroutine ! solveFull


  !------------------------------------------------------------------------------
  !> Convert solution with l-cutoff to format of solution without l-cutoff.
  !>
  !> @param GLLKE1 output: solution in old format
  !> @param mat_X: solution with l-cutoff
  subroutine toOldSolutionFormat(gllke1, mat_X, lmmaxd, kvstr)
    double complex, intent(out) :: gllke1(:,:)
    double complex, intent(in) :: mat_X(:,:)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: kvstr(:)

    integer :: num_atoms, atom_index, start, lm1, lmmax1

    gllke1 = ZERO
    num_atoms = size(kvstr) - 1

    do atom_index = 1, num_atoms

      start = kvstr(atom_index) - 1
      lmmax1 = kvstr(atom_index + 1) - kvstr(atom_index)

      do lm1 = 1, lmmax1
        gllke1((atom_index-1)*lmmaxd+lm1,:) = mat_X(start+lm1,:)
      enddo ! lm1

    enddo ! atom_index

  endsubroutine ! toOldSolutionFormat
  
  
  
  
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

endmodule
