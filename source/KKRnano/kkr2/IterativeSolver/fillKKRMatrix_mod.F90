#ifndef NDEBUG
#define ASSERT(CONDITION) if (.not. (CONDITION)) then; write(*,*) "Assertion ", #CONDITION, " failed: ",  __FILE__, __LINE__; endif
#else
#define ASSERT(CONDITION)
#endif


module fillKKRMatrix_mod

contains

!------------------------------------------------------------------------------
!> Setup of the sparsity pattern of the KKR-Matrix.
  subroutine getKKRMatrixStructure(lmmaxd_array, numn0, indn0, & ! in
                                   sparse) ! out

    use SparseMatrixDescription_mod
    implicit none
    integer, dimension(:), intent(in) :: lmmaxd_array
    integer, dimension(:), intent(in) :: numn0
    integer, dimension(:,:), intent(in) :: indn0
    type (SparseMatrixDescription), intent(inout) :: sparse

    !----- local
    integer :: nnz_blocks
    integer :: nrows
    integer :: ii
    integer :: irow, icol
    integer :: start_address

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
    end do

    ASSERT(nnz_blocks <= size(sparse%ja))
    ASSERT(nnz_blocks + 1 <= size(sparse%ka))

    sparse%kvstr(1) = 1
    do ii = 2, nrows + 1
      sparse%kvstr(ii) = sparse%kvstr(ii-1) + lmmaxd_array(ii - 1)
    end do

    ii = 1
    do irow = 1, nrows
      sparse%ia(irow) = ii
      do icol = 1, numn0(irow)  ! square matrix
        ASSERT(icol <= size(indn0, 2))
        ASSERT(ii <= nnz_blocks)
        sparse%ja(ii) = indn0(irow, icol)
        ii = ii + 1
      end do
    end do
    sparse%ia(nrows + 1) = ii

    start_address = 1
    ii = 1
    do irow = 1, nrows
      do icol = 1, numn0(irow)
      sparse%ka(ii) = start_address

      ASSERT( indn0(irow, icol) >= 1 .and. indn0(irow, icol) <= nrows )

      start_address = start_address + lmmaxd_array(irow)*lmmaxd_array(indn0(irow,icol))
      ii = ii + 1
      end do
    end do
    sparse%ka(ii) = start_address

    sparse%max_blockdim = maxval(lmmaxd_array)
    sparse%max_blocks_per_row = maxval(numn0)

  ! DONE!
  end subroutine

  !---------------------------------------------------------------------
  !> Given G_ref - build (1 - TG_ref) -> coefficent matrix.
  !>
  !> @param smat  on input:  sparse block matrix containing Gref
  !>              on output: coefficient matrix (1 - TG_ref)
  !> @param sparse  sparse matrix description
  !> ia    for each row give index of first non-zero block in ja
  !> ja    column index array of non-zero blocks
  subroutine buildKKRCoeffMatrix(smat, TMATLL, lmmaxd, num_atoms, sparse)
    use SparseMatrixDescription_mod
    implicit none
    double complex, dimension(:), intent(inout) :: smat
    double complex, dimension(lmmaxd,lmmaxd,num_atoms), intent(in) :: TMATLL
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: num_atoms

    type (SparseMatrixDescription), intent(inout) :: sparse
    ! ------- local

    double complex, dimension(lmmaxd) :: temp
    integer :: block_row
    integer :: block_col
    integer :: start
    integer :: lm1
    integer :: lm2
    integer :: lm3
    integer :: lmmax1, lmmax2, lmmax3
    integer :: istart_row,istop_row
    integer :: istart_col,istop_col, ind_ia

    double complex, parameter :: CZERO =(0.0D0,0.0D0)
    double complex :: CONE = (1.0D0,0.0D0)

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
        end if

        if (block_col < 1 .or. block_col > num_atoms) then
          write (*,*) "buildKKRCoeffMatrix: invalid block_col", block_col
          STOP
        end if
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

          temp = CZERO
          do lm1 = 1, lmmax1
            do lm3 = 1, lmmax3
              temp(lm1) = temp(lm1) - TMATLL(lm1,lm3,block_row) * smat(start + (lm2 - 1) * lmmax3 + lm3)  ! -T*G
            end do
          end do

          do lm1 = 1, lmmax1

            if (block_row == block_col .and. lm1 == lm2) then
              temp(lm1) = temp(lm1) + CONE  ! add 1 on the diagonal
            end if

            smat(start + (lm2 - 1) * lmmax1 + lm1) = temp(lm1)
            !smat(start + (lm2 - 1) * lmmax1 + lm1) = (block_row * 100.d0 + block_col) * CONE
          end do

        end do ! lm2

        start = start + lmmax2*lmmax1
      end do ! block columns
    end do ! block rows

  end subroutine

!------------------------------------------------------------------------------
!> Builds the right hand site for the linear KKR matrix equation.
  subroutine buildRightHandSide(mat_B, TMATLL, lmmaxd, atom_index, kvstr)
    implicit none

    double complex, dimension(:,:), intent(inout) :: mat_B
    double complex, dimension(lmmaxd,lmmaxd,*), intent(in) :: TMATLL
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: atom_index
    integer, dimension(:), intent(in) :: kvstr

    double complex, parameter :: CZERO =(0.0D0,0.0D0)
    integer :: start
    integer :: lm1, lm2
    integer :: lmmax1, lmmax2

    mat_B = CZERO

    lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)

#ifndef NDEBUG
        if (lmmax1 /= lmmaxd) then
          write (*,*) "Central atom not treated with highest lmax", atom_index
          STOP
        end if
#endif

    !lmmax2 = lmmaxd
    lmmax2 = lmmax1
    ! use naive truncation: lmmax1 = lmmax2
    ! Note: this is irrelevant, since the central atom
    ! should always be treated with the highest lmax

    start = kvstr(atom_index) - 1
    do lm2 = 1, lmmax2
      do lm1 = 1, lmmax1
                    ! TODO: WHY DO I NEED A MINUS SIGN HERE? CHECK
        mat_B( start + lm1, lm2 ) = - TMATLL(lm1, lm2, atom_index)
      end do
    end do

  end subroutine

!==============================================================================
! Routines for testing
!==============================================================================

  !----------------------------------------------------------------------------
  !> Given the sparse matrix data 'smat' and the sparsity information,
  !> create the dense matrix representation of the matrix.
  subroutine convertToFullMatrix(smat, ia, ja, ka, kvstr, kvstc, full)
    implicit none
    double complex, dimension(:), intent(in) :: smat
    integer, dimension(:), intent(in) :: ia
    integer, dimension(:), intent(in) :: ja
    integer, dimension(:), intent(in) :: ka
    integer, dimension(:), intent(in) :: kvstr
    integer, dimension(:), intent(in) :: kvstc
    double complex, dimension(:,:), intent(out) :: full

    !------------
    integer :: iblockrow
    integer :: iblockcol
    integer :: irow
    integer :: icol
    integer :: nrows
    integer :: ind
    integer :: ind_ia
    integer :: istartcol, istartrow, istopcol, istoprow

    double complex, parameter :: CZERO =(0.0D0,0.0D0)

    full = CZERO

    nrows = size(ia) - 1

    do iblockrow = 1, nrows
      ind = ka(ia(iblockrow))
      do ind_ia = ia(iblockrow), ia(iblockrow + 1) - 1

        iblockcol = ja(ind_ia)

        istartcol = kvstc(iblockcol)
        istopcol  = kvstc(iblockcol + 1) - 1
        istartrow = kvstr(iblockrow)
        istoprow  = kvstr(iblockrow + 1) - 1

        do icol = istartcol, istopcol
          do irow = istartrow, istoprow
            full(irow, icol) = smat(ind)
            ind = ind + 1
          end do
        end do

      end do
    end do

  end subroutine

  !----------------------------------------------------------------------------
  !> Solution of a system of linear equations with multiple right hand sides,
  !> using standard dense matrix LAPACK routines.
  subroutine solveFull(full, mat_B)
    implicit none

    double complex, dimension(:,:), intent(inout) :: full
    double complex, dimension(:,:), intent(inout) :: mat_B

    !------------
    integer, dimension(:), allocatable :: ipvt
    integer :: info
    integer :: ndim
    integer :: num_rhs

    ndim = size(full, 1)

    allocate(ipvt(ndim))

    num_rhs = size(mat_B, 2)

    CALL ZGETRF(ndim,ndim,full,ndim,IPVT,INFO)
    CALL ZGETRS('N',ndim,num_rhs,full,ndim,IPVT,mat_B,ndim,INFO)

    deallocate(ipvt)

  end subroutine


  !------------------------------------------------------------------------------
  !> Convert solution with l-cutoff to format of solution without l-cutoff.
  !>
  !> @param GLLKE1 output: solution in old format
  !> @param mat_X: solution with l-cutoff
  subroutine toOldSolutionFormat(GLLKE1, mat_X, lmmaxd, kvstr)
    implicit none

    double complex, dimension(:,:), intent(out) :: GLLKE1
    double complex, dimension(:,:), intent(in) :: mat_X
    integer, intent(in) :: lmmaxd
    integer, dimension(:), intent(in) :: kvstr
    !---------------------

    integer :: num_atoms, atom_index
    double complex, parameter :: CZERO =(0.0D0,0.0D0)
    integer :: start
    integer :: lm1
    integer :: lmmax1

    GLLKE1 = CZERO
    num_atoms = size(kvstr) - 1


    do atom_index = 1, num_atoms

      start = kvstr(atom_index) - 1
      lmmax1 = kvstr(atom_index + 1) - kvstr(atom_index)

      do lm1 = 1, lmmax1
        GLLKE1((atom_index - 1)*lmmaxd + lm1, : ) = mat_X(start + lm1, :)
      end do

    end do

  end subroutine

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to unformatted file
  !> - useful for testing.
  subroutine dumpSparseMatrixData(smat, filename)
    implicit none
    double complex, dimension(:), intent(in) :: smat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='unformatted')
    write(FILEHANDLE) smat
    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Read sparse matrix data from unformatted file
  !> - useful for testing.
  subroutine readSparseMatrixData(smat, filename)
    implicit none
    double complex, dimension(:), intent(inout) :: smat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='unformatted')
    read(FILEHANDLE) smat
    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write dense matrix to unformatted file
  !> - useful for testing.
  subroutine dumpDenseMatrix(mat, filename)
    implicit none
    double complex, dimension(:,:), intent(in) :: mat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='unformatted')
    write(FILEHANDLE) mat
    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Read dense matrix data from unformatted file
  !> - useful for testing.
  subroutine readDenseMatrix(mat, filename)
    implicit none
    double complex, dimension(:,:), intent(inout) :: mat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='unformatted')
    read(FILEHANDLE) mat
    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write GLLH to unformatted file
  !> - useful for testing.
  subroutine dumpGLLH(mat, filename)
    implicit none
    double complex, dimension(:,:,:), intent(in) :: mat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='unformatted')
    write(FILEHANDLE) mat
    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Read GLLH from unformatted file
  !> - useful for testing.
  subroutine readGLLH(mat, filename)
    implicit none
    double complex, dimension(:,:,:), intent(inout) :: mat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97

    open(FILEHANDLE, file=filename, form='unformatted')
    read(FILEHANDLE) mat
    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write sparse matrix data (without description) to formatted file
  !> - useful for testing.
  subroutine dumpSparseMatrixDataFormatted(smat, filename)
    implicit none
    double complex, dimension(:), intent(in) :: smat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97
    integer :: ii

    open(FILEHANDLE, file=filename, form='formatted')

    do ii = 1, size(smat)
      write(FILEHANDLE, *) real(smat(ii)), imag(smat(ii))
    end do

    close(FILEHANDLE)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write dense matrix to formatted file
  !> First line gives matrix dimension: rows cols
  !> - useful for testing.
  subroutine dumpDenseMatrixFormatted(mat, filename)
    implicit none
    double complex, dimension(:,:), intent(in) :: mat
    character(len = *), intent(in) :: filename
    !--------------

    integer, parameter :: FILEHANDLE = 97
    integer :: ii, jj

    open(FILEHANDLE, file=filename, form='formatted')
    write(FILEHANDLE, *) size(mat,1), size(mat,2)
    do jj = 1, size(mat,2)
      do ii = 1, size(mat,1)
         write(FILEHANDLE, *) real(mat(ii,jj)), imag(mat(ii,jj))
      end do
    end do
    close(FILEHANDLE)
  end subroutine

end module
