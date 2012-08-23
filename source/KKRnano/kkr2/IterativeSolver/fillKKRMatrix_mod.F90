#ifndef NDEBUG
#define ASSERT(CONDITION) if (.not. (CONDITION)) then; write(*,*) "Assertion ", #CONDITION, " failed: ",  __FILE__, __LINE__; endif
#else
#define ASSERT(CONDITION)
#endif


module fillKKRMatrix_mod

!  !> Describes a square, variable block row (VBR) sparse matrix
!  !> (double complex)
!  type VBRSparseMatrixDescZ
!    integer, dimension(:), intent(out) :: ia
!    integer, dimension(:), intent(out) :: ja
!    integer, dimension(:), intent(out) :: kvstr
!  end type

contains

  subroutine getKKRMatrixStructure(lmmaxd_array, numn0, indn0, & ! in
                                   ia, ja, ka, kvstr) ! out
    implicit none
    integer, dimension(:), intent(in) :: lmmaxd_array
    integer, dimension(:), intent(in) :: numn0
    integer, dimension(:,:), intent(in) :: indn0
    integer, dimension(:), intent(out) :: ia
    integer, dimension(:), intent(out) :: ja
    integer, dimension(:), intent(out) :: ka
    integer, dimension(:), intent(out) :: kvstr

    !----- local
    integer :: nnz_blocks
    integer :: nrows
    integer :: ii
    integer :: irow, icol
    integer :: start_address


    nrows = size(lmmaxd_array)

    ASSERT(size(numn0) == nrows)
    ASSERT(size(indn0, 1) == nrows)
    ASSERT(size(ia) == nrows + 1)
    ASSERT(size(kvstr) == nrows + 1)

    ia = 0
    kvstr = 0
    ja = 0
    ka = 0

    nnz_blocks = 0
    do ii = 1, nrows
      nnz_blocks = nnz_blocks + numn0(ii)
    end do

    ASSERT(nnz_blocks <= size(ja))
    ASSERT(nnz_blocks + 1 <= size(ka))

    kvstr(1) = 1
    do ii = 2, nrows + 1
      kvstr(ii) = kvstr(ii-1) + lmmaxd_array(ii - 1)
    end do

    ii = 1
    do irow = 1, nrows
      ia(irow) = ii
      do icol = 1, numn0(irow)  ! square matrix
        ASSERT(icol <= size(indn0, 2))
        ASSERT(ii <= nnz_blocks)
        ja(ii) = indn0(irow, icol)
        ii = ii + 1
      end do
    end do
    ia(nrows + 1) = ii

    start_address = 1
    ii = 1
    do irow = 1, nrows
      do icol = 1, numn0(irow)  ! square matrix
      ka(ii) = start_address
      start_address = start_address + lmmaxd_array(irow)*lmmaxd_array(icol)
      ii = ii + 1
      end do
    end do
    ka(ii) = start_address

  ! DONE!
  end subroutine

  !---------------------------------------------------------------------
  !> Given G_ref - build (1 - TG_ref) -> coefficent matrix
  !> @param smat  on input:  sparse block matrix containing Gref
  !>              on output: coefficient matrix (1 - TG_ref)
  !> @param ia    for each row give index of first non-zero block in ja
  !> @param ja    column index array of non-zero blocks
  !> ASSUMING SQUARE BLOCKS!!! (FOR NOW)
  subroutine buildKKRCoeffMatrix(smat, TMATLL, lmmaxd, num_atoms, ia, ja, kvstr)
    implicit none
    double complex, dimension(:), intent(inout) :: smat
    double complex, dimension(lmmaxd,lmmaxd,num_atoms), intent(in) :: TMATLL
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: num_atoms

    integer, dimension(:), intent(in) :: ia
    integer, dimension(:), intent(in) :: ja
    integer, dimension(:), intent(in) :: kvstr
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
      istart_row = kvstr(block_row)
      istop_row  = kvstr(block_row+1)
      do ind_ia = ia(block_row), ia(block_row+1)-1

        block_col = ja(ind_ia)  !ja gives the block-column indices of non-zero blocks

        istart_col = kvstr(block_col)
        istop_col  = kvstr(block_col+1)

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
    integer :: istartcol, istartrow, istopcol, istoprow

    double complex, parameter :: CZERO =(0.0D0,0.0D0)

    full = CZERO

    nrows = size(ia) - 1

    do iblockrow = 1, nrows
      ind = ka(ia(iblockrow))
      do iblockcol = ja(ia(iblockrow)), ja(ia(iblockrow + 1) - 1)

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

end module
