module fillKKRMatrix_mod

contains

  subroutine getKKRMatrixStructure(lmmaxd_array, numn0, indn0, & ! in
                                   ia, ja, kvstr) ! out
    implicit none
    integer, dimension(:), intent(in) :: lmmaxd_array
    integer, dimension(:), intent(in) :: numn0
    integer, dimension(:,:), intent(in) :: indn0
    integer, dimension(:), intent(out) :: ia
    integer, dimension(:), intent(out) :: ja
    integer, dimension(:), intent(out) :: kvstr

    !----- local
    integer :: nnz_blocks
    integer :: nrows
    integer :: ii
    integer :: irow, icol


    nrows = size(lmmaxd_array)

    !ASSERT(size(numn0) == nrows)
    !ASSERT(size(indn0, 1) == nrows)
    !ASSERT(size(ia) == nrows + 1)
    !ASSERT(size(kvstr == nrows + 1)

    nnz_blocks = 0
    do ii = 1, nrows
      nnz_blocks = nnz_blocks + numn0(ii)
    end do

    !ASSERT(nnz_blocks < size(ja))

    kvstr(1) = 1
    do ii = 2, nrows + 1
      kvstr(ii) = kvstr(ii-1) + lmmaxd_array(ii - 1)
    end do

    ii = 1
    do irow = 1, nrows
      do icol = 1, numn0(irow)  ! square matrix
        !ASSERT(icol <= size(indn0, 2))
        !ASSERT(ii <= nnz_blocks)
        ja(ii) = indn0(irow, icol)
        ii = ii + 1
      end do
      ia(irow) = ii
    end do
    ia(nrows + 1) = ii

  ! DONE!
  end subroutine

end module
