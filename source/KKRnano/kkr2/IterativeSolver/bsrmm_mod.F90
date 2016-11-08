!> Multiplication of two BSR-matrices (block sparse row matrices)

module bsrmm_mod
  implicit none
  private
  public :: bsr_times_bsr
  
  contains

  subroutine bsr_times_bsr(Y, ia, ja, A, ix, jx, X, nFlops)
    double complex, intent(out) :: Y(:,:,:) ! dim(lmsa,lmsd,X%nnzb)
    integer, intent(in) :: ia(:) !> dim(A%nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)      ! column indices
    double complex, intent(in)  :: A(:,:,:) ! dim(lmsa,lmsd,A%nnzb)
    integer, intent(in) :: ix(:) !> dim(X%nRows + 1) !  start indices
    integer, intent(in) :: jx(:) !> dim(X%nnzb)      ! column indices
    double complex, intent(in)  :: X(:,:,:) ! dim(lmsa,lmsd,X%nnzb)
    integer(kind=8), intent(inout) :: nFlops
    ! for the result Y we assume the same BSR structure as X
#define iy ix
#define jy jx

    external :: ZGEMM ! from BLAS

    ! local variables
    integer :: XnRows, leadDim_Y, leadDim_X, leadDim_A, lmsd
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    ! private variables
    integer :: yRow, Yind, yCol, Aind, aCol, Xind
    double complex :: beta

    leadDim_A = size(A, 1)
    leadDim_X = size(X, 1)
    leadDim_Y = size(Y, 1)
    
    lmsd = size(A, 2)

    if (leadDim_Y /= leadDim_X) stop __LINE__
    if (leadDim_A /= leadDim_X) stop __LINE__
    if (size(X, 3) /= size(Y, 3)) stop __LINE__

    XnRows = size(ix) - 1
#define YnRows XnRows
    if (XnRows /= size(ia) - 1) stop __LINE__

    do yRow = 1, YnRows
      ! reuse elements of Y, i.e. keep the accumulator in the cache
      do Yind = iy(yRow), iy(yRow + 1) - 1 
        yCol = jy(Yind) ! update   matrix block element Y_full[yRow,yCol]

        beta = zero
        
        do Aind = ia(yRow), ia(yRow + 1) - 1
          aCol = ja(Aind) ! for each matrix block element A_full[yRow,aCol]

          Xind = BSR_entry_exists(ix, jx, row=aCol, col=yCol) ! find out, if X_full[aCol,yCol] exists
          if (Xind > 0) then ! yes

            ! now: Y[:,:,Yind] += A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += A(m,k)*B(k,n)
            !                    M     N     K          A                       B                             C
            call zgemm('n', 'n', lmsd, lmsd, lmsd, one, A(:,1,Aind), leadDim_A, X(:,1,Xind), leadDim_X, beta, Y(:,1,Yind), leadDim_Y)

            nFlops = nFlops + (8_8*lmsd)*(lmsd*lmsd)

            beta = one
          endif ! X_full[aCol,yCol] exists

        enddo ! Aind
      enddo ! Yind
    enddo ! yRow

  endsubroutine ! bsr_times_bsr

  integer function BSR_entry_exists(RowStart, Colindex, row, col) result(Ind)
    integer, intent(in) :: RowStart(:), ColIndex(:), row, col

    ! ToDo: the ColIndex list should be sorted ascendingly, so bisection search will be faster
    do Ind = RowStart(row), RowStart(row + 1) - 1
      if (ColIndex(Ind) == col) return ! Ind
    enddo ! Ind
    Ind = -1 ! not found
    
  endfunction ! exists

endmodule ! vbrmv_mat_mod
