!> Multiplication of VBR-matrix (variable block row sparse matrix) with dense
!> matrix.

!     @PROCESS HOT=noarraypad:level=1:simd:vector
!     See:
!     Y. Saad, SPARSKIT: a basic tool kit for sparse matrix computations - Version 2 (1994)
!     modified for double complex and for
!     multiplying several vectors at once
!     => matrix-matrix multiplication instead of matrix-vector multiplication


module vbrmv_mat_mod
  use TruncationZone_mod, only: clear_non_existing_entries
#define CLEAR_NON_EXISTING_ENTRIES(B, N) call clear_non_existing_entries(B, N)
  implicit none
  private
  public :: multiply_vbr
  
  contains

  subroutine multiply_vbr(a, x, b, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    double complex, intent(in)  :: a(:), x(:,:)
    double complex, intent(out) :: b(:,:)
    type(SparseMatrixDescription), intent(in) :: sparse

    call vbrmv_mat(sparse%blk_nrows, sparse%ia, sparse%ja, sparse%ka, &
                   a, sparse%kvstr, sparse%kvstr, x, b, &
                   sparse%max_blockdim, sparse%max_blocks_per_row)

  endsubroutine ! multiply_vbr

  !> Heavily modified routine from SPARSKIT
  subroutine vbrmv_mat(blk_nrows, ia, ja, ka, A, kvstr, kvstc, x, b, max_blockdim, max_blocks_per_row)
                       
    integer, intent(in) :: blk_nrows, ia(blk_nrows+1), ja(:), ka(:), kvstr(blk_nrows+1), kvstc(:)
    integer, intent(in) :: max_blockdim, max_blocks_per_row
    double complex, intent(in)  :: A(:), x(:,:)
    double complex, intent(out) :: b(:,:)
    !-----------------------------------------------------------------------
    !     Sparse matrix-full vector product, in VBR format.
    !-----------------------------------------------------------------------
    !     On entry:
    !--------------
    !     blk_nrows      = number of block rows in matrix A
    !     ia,ja,ka,A,kvstr,kvstc = matrix A in variable block row format
    !     x       = multiplier vector in full format

    !     ncols = number of columns of matrix A

    !     On return:
    !---------------
    !     b = product of matrix A times vector x in full format

    !     Algorithm:
    !---------------
    !     Perform multiplication by traversing A in order.

    !-----------------------------------------------------------------------
    !-----local variables

!IBM* ALIGN(32, buffer)
    double complex :: buffer(max_blockdim*max_blocks_per_row,size(x, 2))

    integer :: ibr, ibc, j, k, istart, ncols, icols, nrowbuf, sum_nrowbuf
    integer :: startrow, leaddim_b, nrows, leaddim_buffer

    double complex, parameter :: CZERO = (0.d0, 0.d0), CONE  = (1.d0, 0.d0)

    leaddim_b = size(b, 1)
    ncols = size(b, 2)
    leaddim_buffer = max_blockdim*max_blocks_per_row

    b = CZERO

!     can parallelise this loop

!$OMP PARALLEL PRIVATE(ibr,istart,nrows,sum_nrowbuf,j,startrow,nrowbuf,icols,buffer,k)
!$OMP DO
    do ibr = 1, blk_nrows
      istart = kvstr(ibr)
      nrows  = kvstr(ibr+1) - istart
      sum_nrowbuf = 0

      k = ka(ia(ibr))
      do j = ia(ibr), ia(ibr+1)-1
        ibc = ja(j)
        startrow = kvstc(ibc)
        nrowbuf  = kvstc(ibc+1) - startrow

!IBM* ASSERT(ITERCNT(16))
        do icols = 1, ncols
          !call ZCOPY(nrowbuf, x(startrow, icols), 1, buffer(sum_nrowbuf+1,icols), 1)
          buffer(sum_nrowbuf+1:sum_nrowbuf+nrowbuf,icols) = x(startrow:startrow+nrowbuf-1,icols)
        enddo ! icols

        sum_nrowbuf = sum_nrowbuf + nrowbuf
      enddo ! j

      call ZGEMM('N','N',nrows,ncols,sum_nrowbuf,CONE,A(k:),nrows,buffer,leaddim_buffer,CZERO,b(istart,1),leaddim_b)

    enddo ! ibr
!$OMP endDO
!$OMP endPARALLEL

    CLEAR_NON_EXISTING_ENTRIES(b, max_blockdim)
    
  endsubroutine ! vbrmv_mat

endmodule ! vbrmv_mat_mod
