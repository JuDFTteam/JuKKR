!> Multiplication of VBR-matrix (variable block row sparse matrix) with dense
!> matrix.

!     @PROCESS HOT=noarraypad:level=1:simd:vector
!     See:
!     Y. Saad, SPARSKIT: a basic tool kit for sparse matrix computations - Version 2 (1994)
!     modified for double complex and for
!     multiplying several vectors at once
!     => matrix-matrix multiplication instead of matrix-vector multiplication


module vbrmv_mat_mod
  implicit none
  private
  public :: multiply_vbr
  
  contains

  subroutine multiply_vbr(a, x, b, sparse)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription

    type (SparseMatrixDescription), intent(in) :: sparse
    double complex a(:), x(:,:), b(:,:)

    call vbrmv_mat(sparse%blk_nrows, sparse%ia, sparse%ja, sparse%ka, &
                   a, sparse%kvstr, sparse%kvstr, x, b, &
                   sparse%max_blockdim, sparse%max_blocks_per_row)

  end subroutine


!     Heavily modified routine from SPARSKIT
!-----------------------------------------------------------------------
  subroutine vbrmv_mat(blk_nrows, ia, ja, ka, a, kvstr, kvstc, x, b, &
                       max_blockdim, max_blocks_per_row)
                       
    integer blk_nrows, ia(blk_nrows+1), ja(:), ka(:), kvstr(blk_nrows+1), kvstc(:)
    integer max_blockdim, max_blocks_per_row
    integer ncols
    double complex  a(:), x(:,:), b(:,:)
    !-----------------------------------------------------------------------
    !     Sparse matrix-full vector product, in VBR format.
    !-----------------------------------------------------------------------
    !     On entry:
    !--------------
    !     blk_nrows      = number of block rows in matrix A
    !     ia,ja,ka,a,kvstr,kvstc = matrix A in variable block row format
    !     x       = multiplier vector in full format

    !     ncols = number of columns of matrix A

    !     On return:
    !---------------
    !     b = product of matrix A times vector x in full format

    !     Algorithm:
    !---------------
    !     Perform multiplication by traversing a in order.

    !-----------------------------------------------------------------------
    !-----local variables
    integer i, j, k, istart, istop

    integer icols

!IBM* ALIGN(32, buffer)
    double complex buffer(max_blockdim*max_blocks_per_row, size(x,2))

    integer nrowbuf, rowbuf, sum_nrowbuf
    integer startrow, leaddim_b, num_rows, leaddim_buffer
    integer jlow, jhigh

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    double complex, parameter :: CONE  = (1.0d0, 0.0d0)

    leaddim_b = size(b, 1)
    leaddim_buffer = size(buffer, 1)

    ncols = size(b, 2)

    b = 0.d0

!     can parallelise this loop

!$OMP PARALLEL PRIVATE(i, istart, istop,num_rows,rowbuf, &
!$OMP                  sum_nrowbuf,j,jlow,jhigh, &
!$OMP                  startrow,nrowbuf,icols,buffer,k)

!$OMP DO
    do i = 1, blk_nrows
      istart = kvstr(i)
      istop  = kvstr(i+1)-1
      num_rows = istop - istart + 1
      rowbuf = 1
      sum_nrowbuf = 0

      jlow = ia(i)
      jhigh = ia(i+1)-1

      k = ka(jlow)

      do j = jlow, jhigh
        startrow = kvstc(ja(j))
        nrowbuf = kvstc(ja(j)+1) - startrow

!IBM* ASSERT(ITERCNT(16))
        do icols = 1, ncols
          !call ZCOPY(nrowbuf, x(startrow, icols), 1, &
          !               buffer(rowbuf,   icols), 1)
          buffer(rowbuf:(rowbuf + nrowbuf -1), icols) = &
          x(startrow:(startrow + nrowbuf - 1), icols)
        enddo

        sum_nrowbuf = sum_nrowbuf + nrowbuf
        rowbuf = rowbuf + nrowbuf
      enddo

      !k = ka(ia(i))

      call ZGEMM('N','N',num_rows,ncols,sum_nrowbuf, &
      CONE,a(k),num_rows, &
      buffer,leaddim_buffer, &
      CZERO,b(istart,1),leaddim_b)

    enddo
!$OMP END DO

!$OMP END PARALLEL

  end subroutine

end module
