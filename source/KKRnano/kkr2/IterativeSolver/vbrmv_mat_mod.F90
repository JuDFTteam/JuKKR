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
  public :: vbrmv_mat
  
  contains

  !> Heavily modified routine from SPARSKIT
  subroutine vbrmv_mat(blk_nrows, ia, ja, ka, A, kvstr, kvstc, x, Ax, max_blockdim, max_blocks_per_row)
                       
    integer, intent(in) :: blk_nrows, ia(blk_nrows+1), ja(:), ka(:), kvstr(:), kvstc(:)
    integer, intent(in) :: max_blockdim, max_blocks_per_row
    double complex, intent(in)  :: A(:), x(:,:)
    double complex, intent(out) :: Ax(:,:)
    !-----------------------------------------------------------------------
    !     Sparse matrix-full vector product, in VBR format.
    !-----------------------------------------------------------------------
    !     On entry:
    !--------------
    !     blk_nrows      = number of block rows in matrix A
    !     ia,ja,ka,A,kvstr,kvstc = matrix A in variable block row format
    !     x       = multiplier vector in full format

    !     nRHSs = number of columns of matrix A

    !     On return:
    !---------------
    !     Ax = product of matrix A times vector x in full format

    !     Algorithm:
    !---------------
    !     Perform multiplication by traversing A in order.

    !-----------------------------------------------------------------------
    !-----local variables

!IBM* ALIGN(32, Buffer)
    double complex :: Buffer(max_blockdim*max_blocks_per_row,size(x, 2))

    integer :: ibr, isr, j, k, nRHSs, iRHSs, nblk, nsum, nrows
    integer :: isc, ibc, leaddim_Ax, leaddim_Buffer

    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)

    nRHSs      = size(Ax, 2)
    leaddim_Ax = size(Ax, 1)
    leaddim_Buffer = max_blockdim*max_blocks_per_row

    Ax = ZERO

!$OMP PARALLEL PRIVATE(ibr,ibc,isr,isc,nrows,nsum,nblk,j,iRHSs,Buffer,k)
!$OMP DO
    do ibr = 1, blk_nrows
      isr   = kvstr(ibr)
      nrows = kvstr(ibr+1) - isr
      nsum  = 0

      k = ka(ia(ibr))
      do j = ia(ibr), ia(ibr+1)-1
        ibc  = ja(j)
        isc  = kvstc(ibc)
        nblk = kvstc(ibc+1) - isc

!IBM* ASSERT(ITERCNT(16))
        do iRHSs = 1, nRHSs
          ! call ZCOPY(nblk, x(isc,iRHSs), 1, Buffer(nsum+1,iRHSs), 1)
          Buffer(nsum+1:nsum+nblk,iRHSs) = x(isc:isc+nblk-1,iRHSs)
        enddo ! iRHSs

        nsum = nsum + nblk
      enddo ! j

      call ZGEMM('N', 'N', nrows, nRHSs, nsum, ONE, A(k:), nrows, Buffer, leaddim_Buffer, ZERO, Ax(isr,1), leaddim_Ax)

    enddo ! ibr
!$OMP endDO
!$OMP endPARALLEL

  endsubroutine ! vbrmv_mat

endmodule ! vbrmv_mat_mod
