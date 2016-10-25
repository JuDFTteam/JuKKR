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
  subroutine vbrmv_mat(blk_nrows, ia, ja, A, x, Ax, lmsmax, max_blocks_per_row, nFlops)
    
    integer, intent(in) :: blk_nrows          ! ToDo remove from interface
    integer, intent(in) :: max_blocks_per_row ! ToDo remove from interface
    
    integer, intent(in) :: ia(:) !> dim(blk_nrows + 1)
    integer, intent(in) :: ja(:) !> dim(A%nnzb)
    double complex, intent(in)  ::  A(:,:,:) ! (blockDim,blockDim,nnzb)
    double complex, intent(in)  ::  x(:,:,:) ! (blockDim,nRHSs,nRows)
    double complex, intent(out) :: Ax(:,:,:) ! (blockDim,nRHSs,nRows)
    integer, intent(in) :: lmsmax
    integer(kind=8), intent(inout) :: nFlops
    
    external :: ZGEMM ! from BLAS
    
    !-----------------------------------------------------------------------
    !     Sparse matrix-full vector product, in VBR format.
    !-----------------------------------------------------------------------
    !     On entry:
    !--------------
    !     blk_nrows      = number of block rows in matrix A
    !     ia,ja,ka,A,kvstr = matrix A in variable block row format, kvstc==kvstr since square operator
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

    integer :: ibr, Aind, nRHSs, nRows,leadDim_Ax, leadDim_x, leadDim_A
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    double complex :: beta

    nRows      = size(Ax, 3)
    if (nRows /= size(x, 3)) stop __LINE__
    nRHSs      = size(Ax, 2)
    if (nRHSs /= size(x, 2)) stop __LINE__
    leadDim_Ax = size(Ax, 1)
    leadDim_x  = size(x, 1)
    if (leadDim_Ax /= leadDim_x) stop __LINE__
    leadDim_A  = size(A, 1)
    if (leadDim_Ax /= leadDim_A) stop __LINE__
    
    if (leadDim_Ax < lmsmax) stop __LINE__
    if (leadDim_A  < lmsmax) stop __LINE__
    if (leadDim_x  < lmsmax) stop __LINE__

! #define GENERIC    
  
!$OMP PARALLEL PRIVATE(ibr, beta, Aind) reduction(+:nFlops)
!$OMP DO
    do ibr = 1, nRows
#ifdef GENERIC
      Ax(:,:,ibr) = ZERO
#else      
      beta = ZERO
#endif
      do Aind = ia(ibr), ia(ibr + 1) - 1
#ifdef GENERIC
        Ax(:,:,ibr) = Ax(:,:,ibr) + matmul(A(:,:,Aind), x(:,:,ja(Aind)))
#else
        !     gemm:          N       M,     K       1    A(N,K)       N          B(K,M)           K          1     C(N,M)       N        
        !     here:          N       M,     N       1    A(N,N)       N          B(N,M)           N          1     C(N,M)       N        
        call ZGEMM('n', 'n', lmsmax, nRHSs, lmsmax, ONE, A(1,1,Aind), leadDim_A, x(1,1,ja(Aind)), leadDim_x, beta, Ax(1,1,ibr), leadDim_Ax)
#endif
        nFlops = nFlops + (8_8*nRHSs)*(lmsmax*lmsmax)
        beta = ONE
      enddo ! Aind
    enddo ! ibr
!$OMP END DO
!$OMP END PARALLEL

  endsubroutine ! vbrmv_mat

endmodule ! vbrmv_mat_mod
