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
  public :: bsr_times_mat
  
  contains

  !> Heavily modified routine from SPARSKIT
  subroutine bsr_times_mat(ia, ja, A, x, Ax, nFlops)
    integer, intent(in) :: ia(:) !> dim(A%nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)      ! column indices
    double complex, intent(in)  ::  A(:,:,:) ! dim(blockDim,blockDim,nnzb)
    double complex, intent(in)  ::  x(:,:,:) ! dim(blockDim,nRHSs,nRows)
    double complex, intent(out) :: Ax(:,:,:) ! dim(blockDim,nRHSs,nRows)
    integer(kind=8), intent(inout) :: nFlops

    external :: ZGEMM ! from BLAS
    
    !-----------------------------------------------------------------------
    !     Sparse matrix-full vector product, in VBR format.
    !-----------------------------------------------------------------------
    !     On entry:
    !--------------
    !     blk_nrows      = number of block rows in matrix A
    !     ia,ja,ka,A = matrix A in variable block row format
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

    integer :: ibr, Aind, nRHSs, nRows,leadDim_Ax, leadDim_x, leadDim_A, blocksize
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
    blocksize  = size(A, 2)
    if (leadDim_Ax /= leadDim_A) stop __LINE__
    
    if (leadDim_Ax < blocksize) stop __LINE__
    if (leadDim_A  < blocksize) stop __LINE__
    if (leadDim_x  < blocksize) stop __LINE__

! #define GENERIC_matmul
#ifdef  GENERIC_matmul

    do ibr = 1, nRows
      Ax(:,:,ibr) = ZERO
      do Aind = ia(ibr), ia(ibr + 1) - 1
        Ax(:,:,ibr) = Ax(:,:,ibr) + matmul(A(:,:,Aind), x(:,:,ja(Aind)))
      enddo ! Aind
    enddo ! ibr

#else

!$OMP DO PARALLEL PRIVATE(ibr, beta, Aind) reduction(+:nFlops)
    do ibr = 1, nRows
      beta = ZERO ! instead of Ax(:,:,ibr) = ZERO, we simply set beta = ZERO for the 1st Aind-loop iteration
      do Aind = ia(ibr), ia(ibr + 1) - 1
      
        !     gemm:          N         M,      K          1    A(N,K)       N          B(K,M)           K          1     C(N,M)       N        
        !     here:          N         M,      N          1    A(N,N)       N          B(N,M)           N          1     C(N,M)       N        
        call ZGEMM('n', 'n', blocksize, nRHSs, blocksize, ONE, A(1,1,Aind), leadDim_A, x(1,1,ja(Aind)), leadDim_x, beta, Ax(1,1,ibr), leadDim_Ax)
        
        beta = ONE ! from the 2nd Aind-loop iteration, we need to accumulate onto Ax(:,:,ibr)
      enddo ! Aind
      nFlops = nFlops + (8_8*nRHSs)*(blocksize*blocksize)*(ia(ibr + 1) - ia(ibr))
    enddo ! ibr
!$OMP END PARALLEL DO

#endif
  endsubroutine ! bsr_times_mat

endmodule ! vbrmv_mat_mod
