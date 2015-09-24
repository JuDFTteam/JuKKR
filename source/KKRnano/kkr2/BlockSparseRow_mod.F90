module BlockSparseRow_mod
!! format according to http://docs.nvidia.com/cuda/cusparse/#block-compressed-sparse-row-format-bsr
  implicit none
  private
  
  public :: BlockSparseRow
  
  type BlockSparseRow
    integer :: blockDim = 0 !< = (lmaxd+1)^2
    integer :: mb = 0 !< number of block rows
    integer :: nb = 0 !< number of block columns
    integer :: nnzb = 0 !< number of non-zero blocks
    double complex, allocatable :: bsrValA(:,:,:) !< dim(blockDim,blockDim,nnzb)
    integer, allocatable :: bsrRowPtrA(:) !< dim(mb+1)
    integer, allocatable :: bsrColIndA(:) !< dim(nnzb)
  endtype
  
  contains
  
  subroutine multiply(A, vec, Avec)
    type(BlockSparseRow), intent(in) :: A
#define N A%blockDim    
    double complex, intent(in)  ::  vec(N,N,*)
    double complex, intent(out) :: Avec(N,N,*)
    
    double complex :: beta
    double complex, parameter :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)
    integer :: iRow, iCol, iList

    !$omp parallel do private(iRow, iList, iCol, beta)
    do iRow = 1, A%mb
      beta = zero ! instead of an initialization of avec to zero
      do iList = A%bsrRowPtrA(iRow), A%bsrRowPtrA(iRow+1)-1
        iCol = A%bsrColIndA(iList) ! indirection
        
        ! now: Avec(:,:,iRow) = beta*Avec(:,:,iRow) + matmul(A%bsrValA(:,:,iList), vec(:,:,iCol))
        call zgemm('n', 'n', N, N, N, one, A%bsrValA(:,1,iList), N, vec(:,1,iCol), N, beta, Avec(:,1,iRow), N)
        
        beta = one ! simply add all further contributions in this row
      enddo ! iList
    enddo ! iRow
    !$omp end parallel do
    
#undef N
  endsubroutine ! multiply
  
endmodule ! BlockSparseRow_mod