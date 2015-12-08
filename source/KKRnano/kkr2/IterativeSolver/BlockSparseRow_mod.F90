!! this module is a dummy that test how the GPU solver has to be brought into KKRnano
module BlockSparseRow_mod
!! format according to http://docs.nvidia.com/cuda/cusparse/#block-compressed-sparse-row-format-bsr
  implicit none
  private

  public :: BlockSparseRow, create, destroy

#define BSRX

  type BlockSparseRow
    integer :: blockDim = 0 !< == (lmax+1)^2
    integer :: mb = 0 !< number of block rows
    integer :: nb = 0 !< number of block columns !! not in use...
    integer :: nnzb = 0 !< number of non-zero blocks
    double complex, pointer :: bsrValA(:,:,:) !< dim(blockDim,blockDim,nnzb)
    integer, allocatable :: bsrColIndA(:) !< dim(nnzb)
    integer, allocatable :: bsrRowPtrA(:) !< dim(mb+1)
#ifdef BSRX
!!! in the extended Block Compressed Sparse Row Format (BSRX), the
!!! last value of iList in a row is given by bsrEndPtrA(iRow)
    integer, allocatable :: bsrEndPtrA(:) !< dim(mb)
#define  A_bsrEndPtrA(i)  A%bsrEndPtrA(i)
#else
!!! in the default Block Compressed Sparse Row Format (BSR), the
!!! last value of iList in a row is given by bsrRowPtrA(iRow+1)-1
#define  A_bsrEndPtrA(i)  A%bsrRowPtrA(i+1)-1
#endif
  endtype ! BlockSparseRow

  interface create
    module procedure createBSR_from_VBSR
  endinterface
  
  interface destroy
    module procedure destroyBlockSparseRow
  endinterface
  
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
      beta = zero ! instead of an initialization of Avec to zero
      do iList = A%bsrRowPtrA(iRow), A_bsrEndPtrA(iRow)
        iCol = A%bsrColIndA(iList) ! indirection

        ! now: Avec(:,:,iRow) = beta*Avec(:,:,iRow) + matmul(A%bsrValA(:,:,iList), vec(:,:,iCol))

        call zgemm('n', 'n', N, N, N, one, A%bsrValA(:,1,iList), N, vec(:,1,iCol), N, beta, Avec(:,1,iRow), N)

        beta = one ! simply add all further contributions in this row
      enddo ! iList
    enddo ! iRow
    !$omp end parallel do

#undef N
  endsubroutine ! multiply

  subroutine createBSR_from_VBSR(self, lmax, mb, values, ia, ja)
    type(BlockSparseRow), intent(inout) :: self
    integer, intent(in) :: lmax, mb
    double complex, intent(in) :: values(:) !< dim(N*N*nnzb)
    integer, intent(in) :: ia(:) !< dim(mb+1)
    integer, intent(in) :: ja(:) !< dim(nnzb)
    
    integer :: ist
    
    self%blockDim = (lmax + 1)**2
    self%mb = mb
    self%nb = 0 !! not in use
    self%nnzb = size(ja)
    
    allocate(self%bsrValA(self%blockDim,self%blockDim,self%nnzb), &
      self%bsrColIndA(self%nnzb), self%bsrRowPtrA(mb+1), stat=ist)
      
    ! copy data in
    self%bsrValA = reshape(values, [self%blockDim,self%blockDim,self%nnzb])
    
    ! copy index lists
    self%bsrRowPtrA(:) = ia(:)
    self%bsrColIndA(:) = ja(:)
    
#ifdef BSRX
    allocate(self%bsrEndPtrA(mb), stat=ist)
    self%bsrEndPtrA(:) = self%bsrRowPtrA(2:) - 1
#endif
  endsubroutine ! create

  elemental subroutine destroyBlockSparseRow(self)
    type(BlockSparseRow), intent(inout) :: self
    integer :: ist
    self%blockDim = 0
    self%mb = 0
    self%nb = 0
    self%nnzb = 0
    deallocate(self%bsrColIndA, self%bsrRowPtrA, stat=ist)
#ifdef BSRX
    deallocate(self%bsrEndPtrA, stat=ist)
#endif
    deallocate(self%bsrValA, stat=ist)
  
  endsubroutine ! destroy

endmodule ! BlockSparseRow_mod
