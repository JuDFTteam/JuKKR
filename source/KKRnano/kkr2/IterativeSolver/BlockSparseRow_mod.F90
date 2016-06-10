!! this module is a dummy that tests how the GPU solver has to be brought into KKRnano
module BlockSparseRow_mod
!! format according to http://docs.nvidia.com/cuda/cusparse/#block-compressed-sparse-row-format-bsr
  implicit none
  private

  public :: BlockSparseRow, create, destroy, multiply

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
!!! last value of Aind in a row is given by bsrEndPtrA(iRow)
    integer, allocatable :: bsrEndPtrA(:) !< dim(mb)
#define  _bsrEndPtrA(i)  %bsrEndPtrA(i)
#else
!!! in the default Block Compressed Sparse Row Format (BSR), the
!!! last value of Aind in a row is given by bsrRowPtrA(iRow+1)-1
#define  _bsrEndPtrA(i)  %bsrRowPtrA(i+1)-1
#endif
  endtype ! BlockSparseRow

  interface create
    module procedure createBSR_from_VBSR, createBSR_from_full
  endinterface
  
  interface destroy
    module procedure destroyBlockSparseRow
  endinterface

  interface multiply
    module procedure multiplyBSR_to_Bvec, multiplyBSR_to_BSR
  endinterface
  
  contains

  subroutine multiplyBSR_to_Bvec(A, vec, Avec)
    !! multiply a blocks sparse operator to a (dense) block vectors
    !! can be used for a single atom per MPI process
    type(BlockSparseRow), intent(in) :: A
#define N A%blockDim
    double complex, intent(in)  ::  vec(N,N,*) ! block vector
    double complex, intent(out) :: Avec(N,N,*) ! block vector

    double complex :: beta
    double complex, parameter :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)
    integer :: iRow, jCol, Aind

    !$omp parallel do private(iRow, Aind, jCol, beta)
    do iRow = 1, A%mb
      beta = zero ! instead of an initialization of Avec to zero
      do Aind = A%bsrRowPtrA(iRow), A _bsrEndPtrA(iRow)
        jCol = A%bsrColIndA(Aind) ! indirection

        ! now: Avec(:,:,iRow) = beta*Avec(:,:,iRow) + matmul(A%bsrValA(:,:,Aind), vec(:,:,jCol))

        call zgemm('n', 'n', N, N, N, one, A%bsrValA(:,1,Aind), N, vec(:,1,jCol), N, beta, Avec(:,1,iRow), N)

        beta = one ! simply add all further contributions in this row
      enddo ! Aind
    enddo ! iRow
    !$omp end parallel do

#undef N
  endsubroutine ! multiply

  
  subroutine multiplyBSR_to_BSR(A, G, Gval, Rval, GiFlop)
    type(BlockSparseRow), intent(in) :: A !< Matrix to be inverted
    type(BlockSparseRow), intent(in) :: G !< Green function: Gval is used instead of G%bsrValA
#define N A%blockDim
    double complex, intent(in)  :: Gval(N,N,*) ! input  block vector
    double complex, intent(out) :: Rval(N,N,*) ! result block vector, assume R to have the same structure as G
#define R G
    real, intent(out), optional :: GiFlop

    double complex, parameter :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)
    integer :: iRow, jCol, kCol, Aind, Rind, Gind, indG
    integer(kind=8) :: nBlockOps
    
    nBlockOps = 0

    !$omp parallel
    
    !$omp workshare
    Rval(:,:,:R%nnzb) = zero ! initialization to zero
    !$omp end workshare
    
    !$omp do private(iRow, jCol, kCol, Aind, Rind, Gind, indG) reduction(+:nBlockOps)
    do iRow = 1, A%mb
      do Aind = A%bsrRowPtrA(iRow), A _bsrEndPtrA(iRow)
        jCol = A%bsrColIndA(Aind) ! indirection
        do Rind = R%bsrRowPtrA(iRow), R _bsrEndPtrA(iRow)
          kCol = R%bsrColIndA(Rind) ! indirection

#define   jRow jCol
          ! find out, if G(jRow,kCol) exists
          indG = -1
          do Gind = G%bsrRowPtrA(jRow), G _bsrEndPtrA(jRow)
            if (G%bsrColIndA(Gind) == kCol) indG = Gind
          enddo ! Gind
          if (indG > -1) then

            ! now: Rval(:,:,Rind) = one*Rval(:,:,Rind) + matmul(A%bsrValA(:,:,Aind), Gval(:,:,indG))
            call zgemm('n', 'n', N, N, N, one, A%bsrValA(:,1,Aind), N, Gval(:,1,indG), N, one, Rval(:,1,Rind), N)
            nBlockOps = nBlockOps + 1

          endif ! G(jRow,kCol) exists
#undef    jRow

        enddo ! Rind
      enddo ! Aind
    enddo ! iRow
    !$omp end do
    
    !$omp end parallel
    if(present(GiFlop)) GiFlop = nBlockOps*(N*.5d0**10)**3*8.
#undef R
#undef N
  endsubroutine ! multiply
  
  
  subroutine createBSR_from_VBSR(self, lmax, mb, values, ia, ja)
    !! get the data from a variable block row format
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
    allocate(self%bsrEndPtrA(self%mb), stat=ist)
    self%bsrEndPtrA(:) = self%bsrRowPtrA(2:) - 1
#endif
  endsubroutine ! create

  
  subroutine createBSR_from_full(self, values)
    !! get the data from a dense block-matrix
    type(BlockSparseRow), intent(inout) :: self
    double complex, intent(in) :: values(:,:,:,:) !< dim(N,N,nb,mb)
    
    integer :: ist, iRow, Aind, jCol
    logical(kind=1), allocatable :: nnz(:,:)
    
    self%blockDim = size(values, 1)
    if (self%blockDim /= size(values, 2)) stop 'createBSR_from_full: requested shape values(blockDim,blockDim,mb,*)'
    self%nb = size(values, 3) !! number of columns self%nb is not in use
    self%mb = size(values, 4) !! number of rows
    allocate(nnz(self%nb,self%mb))
    nnz(:,:) = any(any(values /= 0.d0, dim=1), dim=1)
    self%nnzb = count(nnz)
    
    allocate(self%bsrValA(self%blockDim,self%blockDim,self%nnzb), &
      self%bsrColIndA(self%nnzb), self%bsrRowPtrA(self%mb+1), stat=ist)

    Aind = 0
    self%bsrRowPtrA(1) = 1
    do iRow = 1, self%mb
      do jCol = 1, self%nb
        if (nnz(jCol,iRow)) then
          Aind = Aind + 1
          ! copy data in
          self%bsrValA(:,:,Aind) = values(:,:,jCol,iRow)
          ! create index lists
          self%bsrColIndA(Aind) = jCol
          self%bsrRowPtrA(iRow+1) = Aind + 1
        endif ! non-zero
      enddo ! jCol     
    enddo ! iRow
    if (Aind /= self%nnzb) stop 'createBSR_from_full: fatal counting error!'

#ifdef BSRX
    allocate(self%bsrEndPtrA(self%mb), stat=ist)
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


#ifdef TESTMAIN_BlockSparseRow
!+ TESTMAIN_BlockSparseRow

!!!
!!!>  ifort -warn -check -g -mkl -D TESTMAIN_BlockSparseRow IterativeSolver/BlockSparseRow_mod.F90  && ./a.out 1024 64
!!!
program test_bsr
  use BlockSparseRow_mod !, only:
implicit none
  character(len=8) :: CLarg(0:3)
  integer, parameter :: ShowR=0, ShowG=0, ShowH=0, Hfill=16, N = 1 ! BlockSize
  integer :: ilen, ios, iarg, Nrows, Mcols, iRow, jCol, Rind, kCol, ip, nerror(19)=0
! integer :: l, i, j
  double precision, parameter :: Gfill=0.5, point=.5d0**10
  double precision :: col, elem
#define dp 8
#define gemm DGEMM
  external :: gemm ! BLAS matrix matrix multiplication
  real(kind=dp), parameter :: one = 1.0, sero = 0.0
  real(kind=dp), allocatable :: Hfull(:,:), Gfull(:,:), Rfull(:,:)
  double complex, allocatable :: Rval(:,:,:)
  type(BlockSparseRow) :: H, G ! operators
  real :: GiFlop
  
  do iarg = 0, ubound(CLarg, 1)
    call get_command_argument(iarg, CLarg(iarg), ilen, ios)
  enddo ! iarg
  read(unit=CLarg(1), fmt=*) Nrows
  read(unit=CLarg(2), fmt=*) Mcols
  
  allocate(Rfull(Mcols,Nrows), Gfull(Mcols,Nrows), Hfull(Nrows,Nrows)) ! R:=G*H
  Hfull = sero ; Gfull = sero ; Rfull = sero
  
  ! fill H
  do iRow = 1, Nrows
    do while(count(Hfull(:,iRow) /= 0) < min(Hfill, Nrows))
      call random_number(col) ; jCol = ceiling(col*Nrows)
!     if(ShowH) write(*, fmt="(A,9I6)") 'rand H',jCol,iRow
      Hfull(jCol,iRow) = jCol + point * iRow 
    enddo ! while
!   write(*, fmt="(A,I4,99F9.4)") 'H',jCol,Hfull(1:min(Nrows, 8),jCol)
  ! fill G
    do while(count(Gfull(:,iRow) /= 0) < Gfill*Mcols)
      call random_number(col) ; jCol = ceiling(col*Mcols)
!     if(ShowG) write(*, fmt="(A,9I6)") 'rand G',jCol,iRow
      Gfull(jCol,iRow) = jCol + point * iRow 
    enddo ! while
!   write(*, fmt="(A,I4,99F9.4)") 'G',kCol,Gfull(1:min(Nrows, 8),kCol)
  enddo ! iCol

  ! create reference A an m x k matrix, B a k x n matrix and C an m x n matrix
  !                   M      N      K           A             B                   C
  call gemm('n', 'n', Mcols, Nrows, Nrows, one, Gfull, Mcols, Hfull, Nrows, sero, Rfull, Mcols) ! BLAS routine
  !  do j = 1, n ;  do l = 1, k ;  do i = 1, m ;   c(i,j) = c(i,j) + b(l,j)*a(i,l)
! do j = 1, Nrows ;  do l = 1, Nrows ;  do i = 1, Mcols ; Rfull(i,j) = Rfull(i,j) + Hfull(l,j)*Gfull(i,l) ; enddo ; enddo ; enddo

  nerror = 0
  do iRow = 1, Nrows
    if(ShowR) write(*, fmt="(A,I4,99F9.4)") 'R',iRow,Rfull(1:min(Mcols, 8),iRow)/Nrows
    do kCol = 1, Mcols
      elem = dot_product(Hfull(:,iRow), Gfull(kCol,:)) ! simple evaluation of a single element of a matrix-matrix product
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Rfull(kCol,iRow)) > 0.1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! kCol
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  call create(G, dcmplx(reshape(Gfull, [1,1,Mcols,Nrows])))
  call create(H, dcmplx(reshape(Hfull, [1,1,Nrows,Nrows])))

  allocate(Rval(1,1,G%nnzb))

  
  call multiply(H, G, G%bsrValA, Rval, GiFlop=GiFlop)
  write(*,"(9(A,F0.3))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop'

  nerror = 0
  do iRow = 1, G%mb
    do Rind = G%bsrRowPtrA(iRow), G _bsrEndPtrA(iRow)
      kCol = G%bsrColIndA(Rind) ! indirection
      elem = dreal(Rval(1,1,Rind)) ! take only the real part
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Rfull(kCol,iRow)) > 0.1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! Rind
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  call destroy(H)
  call destroy(G)
 
endprogram ! test

!- TESTMAIN_BlockSparseRow
#endif
