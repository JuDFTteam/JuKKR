!! this module is a dummy that tests how the GPU solver has to be brought into KKRnano
module BlockSparseRow_mod
!! format according to http://docs.nvidia.com/cuda/cusparse/#block-compressed-sparse-row-format-bsr
  implicit none
  private

  public :: BlockSparseRow, create, destroy, multiply

! #define BSRX
#define complex_data_t double complex

  type BlockSparseRow
    integer :: blockDim = 0 !< == (lmax+1)^2
    integer :: mb = 0 !< number of block rows
    integer :: nb = 0 !< number of block columns !! not in use...
    integer :: nnzb = 0 !< number of non-zero blocks
    integer, allocatable :: bsrColInd(:) !< dim(nnzb)
    integer, allocatable :: bsrRowPtr(:) !< dim(mb+1)
#ifdef BSRX
!!! in the extended Block Compressed Sparse Row Format (BSRX), the
!!! last value of Aind in a row is given by bsrEndPtr(iRow)
    integer, allocatable :: bsrEndPtr(:) !< dim(mb)
#define  _bsrEndPtr(i)  %bsrEndPtr(i)
#else
!!! in the default Block Compressed Sparse Row Format (BSR), the
!!! last value of Aind in a row is given by bsrRowPtr(iRow+1)-1
#define  _bsrEndPtr(i)  %bsrRowPtr(i+1)-1
#endif

!!! !!! naming in cuSPARSE
!!!       bsrColIndA: bsrColInd
!!!       bsrRowPtrA: bsrRowPtr
!!!       bsrEndPtrA: bsrEndPtr
!!!  complex_data_t, pointer :: bsrValA(:,:,:) !< dim(blockDim,blockDim,nnzb)
!!!       bsrValA:    Val       (is not included into the BlockSparseRow descriptor)
  endtype ! BlockSparseRow

  interface create
    module procedure createBSR_from_VBSR, createBSR_from_full
  endinterface
  
  interface destroy
    module procedure destroyBlockSparseRow
  endinterface

  interface multiply
    module procedure multiplyBSR_to_Bvec, multiplyBSR_to_BSR_truncate
  endinterface
  
  contains

  subroutine multiplyBSR_to_Bvec(A, Aval, vec, Avec)
    !! multiply a blocks sparse operator to a (dense) block vectors
    !! can be used for a single atom per MPI process
    !! could call cuSPARSE GPU kernel
    type(BlockSparseRow), intent(in) :: A
#define M A%blockDim
#define N A%blockDim
#define K A%blockDim
    complex_data_t, intent(in)  :: Aval(K,N,*) ! block list of operator A
    complex_data_t, intent(in)  ::  vec(M,K,*) ! block vector input
    complex_data_t, intent(out) :: Avec(M,N,*) ! block vector result

    ! .. locals ..
    complex_data_t :: beta
    complex_data_t, parameter :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)
    integer :: iRow, jCol, Aind

    !$omp parallel do private(iRow, jCol, Aind, beta)
    do iRow = 1, A%mb
      beta = zero ! instead of an initialization of Avec to zero
      do Aind = A%bsrRowPtr(iRow), A _bsrEndPtr(iRow)
        jCol = A%bsrColInd(Aind) ! indirection

        ! now: Avec(:,:,iRow) = beta*Avec(:,:,iRow) + Aval(:,:,Aind) .times. vec(:,:,jCol)
        ! use gemm('n', 'n', M, N, K, 1.,  A(1:M,1:K),    M, B(1:K,1:N),     K, 1.,   C(1:M,1:N),     M)
        call zgemm('n', 'n', M, N, K, one, vec(:,1,jCol), M, Aval(:,1,Aind), K, beta, Avec(:,1,iRow), M)

        beta = one ! simply add all further contributions in this row
      enddo ! Aind
    enddo ! iRow
    !$omp end parallel do

#undef K
#undef N
#undef M
  endsubroutine ! multiply

  
  subroutine multiplyBSR_to_BSR_truncate(A, Aval, G, Gval, Rval, GiFlop)
    type(BlockSparseRow), intent(in) :: A !< Matrix to be inverted
    type(BlockSparseRow), intent(in) :: G !< Green function
!!! type(BlockSparseRow), intent(in) :: R !< Result operator (argument is not needed since the structure of R is that of G)
#define R G
#define M G%blockDim
#define N A%blockDim
#define K A%blockDim
    complex_data_t, intent(in)  :: Aval(K,N,*) ! input  block list of (square) operator A
    complex_data_t, intent(in)  :: Gval(M,K,*) ! input  block list 
    complex_data_t, intent(out) :: Rval(M,N,*) ! result block list, assume R to have the same structure as G
    real, intent(out), optional :: GiFlop

    ! .. locals ..
    complex_data_t, parameter :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)
    integer :: iRow, jCol, kCol, Aind, Rind, Gind, indG
    integer(kind=8) :: nBlockOps

    nBlockOps = 0

    !$omp parallel
    
    !$omp workshare
    Rval(:,:,:R%nnzb) = zero ! initialization to zero
    !$omp end workshare
    
    !$omp do private(iRow, jCol, kCol, Aind, Rind, Gind, indG) reduction(+:nBlockOps)
    do iRow = 1, A%mb
      do Aind = A%bsrRowPtr(iRow), A _bsrEndPtr(iRow)
        jCol = A%bsrColInd(Aind) ! indirection
        ! for each matrix block element A_full(jCol,iRow)
        do Rind = R%bsrRowPtr(iRow), R _bsrEndPtr(iRow)
          kCol = R%bsrColInd(Rind) ! indirection
          ! update matrix block element R_full(kCol,iRow)

#define   jRow jCol
          ! find out, if G_full(jRow,kCol) exists
          indG = -1 ! init as invalid index
          do Gind = G%bsrRowPtr(jRow), G _bsrEndPtr(jRow) ! search
            if (G%bsrColInd(Gind) == kCol) indG = Gind
          enddo ! Gind
          if (indG > -1) then ! yes

            ! now: Rval(:,:,Rind) = one*Rval(:,:,Rind) + Aval(:,:,Aind) .times. Gval(:,:,indG)
            ! use gemm('n', 'n', M, N, K, 1.,  A(1:M,1:K),     M, B(1:K,1:N),     K, 1.,  C(1:M,1:N),     M)
            call zgemm('n', 'n', M, N, K, one, Gval(:,1,indG), M, Aval(:,1,Aind), K, one, Rval(:,1,Rind), M)

            nBlockOps = nBlockOps + 1

          endif ! G(jRow,kCol) exists
#undef    jRow

        enddo ! Rind
      enddo ! Aind
    enddo ! iRow
    !$omp end do
    
    !$omp end parallel
    
    if(present(GiFlop)) GiFlop = nBlockOps*(M*8.*N*.5d0**30*K) ! assume a complex data_t, so each FMA has 8 Flop
#undef K
#undef M
#undef N
#undef R
  endsubroutine ! multiply

  
  subroutine createBSR_from_VBSR(self, lmax, mb, values, ia, ja, bsrVal)
    !! get the data from a variable block row format
    type(BlockSparseRow), intent(inout) :: self
    integer, intent(in) :: lmax, mb
    complex_data_t, intent(in) :: values(:) !< dim(N*N*nnzb)
    integer, intent(in) :: ia(:) !< dim(mb+1)
    integer, intent(in) :: ja(:) !< dim(nnzb)
    complex_data_t, allocatable, intent(inout), optional :: bsrVal(:,:,:) !< dim(N,N,nnz) ! warning: allocation inside this routine
    
    integer :: ist
    
    self%blockDim = (lmax + 1)**2
    self%mb = mb
    self%nb = 0 !! not in use
    self%nnzb = size(ja)
    
    if (present(bsrVal)) then
      deallocate(bsrVal, stat=ist)
      allocate(bsrVal(self%blockDim,self%blockDim,self%nnzb), stat=ist)
      ! copy data in
      bsrVal = reshape(values, [self%blockDim,self%blockDim,self%nnzb])
    endif

    allocate(self%bsrColInd(self%nnzb), self%bsrRowPtr(mb+1), stat=ist)
    
    ! copy index lists
    self%bsrRowPtr(:) = ia(:)
    self%bsrColInd(:) = ja(:)
    
#ifdef BSRX
    allocate(self%bsrEndPtr(self%mb), stat=ist)
    self%bsrEndPtr(:) = self%bsrRowPtr(2:) - 1
#endif
  endsubroutine ! create

  
  subroutine createBSR_from_full(self, values, bsrVal)
    !! get the data from a dense block-matrix
    type(BlockSparseRow), intent(inout) :: self
    complex_data_t, intent(in) :: values(:,:,:,:) !< dim(N,N,nb,mb)
    complex_data_t, allocatable, intent(inout), optional :: bsrVal(:,:,:) !< dim(N,N,nnz) ! warning: allocation inside this routine
    
    integer :: ist, iRow, Aind, jCol
    logical(kind=1), allocatable :: nnz(:,:)
    
    self%blockDim = size(values, 1)
    if (self%blockDim /= size(values, 2)) stop 'createBSR_from_full: requested shape values(blockDim,blockDim,mb,*)'
    self%nb = size(values, 3) !! number of columns self%nb is not in use
    self%mb = size(values, 4) !! number of rows
    allocate(nnz(self%nb,self%mb))
    nnz(:,:) = any(any(values /= 0.d0, dim=1), dim=1)
    self%nnzb = count(nnz)
    
    if (present(bsrVal)) then
      deallocate(bsrVal, stat=ist)
      allocate(bsrVal(self%blockDim,self%blockDim,self%nnzb), stat=ist)
    endif
    
    allocate(self%bsrColInd(self%nnzb), self%bsrRowPtr(self%mb+1), stat=ist)

    Aind = 0
    self%bsrRowPtr(1) = 1
    do iRow = 1, self%mb
      do jCol = 1, self%nb
        if (nnz(jCol,iRow)) then
          Aind = Aind + 1
          ! copy data in
          if (present(bsrVal)) bsrVal(:,:,Aind) = values(:,:,jCol,iRow)
          ! create index list
          self%bsrColInd(Aind) = jCol
          self%bsrRowPtr(iRow+1) = Aind + 1
        endif ! non-zero
      enddo ! jCol     
    enddo ! iRow
    if (Aind /= self%nnzb) stop 'createBSR_from_full: fatal counting error!'

#ifdef BSRX
    allocate(self%bsrEndPtr(self%mb), stat=ist)
    self%bsrEndPtr(:) = self%bsrRowPtr(2:) - 1
#endif
  endsubroutine ! create
  
  
  elemental subroutine destroyBlockSparseRow(self)
    type(BlockSparseRow), intent(inout) :: self
    integer :: ist
    self%blockDim = 0
    self%mb = 0
    self%nb = 0
    self%nnzb = 0
    deallocate(self%bsrColInd, self%bsrRowPtr, stat=ist)
#ifdef BSRX
    deallocate(self%bsrEndPtr, stat=ist)
#endif

  endsubroutine ! destroy

endmodule ! BlockSparseRow_mod


#ifdef TESTMAIN_BlockSparseRow
!+ TESTMAIN_BlockSparseRow

!!!
!!!>  ifort -warn -check all -O0 -g -mkl -D TESTMAIN_BlockSparseRow IterativeSolver/BlockSparseRow_mod.F90  && ./a.out 1024 123
!!!>  gfortran -ffree-line-length-0 -g -D TESTMAIN_BlockSparseRow IterativeSolver/BlockSparseRow_mod.F90 -lblas && ./a.out 1024 123
!!!
program test_bsr
  use BlockSparseRow_mod !, only:
implicit none
  character(len=8) :: CLarg(0:3)
  integer, parameter :: ShowR=0, ShowG=0, ShowH=0, Hfill=16, N = 1 ! BlockSize
  integer :: ilen, ios, iarg, Nrows, Mcols, iRow, jCol, Rind, kCol, ip, nerror(19)=0
  double precision, parameter :: Gfill=0.5, point=.5d0**10
  double precision :: col, elem
#define real_data_t double precision
#define gemm DGEMM
  external :: gemm ! BLAS matrix matrix multiplication
  real_data_t, parameter :: one = 1.0, sero = 0.0
  real_data_t, allocatable :: Hfull(:,:), Gfull(:,:), Rfull(:,:)
  complex_data_t, allocatable :: Rval(:,:,:), Gval(:,:,:), Hval(:,:,:)
  type(BlockSparseRow) :: H, G ! operators
  real :: GiFlop

  do iarg = 0, ubound(CLarg, 1)
    call get_command_argument(iarg, CLarg(iarg), ilen, ios)
  enddo ! iarg
  read(unit=CLarg(1), fmt=*) Nrows
  read(unit=CLarg(2), fmt=*) Mcols

#define Krows Nrows
  allocate(Rfull(Mcols,Nrows), Gfull(Mcols,Krows), Hfull(Krows,Nrows)) ! R:=G*H
  Hfull = sero ; Gfull = sero ; Rfull = sero
  
  do iRow = 1, Nrows
    ! fill H
    do while(count(Hfull(:,iRow) /= 0) < min(Hfill, Nrows))
      call random_number(col) ; jCol = ceiling(col*Krows)
!     if(ShowH>0) write(*, fmt="(A,9I6)") 'rand H',jCol,iRow
      Hfull(jCol,iRow) = jCol + point * iRow 
    enddo ! while
    if(ShowH>0) write(*, fmt="(A,I4,99F9.4)") 'H',jCol,Hfull(1:min(Krows, 8),jCol)
    ! fill G
    do while(count(Gfull(:,iRow) /= 0) < Gfill*Mcols)
      call random_number(col) ; jCol = ceiling(col*Mcols)
!     if(ShowG>0) write(*, fmt="(A,9I6)") 'rand G',jCol,iRow
      Gfull(jCol,iRow) = jCol + point * iRow 
    enddo ! while
    if(ShowG>0) write(*, fmt="(A,I4,99F9.4)") 'G',iRow,Gfull(1:min(Mcols, 8),iRow)
  enddo ! iCol

  ! create reference A an m x k matrix, B a k x n matrix and C an m x n matrix
  !                   M      N      K           A             B                   C
  call gemm('n', 'n', Mcols, Nrows, Krows, one, Gfull, Mcols, Hfull, Krows, sero, Rfull, Mcols) ! BLAS routine
  !!! do j = 1, n ;  do l = 1, k ;  do i = 1, m ;   c(i,j) = c(i,j) + b(l,j)*a(i,l)   ; enddo ; enddo ; enddo
  ! do j = 1, Nrows ;  do l = 1, Nrows ;  do i = 1, Mcols ; Rfull(i,j) = Rfull(i,j) + Hfull(l,j)*Gfull(i,l) ; enddo ; enddo ; enddo

  nerror = 0
  do iRow = 1, Nrows
    if(ShowR>0) write(*, fmt="(A,I4,99F9.4)") 'R',iRow,Rfull(1:min(Mcols, 8),iRow)/Nrows
    do kCol = 1, Mcols
      elem = dot_product(Hfull(:,iRow), Gfull(kCol,:)) ! simple evaluation of a single element of a matrix-matrix product
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Rfull(kCol,iRow)) > .1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! kCol
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  call create(G, dcmplx(reshape(Gfull, [1,1,Mcols,Krows])), bsrVal=Gval)
  call create(H, dcmplx(reshape(Hfull, [1,1,Krows,Nrows])), bsrVal=Hval)

  !! use the sparse structure of G for R
#define R G
  allocate(Rval(1,1,R%nnzb))

  ! API: multiply(A, Aval, G, Gval, Rval, GiFlop)
  call multiply(H, Hval, G, Gval, Rval, GiFlop=GiFlop)
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop'

  nerror = 0
  do iRow = 1, R%mb
    do Rind = R%bsrRowPtr(iRow), R _bsrEndPtr(iRow)
      kCol = R%bsrColInd(Rind) ! indirection
      elem = dreal(Rval(1,1,Rind)) ! take only the real part
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Rfull(kCol,iRow)) > .1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! Rind
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  deallocate(Rval)
#undef R  
  call destroy(H)
  call destroy(G)
#undef Krows
endprogram ! test

!- TESTMAIN_BlockSparseRow
#endif
