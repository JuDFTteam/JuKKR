!! this module is a dummy that tests how the GPU solver has to be brought into KKRnano
module BlockSparseRow_mod
!! format according to http://docs.nvidia.com/cuda/cusparse/#block-compressed-sparse-row-format-bsr
  implicit none
  private

  public :: BlockSparseRow, create, destroy, multiply
  public :: MultBSRplan

! #define BSRX
#define complex_data_t double complex

  type BlockSparseRow
    integer :: fastBlockDim = 0 !< == (lmax+1)^2
    integer :: slowBlockDim = 0 !< == (lmax+1)^2
    integer :: mb = 0 !< number of block rows
    integer :: nb = 0 !< number of block columns !! not in use...
    integer :: nnzb = 0 !< number of non-zero blocks
    integer, allocatable :: bsrRowPtr(:) !< dim(mb+1)
    integer, allocatable :: bsrColInd(:) !< dim(nnzb)
#ifdef BSRX
!!! in the extended Block Compressed Sparse Row Format (BSRX), the
!!! last value of Bind in a row is given by bsrEndPtr(iRow)
    integer, allocatable :: bsrEndPtr(:) !< dim(mb)
#define  _bsrEndPtr(i)  %bsrEndPtr(i)
#else
!!! in the default Block Compressed Sparse Row Format (BSR), the
!!! last value of Bind in a row is given by bsrRowPtr(iRow+1)-1
#define  _bsrEndPtr(i)  %bsrRowPtr(i+1)-1
#endif

!!! !!! naming in cuSPARSE
!!!       bsrColIndA: bsrColInd
!!!       bsrRowPtrA: bsrRowPtr
!!!       bsrEndPtrA: bsrEndPtr
!!!  complex_data_t, pointer :: bsrValA(:,:,:) !< dim(fastBlockDim,slowBlockDim,nnzb)
!!!       bsrValA:    Val       (is not included into the BlockSparseRow descriptor)
  endtype ! BlockSparseRow

  interface create
    module procedure createBSR_from_VBSR, createBSR_from_full, create_multBRSplan
  endinterface
  
  interface multiply
    module procedure multiplyBSR_to_Bvec, multiplyBSR_to_BSR_truncate, multiplyBSR_to_BSR_planned
  endinterface

  interface destroy
    module procedure destroyBlockSparseRow, destroyMultBRSplan
  endinterface

  type MultBSRplan
    integer :: N, M, K, nCij, nBlockOps
    integer(kind=4), allocatable :: CindNelemStart(:,:) ! dim(4,nCij)
    integer(kind=4), allocatable :: AindBind(:,:) ! dim(2,nBlockOps)
  endtype

  complex_data_t, parameter, private :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)

  contains

  subroutine multiplyBSR_to_Bvec(B, Bval, M, Avec, Cvec)
    !! multiply a blocks sparse operator to a (dense) block vectors
    !! can be used for a single atom per MPI process
    !! could call cuSPARSE GPU kernel
    type(BlockSparseRow), intent(in) :: B
#define K B%fastBlockDim
#define N B%slowBlockDim
    complex_data_t, intent(in)  :: Bval(K,N,*) ! block list of block sparse operator B
    integer, intent(in) :: M
    complex_data_t, intent(in)  :: Avec(M,K,*) ! block vector input
    complex_data_t, intent(out) :: Cvec(M,N,*) ! block vector result

    ! .. locals ..
    complex_data_t :: beta
    integer :: iRow, jCol, Bind

    !$omp parallel do private(iRow, jCol, Bind, beta)
    do iRow = 1, B%mb
      beta = zero ! instead of an initialization of Avec to zero
      do Bind = B%bsrRowPtr(iRow), B _bsrEndPtr(iRow)
        jCol = B%bsrColInd(Bind)

        ! now: Cvec(:,:,iRow) = beta*Cvec(:,:,iRow) + Avec(:,:,jCol) .times. Aval(:,:,Bind) 
        ! use gemm('n', 'n', M, N, K, 1.,  A(1:M,1:K),     M, B(1:K,1:N),     K, 1.,   C(1:M,1:N),     M)
        call zgemm('n', 'n', M, N, K, one, Avec(:,1,jCol), M, Bval(:,1,Bind), K, beta, Cvec(:,1,iRow), M)

        beta = one ! simply add all further contributions in this row
      enddo ! Bind
    enddo ! iRow
    !$omp end parallel do

#undef K
#undef N
  endsubroutine ! multiply

  
  subroutine multiplyBSR_to_BSR_truncate(A, Aval, B, Bval, Cval, GiFlop)
    type(BlockSparseRow), intent(in) :: A !< Green function
    type(BlockSparseRow), intent(in) :: B !< Matrix to be inverted
!!! type(BlockSparseRow), intent(in) :: C !< Result operator (argument is not needed since the structure of C is that of A)
#define C A
#define M A%fastBlockDim
#define N B%fastBlockDim
#define K B%fastBlockDim
    complex_data_t, intent(in)  :: Aval(M,K,*) ! input  block list 
    complex_data_t, intent(in)  :: Bval(K,N,*) ! input  block list of (square) operator B
    complex_data_t, intent(out) :: Cval(M,N,*) ! result block list, assume C to have the same structure as A
    real, intent(out), optional :: GiFlop

    ! .. locals ..
    integer :: iRow, jCol, kCol, Aind, Bind, Cind
    integer(kind=8) :: nBlockOps

    nBlockOps = 0

    !$omp parallel
    
    !$omp workshare
    Cval(:,:,:C%nnzb) = zero ! initialization to zero
    !$omp end workshare
    
    !$omp do private(iRow, jCol, kCol, Bind, Cind, Aind) reduction(+:nBlockOps)
    do iRow = 1, B%mb

#ifndef REUSE_B
      ! reuse elements of C, i.e. keep the accumulator in the cache
      do   Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow) ; kCol = C%bsrColInd(Cind) ! update   matrix block element C_full(kCol,iRow)
#endif
      do   Bind = B%bsrRowPtr(iRow), B _bsrEndPtr(iRow) ; jCol = B%bsrColInd(Bind) ! for each matrix block element B_full(jCol,iRow)
#ifdef  REUSE_B
        do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow) ; kCol = C%bsrColInd(Cind) !   update matrix block element C_full(kCol,iRow)
#endif

#define   jRow jCol
          Aind = exists(A, jRow, kCol) ! find out, if A_full(jRow,kCol) exists
          if (Aind > -1) then ! yes

            ! now: Cval(:,:,Cind) = one*Cval(:,:,Cind) + Aval(:,:,Aind) .times. Aval(:,:,Bind) 
            ! use gemm('n', 'n', M, N, K, 1.,  A(1:M,1:K),     M, B(1:K,1:N),     K, 1.,  C(1:M,1:N),     M)
            call zgemm('n', 'n', M, N, K, one, Aval(:,1,Aind), M, Bval(:,1,Bind), K, one, Cval(:,1,Cind), M)

            nBlockOps = nBlockOps + 1

          endif ! A(jRow,kCol) exists
#undef    jRow

        enddo ! Bind or Cind
      enddo ! Cind or Bind

    enddo ! iRow
    !$omp end do
    
    !$omp end parallel
    
    if (present(GiFlop)) GiFlop = nBlockOps*(M*8.*N*.5d0**30*K) ! assume a complex data_t, so each FMA has 8 Flop
#undef K
#undef M
#undef N
#undef C
  endsubroutine ! multiply

  
  integer function exists(bsr, row, col) result(ind)
    type(BlockSparseRow), intent(in) :: bsr
    integer, intent(in) :: row, col
    integer :: i ! local
    ind = -1
    do i = bsr%bsrRowPtr(row), bsr _bsrEndPtr(row)
      if (bsr%bsrColInd(i) /= col) cycle
      ind = i 
      return
    enddo ! i
  endfunction

  
  subroutine multiplyBSR_to_BSR_planned(p, Aval, Bval, Cval)
    type(MultBSRplan), intent(in) :: p !> plan
    complex_data_t, intent(in)  :: Aval(p%M,p%K,*) ! input  block list 
    complex_data_t, intent(in)  :: Bval(p%K,p%N,*) ! input  block list of (square) operator B
    complex_data_t, intent(out) :: Cval(p%M,p%N,*) ! result block list

    ! .. locals ..
    complex_data_t :: beta
    integer :: iCij, start, ielem, Aind, Bind, Cind, Nelem

    ! these asserts work only with deferred shape arrays
    ! assert(all(shape(Aval) == [p%M,p%K,A%nnzb]))
    ! assert(all(shape(Bval) == [p%K,p%N,B%nnzb]))
    ! assert(all(shape(Cval) == [p%M,p%N,C%nnzb]))
    
    !$omp parallel
    
    !$omp do private(iCij, start, ielem, Aind, Bind, Cind, Nelem, beta)
    do iCij = 1, p%nCij
      beta = zero ! initialization to zero omitted
      Cind  = p%CindNelemStart(1,iCij) ! if we do not want to reorder, Cind could be equal to iCij
      Nelem = p%CindNelemStart(2,iCij)
      start = p%CindNelemStart(3,iCij) ! == exclusive_scan on Nelem(iCij) 
      do ielem = 1, Nelem
        Aind = p%AindBind(1,start + ielem)
        Bind = p%AindBind(2,start + ielem)

        ! now: Cval(:,:,Cind) = one*Cval(:,:,Cind) + Aval(:,:,Aind) .times. Bval(:,:,Bind)
        ! use gemm('n', 'n',   M,   N,   K,  1., A(1:M,1:K),       M, B(1:K,1:N),       K, beta, C(1:M,1:N),       M)
        call zgemm('n', 'n', p%M, p%N, p%K, one, Aval(:,1,Aind), p%M, Bval(:,1,Bind), p%K, beta, Cval(:,1,Cind), p%M)

        beta = one
      enddo ! ielem
    enddo ! iCij
    !$omp end do
    
    !$omp end parallel
    
  endsubroutine ! multiply

  
  
  
  
  
  subroutine fusedMultiplyAdd_BSR(C, Cval, Aval, diag, Bval, GiFlop) ! C = A*diag + B 
    type(BlockSparseRow), intent(in) :: C !< operator structure, same for A and B assumed
    complex_data_t, intent(out) :: Cval(:,:,:) ! result block list dim(fastBlockDim,slowBlockDim,nnzb)
    complex_data_t, intent(in)  :: Aval(:,:,:) ! input  block list dim(fastBlockDim,slowBlockDim,nnzb)
    complex_data_t, intent(in)  :: diag(:,:)   ! diagonal matrix,  dim(fastBlockDim,nb)
    complex_data_t, intent(in), optional :: Bval(:,:,:) ! input    dim(fastBlockDim,slowBlockDim,nnzb)
    real, intent(out), optional :: GiFlop

    ! .. locals ..
    integer :: iRow, jCol, Cind, islow
    complex_data_t :: bvec(C%fastBlockDim)

      if (any(shape(Cval) /= [C%fastBlockDim, C%slowBlockDim, C%nnzb])) stop __LINE__
    if (present(Bval)) then
      if (any(shape(Cval) /= shape(Bval))) stop __LINE__
    endif
      if (any(shape(Cval) /= shape(Aval))) stop __LINE__
 
      if (any(shape(diag) /= [C%fastBlockDim, C%nb]) stop __LINE__ 
      ! there could also be other variants: C%slowBlockDim, C%mb, i.e. possible 4 combination
      ! for those that use chack against C%nb and use jCol, we can flatten the loop into:
      !     do Cind = 1, C _bsrEndPtr(C%mb) ; jCol = C%bsrColInd(Cind) ; doIt ; enddo
      ! or  do Cind = 1, C%nnzb ; jCol = C%bsrColInd(Cind) ; doIt ; enddo

    !$omp parallel

    !$omp do private(iRow, jCol, Cind, islow, bvec)
    do iRow = 1, C%mb
      bvec(:) = 0
      do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow) ; jCol = C%bsrColInd(Cind)

          do islow = 1, C%slowBlockDim
            if (present(Bval)) bvec(:) = Bval(:,islow,Cind)
            Cval(:,islow,Cind) = Aval(:,islow,Cind)*diag(:,jCol) + bvec(:)
          enddo ! islow

      enddo ! Cind
    enddo ! iRow
    !$omp end do
    
    !$omp end parallel
    
    if (present(GiFlop)) GiFlop = C%fastBlockDim*8.*C%slowBlockDim*.5d0**30*C%nnzb ! assume a complex data_t, so each FMA has 8 Flop
  endsubroutine ! fusedMultiplyAdd
 
  
  
  
  subroutine create_multBRSplan(p, A, B, C, GiFlop)
    type(MultBSRplan), intent(inout) :: p !> plan
    type(BlockSparseRow), intent(in) :: A !< Green function
    type(BlockSparseRow), intent(in) :: B !< Matrix to be inverted 
    type(BlockSparseRow), intent(in) :: C !< Result operator
    real, intent(out), optional :: GiFlop

    ! .. locals ..
    integer :: iRow, jCol, kCol, Bind, Cind, Aind, i01, ist, Nelem, nCij, Start
    integer(kind=8) :: nBlockOps

    p%M = A%fastBlockDim
    p%N = B%fastBlockDim
    p%K = p%N

    deallocate(p%CindNelemStart, stat=ist)
    allocate(p%CindNelemStart(4,C%nnzb))
    p%CindNelemStart(:,:) = 0 ! init
    
    do i01 = 0, 1 ! two times

      nBlockOps = 0
      Start = 0
      nCij = 0

      do iRow = 1, C%mb
        do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow) ; kCol = C%bsrColInd(Cind)
          ! update matrix block element C_full(kCol,iRow)
          
          nCij = nCij + 1
          Nelem = 0
          
          do Bind = B%bsrRowPtr(iRow), B _bsrEndPtr(iRow) ; jCol = B%bsrColInd(Bind)
            ! for each matrix block element B_full(jCol,iRow)

#define     jRow jCol
            ! find out, if A_full(jRow,kCol) exists
            Aind = exists(A, jRow, kCol) ! find out, if A_full(jRow,kCol) exists
            if (Aind > -1) then ! yes

              if(Aind > A%nnzb) stop 'Aind too large' ! DEBUG
              
              ! now: Cval(:,:,Cind) = one*Cval(:,:,Cind) + Aval(:,:,Aind) .times. Aval(:,:,Bind) 
              ! use gemm('n', 'n', M, N, K, 1., A(1:M,1:K,Aind), M, B(1:K,1:N,Bind), K, 1., C(1:M,1:N,Cind), M)

              Nelem = Nelem + 1
              if (i01 > 0) p%AindBind(:,Start + Nelem) = [Aind, Bind]
              nBlockOps = nBlockOps + 1

            endif ! A(jRow,kCol) exists
#undef      jRow

          enddo ! Bind
          
          p%CindNelemStart(:,nCij) = [nCij, Nelem, Start, 0]
          Start = Start + Nelem

        enddo ! Cind
      enddo ! iRow

      if (i01 == 0) then
        p%nCij = nCij
        if (p%nCij > C%nnzb) stop 'too many result blocks found!'
        p%nBlockOps = nBlockOps
        deallocate(p%AindBind, stat=ist)
        allocate(p%AindBind(2,p%nBlockOps))
        p%AindBind(:,:) = 0 ! init
      else
        if (p%nBlockOps /= nBlockOps) stop 'fatal counting error!'
      endif

    enddo ! i01
    if (present(GiFlop)) GiFlop = nBlockOps*(p%M*8.*p%N*.5d0**30*p%K) ! assume a complex data_t, so each FMA has 8 Flop
  endsubroutine ! create

  
  elemental subroutine destroyMultBRSplan(p)
    type(MultBSRplan), intent(inout) :: p !> plan

    integer :: ist
    deallocate(p%CindNelemStart, p%AindBind, stat=ist)
  endsubroutine ! destroy
  
  
  
  
  
  subroutine createBSR_from_VBSR(self, lmax, mb, values, ia, ja, bsrVal)
    !! get the data from a variable block row format
    type(BlockSparseRow), intent(inout) :: self
    integer, intent(in) :: lmax, mb
    complex_data_t, intent(in) :: values(:) !< dim(N*N*nnzb)
    integer, intent(in) :: ia(:) !< dim(mb+1)
    integer, intent(in) :: ja(:) !< dim(nnzb)
    complex_data_t, allocatable, intent(inout), optional :: bsrVal(:,:,:) !< dim(N,N,nnz) ! warning: allocation inside this routine
    
    integer :: ist
    
    self%fastBlockDim = (lmax + 1)**2
    self%slowBlockDim = (lmax + 1)**2
    self%mb = mb
    self%nnzb = size(ja)
    
    if (present(bsrVal)) then
      deallocate(bsrVal, stat=ist)
      allocate(bsrVal(self%fastBlockDim,self%fastBlockDim,self%nnzb), stat=ist)
      ! copy data in
      bsrVal = reshape(values, [self%fastBlockDim,self%fastBlockDim,self%nnzb])
    endif

    allocate(self%bsrColInd(self%nnzb), self%bsrRowPtr(mb+1), stat=ist)
    
    ! copy index lists
    self%bsrRowPtr(:) = ia(:)
    self%bsrColInd(:) = ja(:)

    self%nb = maxval(self%bsrColInd)
    
#ifdef BSRX
    allocate(self%bsrEndPtr(self%mb), stat=ist)
    self%bsrEndPtr(:) = self%bsrRowPtr(2:) - 1
#endif
  endsubroutine ! create

  
  subroutine createBSR_from_full(self, values, bsrVal)
    !! get the data from a dense block-matrix
    type(BlockSparseRow), intent(inout) :: self
    complex_data_t, intent(in) :: values(:,:,:,:) !< dim(N,nb,M,mb)
    complex_data_t, allocatable, intent(inout), optional :: bsrVal(:,:,:) !< dim(N,M,nnz) ! warning: allocation inside this routine
    
    integer :: ist, iRow, Bind, jCol
    logical(kind=1), allocatable :: nz(:,:)
    
    self%fastBlockDim = size(values, 1)
    self%nb           = size(values, 2) !! number of columns self%nb is not in use
    self%slowBlockDim = size(values, 3)
    self%mb           = size(values, 4) !! number of rows
    allocate(nz(self%nb,self%mb))
    nz(:,:) = any(any(values /= 0.d0, dim=1), dim=1)
    self%nnzb = count(nz)

    if (present(bsrVal)) then
      deallocate(bsrVal, stat=ist)
      allocate(bsrVal(self%fastBlockDim,self%fastBlockDim,self%nnzb), stat=ist)
    endif
    
    allocate(self%bsrColInd(self%nnzb), self%bsrRowPtr(self%mb+1), stat=ist)

    Bind = 0
    self%bsrRowPtr(1) = 1
    do iRow = 1, self%mb
      do jCol = 1, self%nb
        if (nz(jCol,iRow)) then
          Bind = Bind + 1
          ! copy data in
          if (present(bsrVal)) bsrVal(:,:,Bind) = values(:,jCol,:,iRow) ! strided copy
          ! create index list
          self%bsrColInd(Bind) = jCol
          self%bsrRowPtr(iRow+1) = Bind + 1
        endif ! non-zero
      enddo ! jCol     
    enddo ! iRow
    if (Bind /= self%nnzb) stop 'createBSR_from_full: fatal counting error!'

#ifdef BSRX
    allocate(self%bsrEndPtr(self%mb), stat=ist)
    self%bsrEndPtr(:) = self%bsrRowPtr(2:) - 1
#endif
  endsubroutine ! create
  
  
  elemental subroutine destroyBlockSparseRow(self)
    type(BlockSparseRow), intent(inout) :: self
    integer :: ist
    self%fastBlockDim = 0
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
!!!>  ifort -warn -openmp -openmp-report -check all -check-bounds -traceback -O0 -g -mkl -D TESTMAIN_BlockSparseRow IterativeSolver/BlockSparseRow_mod.F90  && ./a.out 133 35
!!!>  gfortran -ffree-line-length-0 -g -D TESTMAIN_BlockSparseRow IterativeSolver/BlockSparseRow_mod.F90 -lblas && ./a.out 1024 123
!!!
program test_bsr
  use BlockSparseRow_mod !, only:
implicit none
  character(len=8) :: CLarg(0:3)
  integer, parameter :: ShowC=0, ShowA=0, ShowB=0, Bfill=16, N = 1 ! BlockSize
  integer :: ilen, ios, iarg, Nrows, Mcols, iRow, jCol, Cind, kCol, ip, nerror(19)=0
  double precision, parameter :: Afill=0.5, point=.5d0**10
  double precision :: col, elem
#define real_data_t double precision
#define gemm DGEMM
  external :: gemm ! BLAS matrix matrix multiplication
  real_data_t, parameter :: one = 1.0, zero = 0.0
  real_data_t, allocatable :: Bfull(:,:), Afull(:,:), Cfull(:,:)
  complex_data_t, allocatable :: Cval(:,:,:), Aval(:,:,:), Bval(:,:,:)
  type(BlockSparseRow) :: A, B!,C==A operators
  real :: GiFlop
  
  type(MultBSRplan) :: plan

  do iarg = 0, ubound(CLarg, 1)
    call get_command_argument(iarg, CLarg(iarg), ilen, ios)
  enddo ! iarg
  read(unit=CLarg(1), fmt=*) Nrows
  read(unit=CLarg(2), fmt=*) Mcols

#define Krows Nrows
  allocate(Cfull(Mcols,Nrows), Afull(Mcols,Krows), Bfull(Krows,Nrows)) ! C:=A*B
  Bfull = zero ; Afull = zero ; Cfull = zero
  
  do iRow = 1, Nrows
    ! fill H
    do while(count(Bfull(:,iRow) /= 0) < min(Bfill, Nrows))
      call random_number(col) ; jCol = ceiling(col*Krows)
!     if(ShowB>0) write(*, fmt="(A,9I6)") 'rand B',jCol,iRow
      Bfull(jCol,iRow) = jCol + point * iRow
    enddo ! while
    if(ShowB>0) write(*, fmt="(A,I4,99F9.4)") 'B',jCol,Bfull(1:min(Krows, 8),jCol)
    ! fill G
    do while(count(Afull(:,iRow) /= 0) < Afill*Mcols)
      call random_number(col) ; jCol = ceiling(col*Mcols)
!     if(ShowA>0) write(*, fmt="(A,9I6)") 'rand G',jCol,iRow
      Afull(jCol,iRow) = jCol + point * iRow
    enddo ! while
    if(ShowA>0) write(*, fmt="(A,I4,99F9.4)") 'A',iRow,Afull(1:min(Mcols, 8),iRow)
  enddo ! iCol

  ! create reference A an m x k matrix, B a k x n matrix and C an m x n matrix
  !                   M      N      K           A             B                   C
  call gemm('n', 'n', Mcols, Nrows, Krows, one, Afull, Mcols, Bfull, Krows, zero, Cfull, Mcols) ! BLAS routine
  !!! do j = 1, n ;  do l = 1, k ;  do i = 1, m ;   c(i,j) = c(i,j) + b(l,j)*a(i,l)   ; enddo ; enddo ; enddo
  ! do j = 1, Nrows ;  do l = 1, Nrows ;  do i = 1, Mcols ; Cfull(i,j) = Cfull(i,j) + Bfull(l,j)*Afull(i,l) ; enddo ; enddo ; enddo

  nerror = 0
  do iRow = 1, Nrows
    if(ShowC>0) write(*, fmt="(A,I4,99F9.4)") 'R',iRow,Cfull(1:min(Mcols, 8),iRow)/Nrows
    do kCol = 1, Mcols
      elem = dot_product(Bfull(:,iRow), Afull(kCol,:)) ! simple evaluation of a single element of a matrix-matrix product
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Cfull(kCol,iRow)) > .1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! kCol
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  call create(A, dcmplx(reshape(Afull, [1,Mcols,1,Krows])), bsrVal=Aval)
  call create(B, dcmplx(reshape(Bfull, [1,Krows,1,Nrows])), bsrVal=Bval)

  !! use the sparse structure of A for C
#define C A
  allocate(Cval(1,1,C%nnzb))
  Cval = 0

  ! API: multiply(A, Aval, G, Aval, Cval, GiFlop)
  call multiply(A, Aval, B, Bval, Cval, GiFlop=GiFlop)
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop'

  nerror = 0
  do iRow = 1, C%mb
    do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow)
      kCol = C%bsrColInd(Cind)
      elem = dreal(Cval(1,1,Cind)) ! take only the real part
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Cfull(kCol,iRow)) > .1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! Cind
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  !============================
  Cval = 0
  
  ! API: multiply(p, A, B, C, GiFlop)
  call create(plan, A, B, A, GiFlop=GiFlop)
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop (planned)'
  ! API: multiply(p, Aval, Bval, Cval)
  call multiply(plan, Aval, Bval, Cval)

  nerror = 0
  do iRow = 1, C%mb
    do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow)
      kCol = C%bsrColInd(Cind)
      elem = dreal(Cval(1,1,Cind)) ! take only the real part
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - Cfull(kCol,iRow)) > .1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    enddo ! Cind
  enddo ! iRow
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(Mcols*.01*Nrows)

  call destroy(plan)
  
  deallocate(Cval)
#undef C
  call destroy(B)
  call destroy(A)
#undef Krows
endprogram ! test

!- TESTMAIN_BlockSparseRow
#endif
