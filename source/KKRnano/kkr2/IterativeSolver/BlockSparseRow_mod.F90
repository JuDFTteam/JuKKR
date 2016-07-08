!! this module is a dummy that tests how the GPU solver has to be brought into KKRnano
module BlockSparseRow_mod
!! format according to http://docs.nvidia.com/cuda/cusparse/#block-compressed-sparse-row-format-bsr
! #include "macros.h"
!   use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private

  public :: BlockSparseRow, create, destroy, multiply
  public :: MultBSRplan

  public :: bsr_shape

! #define BSRX
#define complex_data_t double complex

  type BlockSparseRow
    integer :: fastBlockDim = 0 !< == (lmax+1)^2
    integer :: slowBlockDim = 0 !< == (lmax+1)^2
    integer :: mb = 0 !< number of block rows, slow block dim
    integer :: nb = 0 !< number of block cols, fast block dim
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

!!! !!! WARNING: here, we use column-major data layout within the blocks!
!!! !!! naming in cuSPARSE, 
!!!       blockDim:   fastBlockDim
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
    module procedure multiplyBSR_to_Bmat, multiplyBSR_to_BSR, multiplyBSR_to_BSR_planned
  endinterface

  interface destroy
    module procedure destroyBlockSparseRow, destroyMultBRSplan
  endinterface

  type MultBSRplan
    integer :: N, M, K, nCij, nBlockOps, Annzb, BnnzB, Cnnzb
    integer(kind=4), allocatable :: CindNelemStart(:,:) ! dim(4,nCij)
    integer(kind=4), allocatable :: AindBind(:,:) ! dim(2,nBlockOps)
  endtype

  complex_data_t, parameter, private :: one = (1.d0, 0.d0), zero = (0.d0, 0.d0)

  contains

  subroutine multiplyBSR_to_Bmat(A, Aval, M, Bmat, Cmat)
    !! multiply a blocks sparse operator to a (dense) block vectors
    !! can be used for a single atom per MPI process
    !! could call cuSPARSE GPU kernel
    type(BlockSparseRow), intent(in) :: A
#define K A%slowBlockDim
#define N A%fastBlockDim
    complex_data_t, intent(in)  :: Aval(M,K,*) ! block list of block sparse operator A
    integer, intent(in) :: M
    complex_data_t, intent(in)  :: Bmat(K,N,*) ! block vector input
    complex_data_t, intent(out) :: Cmat(M,N,*) ! block vector result

    ! .. locals ..
    complex_data_t :: beta
    integer :: iRow, jCol, Aind

!   if (dim_error(A, Aval)) stop __LINE__ ! dim_error works only for assumed shape arrays
    
    !$omp parallel do private(iRow, jCol, Aind, beta)
    do iRow = 1, A%mb
      beta = zero ! instead of an initialization of Avec to zero
      do Aind = A%bsrRowPtr(iRow), A _bsrEndPtr(iRow)
        jCol = A%bsrColInd(Aind)

        ! now: Cval[:,:,Cind] += Aval[:,:,Aind] .times. Bval[:,:,Bind] ! GEMM:  C(m,n) += A(m,k)*B(k,n)
        !                    M  N  K       A                  B                        C
        call zgemm('n', 'n', M, N, K, one, Aval(:,1,jCol), M, Bmat(:,1,Aind), K, beta, Cmat(:,1,iRow), M)

        beta = one ! simply add all further contributions in this row
      enddo ! Aind
    enddo ! iRow
    !$omp end parallel do

#undef K
#undef N
  endsubroutine ! multiply

  function bsr_shape(bsr) result(s)
    type(BlockSparseRow), intent(in) :: bsr
    integer :: s(3) ! result
    s = [bsr%fastBlockDim, bsr%slowBlockDim, bsr%nnzb]    
  endfunction ! bsr_shape
    
  logical function dim_error(bsr, val) result(wrong)
    type(BlockSparseRow), intent(in) :: bsr
    complex_data_t, intent(in) :: val(:,:,:)
    
!   wrong = any(shape(val) /= [bsr%fastBlockDim, bsr%slowBlockDim, bsr%nnzb])
    wrong = any(shape(val) /= bsr_shape(bsr))
    
  endfunction ! dim_error
  
  
  integer function exists(bsr, row, col) result(ind)
    type(BlockSparseRow), intent(in) :: bsr
    integer, intent(in) :: row, col
    integer :: i ! local
    
    ind = -1
    ! ToDo: the ColInd list is ordered ascendingly, so bisection search will be faster
    do i = bsr%bsrRowPtr(row), bsr _bsrEndPtr(row)
      if (bsr%bsrColInd(i) /= col) cycle
      ind = i
      return
    enddo ! i
    
  endfunction ! exists


  
  subroutine multiplyBSR_to_BSR(A, Aval, B, Bval, C, Cval, GiFlop) ! unplanned
    type(BlockSparseRow), intent(in) :: A !< Matrix to be inverted
    type(BlockSparseRow), intent(in) :: B !< Green function
    type(BlockSparseRow), intent(in) :: C !< Result operator
    complex_data_t, intent(in)  :: Aval(:,:,:) ! (M,K,*) ! input  block list 
    complex_data_t, intent(in)  :: Bval(:,:,:) ! (K,N,*) ! input  block list
    complex_data_t, intent(out) :: Cval(:,:,:) ! (M,N,*) ! result block list
    real, intent(out), optional :: GiFlop

    ! .. locals ..
    integer :: cRow, aCol, cCol, Bind, Cind, Aind, M, N, K
    integer(kind=8) :: nBlockOps
   
    nBlockOps = 0

    if (dim_error(A, Aval)) stop __LINE__
    if (dim_error(B, Bval)) stop __LINE__
    if (dim_error(C, Cval)) stop __LINE__
 
    M = A%fastBlockDim ; if (M /= C%fastBlockDim) stop __LINE__ ! inner
    K = B%fastBlockDim ; if (K /= A%slowBlockDim) stop __LINE__ ! inner
    N = C%slowBlockDim ; if (N /= B%slowBlockDim) stop __LINE__ ! inner
    
    if (A%nb /= B%mb) then ; write(*,*) "error A%nb /= B%mb", A%nb, B%mb ; stop __LINE__ ; endif ! outer
    if (B%nb /= C%nb) then ; write(*,*) "error B%nb /= C%nb", B%nb, C%nb ; stop __LINE__ ; endif ! outer
    if (C%mb /= A%mb) then ; write(*,*) "error C%mb /= A%mb", C%mb, A%mb ; stop __LINE__ ; endif ! outer
    
    !$omp parallel
    
    !$omp workshare
    Cval(:,:,:) = zero ! initialization to zero
    !$omp end workshare
    
    !$omp do private(cRow, aCol, cCol, Bind, Cind, Aind) reduction(+:nBlockOps)
    do cRow = 1, C%mb
      ! reuse elements of C, i.e. keep the accumulator in the cache
      do Cind = C%bsrRowPtr(cRow), C _bsrEndPtr(cRow) ; cCol = C%bsrColInd(Cind) ! update   matrix block element C_full[cRow,cCol]
        do Aind = A%bsrRowPtr(cRow), A _bsrEndPtr(cRow) ; aCol = A%bsrColInd(Aind) ! for each matrix block element A_full[cRow,aCol]

          Bind = exists(B, row=aCol, col=cCol) ! find out, if B_full[aCol,cCol] exists
          if (Bind > -1) then ! yes

            ! now: Cval[:,:,Cind] += Aval[:,:,Aind] .times. Bval[:,:,Bind] ! GEMM:  C(m,n) += A(m,k)*B(k,n)
            !                    M  N  K       A                  B                       C
            call zgemm('n', 'n', M, N, K, one, Aval(:,1,Aind), M, Bval(:,1,Bind), K, one, Cval(:,1,Cind), M)

            nBlockOps = nBlockOps + 1

          endif ! B_full[aCol,cCol] exists

        enddo ! Aind
      enddo ! Cind
    enddo ! cRow
    !$omp end do

    !$omp end parallel
    
    if (present(GiFlop)) GiFlop = nBlockOps*(M*8.*N*.5d0**30*K) ! assume a complex data_t, so each FMA has 8 Flop
  endsubroutine ! multiply

  
  subroutine multiplyBSR_to_BSR_planned(p, Aval, Bval, Cval)
    type(MultBSRplan), intent(in) :: p !> plan
    complex_data_t, intent(in)  :: Aval(:,:,:) ! (p%M,p%K,*) ! input  block list 
    complex_data_t, intent(in)  :: Bval(:,:,:) ! (p%K,p%N,*) ! input  block list of (square) operator B
    complex_data_t, intent(out) :: Cval(:,:,:) ! (p%M,p%N,*) ! result block list

    ! .. locals ..
    complex_data_t :: beta
    integer :: iCij, start, ielem, Aind, Bind, Cind, Nelem

    ! these asserts work only with deferred shape arrays
    if(any(shape(Aval) /= [p%M,p%K,p%Annzb])) stop __LINE__
    if(any(shape(Bval) /= [p%K,p%N,p%Bnnzb])) stop __LINE__
    if(any(shape(Cval) /= [p%M,p%N,p%Cnnzb])) stop __LINE__
    
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

        ! now: Cval[:,:,Cind] += Aval[:,:,Aind] .times. Bval[:,:,Bind] ! GEMM:  C(m,n) += A(m,k)*B(k,n)
        !                      M    N    K       A                    B                          C
        call zgemm('n', 'n', p%M, p%N, p%K, one, Aval(:,1,Aind), p%M, Bval(:,1,Bind), p%K, beta, Cval(:,1,Cind), p%M)

        beta = one
      enddo ! ielem
    enddo ! iCij
    !$omp end do
    
    !$omp end parallel
    
  endsubroutine ! multiply
  
  
  subroutine create_multBRSplan(p, A, B, C, GiFlop)
    type(MultBSRplan), intent(inout) :: p !> plan
    type(BlockSparseRow), intent(in) :: A !< Matrix to be inverted 
    type(BlockSparseRow), intent(in) :: B !< Green function
    type(BlockSparseRow), intent(in) :: C !< Result operator
    real, intent(out), optional :: GiFlop

    ! .. locals ..
    integer :: iRow, kCol, jCol, Bind, Cind, Aind, i01, ist, Nelem, nCij, Start
    integer(kind=8) :: nBlockOps

    p%M = A%fastBlockDim
    p%N = B%slowBlockDim
    p%K = p%N
    
    p%Annzb = A%nnzb
    p%Bnnzb = B%nnzb
    p%Cnnzb = C%nnzb

    deallocate(p%CindNelemStart, stat=ist)
    allocate(p%CindNelemStart(4,C%nnzb))
    p%CindNelemStart(:,:) = 0 ! init
    
    do i01 = 0, 1 ! two times

      nBlockOps = 0
      Start = 0
      nCij = 0

      do iRow = 1, C%mb
        do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow) ; jCol = C%bsrColInd(Cind)
          ! update matrix block element C_full[iRow,jCol]
          
          nCij = nCij + 1
          Nelem = 0
          
          do Aind = A%bsrRowPtr(iRow), A _bsrEndPtr(iRow) ; kCol = A%bsrColInd(Aind)
            ! for each matrix block element A_full[iRow,kCol]

#define     kRow kCol
            Bind = exists(B, row=kRow, col=jCol) ! find out, if B_full[kRow,jCol] exists
            if (Bind > -1) then ! yes

              if(Bind > B%nnzb) stop 'Bind too large' ! DEBUG
              
              ! now: Cval(:,:,Cind) += Aval(:,:,Aind) .times. Bval(:,:,Bind) 
              ! use gemm('n', 'n', M, N, K, 1., A(m,k,Aind), M, B(k,n,Bind), K, 1., C(m,n,Cind), M)

              Nelem = Nelem + 1
              if (i01 > 0) p%AindBind(:,Start + Nelem) = [Aind, Bind]
              nBlockOps = nBlockOps + 1

            endif ! B_full[kRow,jCol] exists
#undef      kRow

          enddo ! Aind
          
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
      endif

    enddo ! i01
    if (p%nCij      /= nCij     ) stop 'fatal counting error! (nCij)'
    if (p%nBlockOps /= nBlockOps) stop 'fatal counting error! (nBlockOps)'
    if (present(GiFlop)) GiFlop = nBlockOps*(p%M*8.*p%N*.5d0**30*p%K) ! assume a complex data_t, so each FMA has 8 Flop
  endsubroutine ! create

  
  
  
  subroutine fusedMultiplyAdd_BSR(C, Cval, Aval, diag, Bval, GiFlop) ! C = A*diag + B 
    type(BlockSparseRow), intent(in) :: C !< operator structure, same for A and B assumed
    complex_data_t, intent(out) :: Cval(:,:,:) ! result block list dim(fastBlockDim,slowBlockDim,nnzb)
    complex_data_t, intent(in)  :: Aval(:,:,:) ! input  block list dim(fastBlockDim,slowBlockDim,nnzb)
    complex_data_t, intent(in)  :: diag(:,:)   ! diagonal matrix,  dim(slowBlockDim,nb)
    complex_data_t, intent(in), optional :: Bval(:,:,:) ! input    dim(fastBlockDim,slowBlockDim,nnzb)
    real, intent(out), optional :: GiFlop

    ! .. locals ..
!   integer :: iRow
    integer :: jCol, Cind, islow
    complex_data_t :: bvec(C%fastBlockDim)
    
    if (dim_error(C, Aval)) stop __LINE__
    if (present(Bval)) then
    if (dim_error(C, Bval)) stop __LINE__
    endif
    if (dim_error(C, Cval)) stop __LINE__

      ! there could also be other variants: C%slowBlockDim, C%mb, i.e. possible 4 combination
      ! for those that use chack against C%nb and use jCol, we can flatten the loop into:
      !     do Cind = 1, C _bsrEndPtr(C%mb) ; jCol = C%bsrColInd(Cind) ; doIt ; enddo
      ! or  do Cind = 1, C%nnzb ; jCol = C%bsrColInd(Cind) ; doIt ; enddo
      
    !$omp parallel private(bvec)
    
      if (any(shape(diag) /= [C%slowBlockDim, C%nb])) stop __LINE__ ! multiply diag(::) to the fast dimensions
      
      bvec(:) = 0
      !$omp do private(iRow, jCol, Cind, islow, bvec)
      do Cind = 1, C _bsrEndPtr(C%mb) ; jCol = C%bsrColInd(Cind)
#ifdef BSRX
#warning   'Gaps in the BSRX format be ignored!'
#endif
          do islow = 1, C%slowBlockDim
            if (present(Bval)) bvec(:) = Bval(:,islow,Cind)
            Cval(:,islow,Cind) = Aval(:,islow,Cind)*diag(islow,jCol) + bvec(:)
          enddo ! islow

      enddo ! Cind
      !$omp end do

!     !!! multiplying diag from the other side:
!
!     if (any(shape(diag) /= [C%fastBlockDim, C%mb])) stop __LINE__ ! multiply diag(::) to the slow dimensions
!    
!     !$omp do private(iRow, Cind, islow, bvec)
!     do iRow = 1, C%mb
!       bvec(:) = 0
!       do Cind = C%bsrRowPtr(iRow), C _bsrEndPtr(iRow)
! 
!         do islow = 1, C%slowBlockDim
!           if (present(Bval)) bvec(:) = Bval(:,islow,Cind)
!           Cval(:,islow,Cind) = Aval(:,islow,Cind)*diag(:,iRow) + bvec(:)
!         enddo ! islow
! 
!       enddo ! Cind
!     enddo ! iRow
    !$omp end do

    !$omp end parallel

    if (present(GiFlop)) GiFlop = C%fastBlockDim*8.*C%slowBlockDim*.5d0**30*C%nnzb ! assume a complex data_t, so each FMA has 8 Flop
  endsubroutine ! fusedMultiplyAdd
 
  
  
  
  
  
  
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
    complex_data_t, intent(in) :: values(:,:,:,:) !< dim(fast,mb=nrows,slow,nb=ncols)
    complex_data_t, allocatable, intent(inout), optional :: bsrVal(:,:,:) !< dim(fast,slow,nnzb) ! warning: allocation inside this routine

    integer :: ist, iRow, Bind, jCol
    logical(kind=1), allocatable :: nz(:,:)
    
    self%fastBlockDim = size(values, 1)
    self%mb           = size(values, 2) !! number of rows
    self%slowBlockDim = size(values, 3)
    self%nb           = size(values, 4) !! number of columns
    allocate(nz(self%mb,self%nb), stat=ist) ! ToDo: catch status
    nz(:,:) = any(any(values /= 0, dim=1), dim=2) ! inner any: reduce along the innermost dim = fast, outer any: reduce along the middle dim = slow dim
    self%nnzb = count(nz)

    if (present(bsrVal)) then
      deallocate(bsrVal, stat=ist)
      allocate(bsrVal(self%fastBlockDim,self%slowBlockDim,self%nnzb), stat=ist)
    endif

    allocate(self%bsrColInd(self%nnzb), self%bsrRowPtr(self%mb+1), stat=ist)

    Bind = 0
    self%bsrRowPtr(1) = 1
    do iRow = 1, self%mb
      do jCol = 1, self%nb
        if (nz(iRow,jCol)) then
          Bind = Bind + 1
          ! copy data in
          if (present(bsrVal)) bsrVal(:,:,Bind) = values(:,iRow,:,jCol) ! strided copy
          ! create index list
          self%bsrColInd(Bind) = jCol
          self%bsrRowPtr(iRow+1) = Bind + 1
        endif ! non-zero
      enddo ! jCol     
    enddo ! iRow
    if (Bind /= self%nnzb) stop 'createBSR_from_full: fatal counting error!'

!   write(*, '(A,9999(" ",i0))') 'create BSR: RowPtr',self%bsrRowPtr
!   write(*, '(A,9999(" ",i0))') 'create BSR: ColInd',self%bsrColInd
    
#ifdef BSRX
    allocate(self%bsrEndPtr(self%mb), stat=ist)
    self%bsrEndPtr(:) = self%bsrRowPtr(2:) - 1
#endif
    deallocate(nz, stat=ist)
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
  integer, parameter :: ShowR=0, ShowH=0, ShowG=0, Hfill=16, bs=16 ! BlockSize
  integer :: ilen, ios, iarg, mb, nb, kb, M, N, K, ii, jj, jb, ib, Rind, nerror(19)=0, fi, si
  double precision, parameter :: Gfill=0.5, point=.5d0**10
  double precision :: elem, eref
#define real_data_t double precision
#define gemm DGEMM
  external :: gemm ! BLAS matrix matrix multiplication
  real_data_t, parameter :: one = 1.0, zero = 0.0
  real_data_t,    allocatable :: Hfull(:,:),  Gfull(:,:),  Rfull(:,:),  Rfill(:,:)
  complex_data_t, allocatable :: Hval(:,:,:), Gval(:,:,:), Rval(:,:,:)
  type(BlockSparseRow) :: H, G!,R==G operators
  real :: GiFlop
  type(MultBSRplan) :: plan

  do iarg = 0, ubound(CLarg, 1)
    call get_command_argument(iarg, CLarg(iarg), ilen, ios)
  enddo ! iarg
  read(unit=CLarg(1), fmt=*) mb ! number of target blocks
  read(unit=CLarg(2), fmt=*) nb ! number of RHS blocks
  
  kb = mb ! H is a square operator
  
  M = mb*bs
  N = nb*bs
  K = kb*bs
  
  allocate(Hfull(M,K), Gfull(K,N), Rfull(M,N)) ! H*G=R
  Hfull = zero ; Gfull = zero ; Rfull = zero

  do ii = 1, size(Hfull, 2) ! fill H (Hamiltonian)
    do while(count(Hfull(:,ii) /= 0) < min(Hfill, size(Hfull, 1)))
      call random_number(elem) ; jj = ceiling(elem*size(Hfull, 1))
!     if(ShowH>0) write(*, fmt="(A,9I6)") 'rand H',ii,jj
      Hfull(jj,ii) = jj + point * ii
    enddo ! while
    if(ShowH>0) write(*, fmt="(A,I4,99F9.4)") 'H',ii,Hfull(1:min(size(Hfull, 1), 8),ii)
  enddo ! ii

  do ii = 1, size(Gfull, 2) ! fill G (Green function)
    do while(count(Gfull(:,ii) /= 0) < Gfill*size(Gfull, 1))
      call random_number(elem) ; jj = ceiling(elem*size(Gfull, 1))
!     if(ShowG>0) write(*, fmt="(A,9I6)") 'rand G',ii,jj
      Gfull(jj,ii) = jj + point * ii
    enddo ! while
    if(ShowG>0) write(*, fmt="(A,I4,99F9.4)") 'G',ii,Gfull(1:min(size(Gfull, 1), 8),ii)
  enddo ! ii

  ! create reference A an M x K matrix, B a K x N matrix and C an M x N matrix using the BLAS routine
  !                   M  N  K       A         B               C         ! GEMM:  C(m,n) += A(m,k) * B(k,n) ! Fortran style
  call gemm('n', 'n', M, N, K, one, Hfull, M, Gfull, K, zero, Rfull, M) ! here:  R(m,n) += H(m,k) * G(k,n) ! Fortran style
  GiFlop = M*8.*K*.5d0**30*N                                            ! or:   R[n][m] += G[n][k] * H[k][m]  !  C - style
  write(*,"(9(A,F0.6))") "Dense  matrix-matrix  multiply: ",GiFlop,' GiFlop (when complex)'

#ifdef FULL_DEBUG
!+full_debug

  nerror = 0
  do ii = 1, size(Rfull, 2)
    if(ShowR>0) write(*, fmt="(A,I4,99F9.4)") 'R',ii,Rfull(1:min(size(Rfull, 1), 8),ii)/K
    do jj = 1, size(Rfull, 1)
      elem = dot_product(Hfull(jj,:), Gfull(:,ii)) ! simple evaluation of a single element of a matrix-matrix product, but VERY SLOW
      eref = Rfull(jj,ii)
      call compare
    enddo ! jj
  enddo ! ii
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(size(Gfull)*.01)
  
!-full_debug
#endif 

  call create(H, dcmplx(reshape(Hfull, [bs,mb,bs,kb])), bsrVal=Hval) ! reshape to dim(fast,nb=ncols,slow,mb=nrows)
  write(*,"(A,9(' ',i0))") "BlockSparseRow H: ",bsr_shape(H)
  call create(G, dcmplx(reshape(Gfull, [bs,kb,bs,nb])), bsrVal=Gval) ! reshape to dim(fast,nb=ncols,slow,mb=nrows)
  write(*,"(A,9(' ',i0))") "BlockSparseRow G: ",bsr_shape(G)

  !! use the sparse structure of G for R
#define R G
  allocate(Rval(bs,bs,R%nnzb)) ; Rval = 0

  ! API: multiply(A, Aval, B, Bval, C, Cval, GiFlop)
  call multiply(H, Hval, G, Gval, R, Rval, GiFlop=GiFlop)
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop'

  nerror = 0
  do ib = 1, R%mb
    do Rind = R%bsrRowPtr(ib), R _bsrEndPtr(ib) ; jb = R%bsrColInd(Rind)
      do si = 1, bs ; do fi = 1, bs
          elem = dreal(Rval(fi,si,Rind)) ! take only the real part
          eref = Rfull(ib*bs+fi-bs,jb*bs-bs+si)
          call compare
        enddo ; enddo ! fi si
    enddo ! Rind
  enddo ! ib
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(size(Gfull)*.01)
  
  !============================
  Rval = 0
  
  ! API: multiply(p, A, B, C, GiFlop)
  call create(plan, H, G, R, GiFlop=GiFlop)
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop (planned)'
  ! API: multiply(p, Aval, Bval, Cval)
  call multiply(plan, Hval, Gval, Rval)

  nerror = 0
  do ib = 1, R%mb
    do Rind = R%bsrRowPtr(ib), R _bsrEndPtr(ib) ; jb = R%bsrColInd(Rind)
      do si = 1, bs ; do fi = 1, bs
          elem = dreal(Rval(fi,si,Rind)) ! take only the real part
          eref = Rfull(ib*bs-bs+fi,jb*bs-bs+si)
          call compare
      enddo ; enddo ! fi si
    enddo ! Rind
  enddo ! ib
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(size(Gfull)*.01)

  call destroy(plan)
  
  deallocate(Rval)
  

  allocate(Rfill(M,N))
!   !! test the multiplication with dense matrices  --> ToDo
  deallocate(Rfill)
  
#undef R
  call destroy(H)
  call destroy(G)
  
  contains

    subroutine compare()
      integer :: ip
      do ip = 1, ubound(nerror, 1)
        if (abs(elem - eref) > .1d0**ip) nerror(ip) = nerror(ip) + 1  
      enddo ! ip
    endsubroutine
  
endprogram ! test

!- TESTMAIN_BlockSparseRow
#endif
