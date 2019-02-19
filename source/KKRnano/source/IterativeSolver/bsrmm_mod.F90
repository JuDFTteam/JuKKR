!> Multiplication of two BSR-matrices (block sparse row matrices)

!!!!!!! KERNEL == 0: always use zgemm
#define KERNEL 0
#define DEBUG

module bsrmm_mod
  implicit none
  private
  public :: bsr_times_bsr, bsrMultPlan, destroy

  type bsrMultPlan
    integer :: nthreads = 1 ! default=1
    integer(kind=8) :: mtasks = 0 !< == maxval(ntasks)
    integer, allocatable :: ntasks(:) !< dim(nthreads)
    
! ! Structure of Arrays
! #define task_AoSoA(a,b,c) task(b,a,c) 

! Array of Structures
#define task_AoSoA(a,b,c) task(a,b,c) 

    integer(kind=4), allocatable :: task(:,:,:) !< dim(0:3,mtasks,nthreads) ! data layout SoA, ToDo: check if AoS is better
    integer(kind=8) :: nFlop = 0, nByte = 0 ! stats
    integer :: kernel = 0 ! which block multiplication implementation should be used
  endtype

  interface bsr_times_bsr
    module procedure bsr_times_bsr_spontaneous, bsr_times_bsr_planned, create_bsr_time_bsr_plan
  endinterface

  interface destroy
    module procedure destroy_plan
  endinterface

  contains

  elemental subroutine destroy_plan(self)
    type(bsrMultPlan), intent(inout) :: self
    integer :: ist
    deallocate(self%task, self%ntasks, stat=ist) ! ignore status
    self%nthreads = 1
    self%mtasks = 0
    self%nFlop = 0 ! total compute
    self%nByte = 0 ! data transfer
  endsubroutine ! destroy

  subroutine create_bsr_time_bsr_plan(self, ia, ja, shapeA, ix, jx, shapeX)
! #define OVERLAP_PERMUATION
#ifdef  OVERLAP_PERMUATION
    use CacheOverlap_mod, only: maximize_overlap
    integer, allocatable :: perm(:), weights(:)
#else
#define PERM(x) x
#endif

    type(bsrMultPlan), intent(inout) :: self
    integer, intent(in) :: ia(:) !> dim(A%nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)      ! column indices
    integer, intent(in) :: shapeA(:) ! == [lmsa, lmsd, A%nnzb, ...]
    integer, intent(in) :: ix(:) !> dim(X%nRows + 1) !  start indices
    integer, intent(in) :: jx(:) !> dim(X%nnzb)      ! column indices
    integer, intent(in) :: shapeX(:) ! == [lmsa, nRHS, X%nnzb]
    ! for the result Y we assume the same BSR structure and shape as X
#define iy ix
#define jy jx

!$  integer, external :: omp_get_num_threads

    ! local variables
    integer :: leadDim, lmsa, lmsd, nRHS, nRows
    ! private variables
    integer :: Yind, Aind, Xind
    integer :: iRow, jCol, kCol
    integer :: ist, i01, beta, ithread
    integer, parameter :: one=1, zero=0
    character(len=32) :: name
    integer :: iperm

    lmsa    = shapeA(1)
    leadDim = shapeX(1)

    lmsd = shapeA(2)
    
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    nRHS = shapeX(1)
#else
    nRHS = shapeX(2)
    if (lmsa /= leadDim) stop __LINE__
#endif
    
    nRows = size(ix) - 1
    if (nRows /= size(ia) - 1) stop __LINE__

#ifdef  OVERLAP_PERMUATION    
    allocate(perm(nRows), weights(nRows))
    write(*, '(9(a,i0))') "permute order of rows"
!   weights(:) = 1 ! default
    weights(:) = ix(2:nRows+1) - ix(1:nRows) ! how many elements are non-zero in this row of X or Y
    call maximize_overlap(perm, ia, ja, weights)
#endif

  do i01 = 0, 1 ! run two iterations, 1st and 2nd

    if (i01 == 0) then
      ! before 1st iteration
      call destroy(self) ! deallocate everything
      self%nthreads = 1
!$omp parallel
!$      self%nthreads = omp_get_num_threads()
!$omp end parallel
!      if (self%nthreads > 1) write(*, '(9(a,i0))') "create a multiplication plan for  ",self%nthreads," threads"
      allocate(self%ntasks(0:self%nthreads-1)) ! allocate a dummy

    else
      ! after 1st iteration and before 2nd iteration
      self%mtasks = maxval(self%ntasks) ! determine the maximum number of tasks per thread
!      write(*, '(3(a,i0),999(" ",i0))') "multiplication plan for  ",sum(self%ntasks)," tasks on  ",self%nthreads," threads is balanced as  ",self%ntasks

      allocate(self%task_AoSoA(0:3,self%mtasks,0:self%nthreads-1), stat=ist)
      if (ist /= 0) stop __LINE__ ! allocation failed

    endif
  
    self%nFlop = 0
    self%nByte = 0
    self%ntasks(:) = 0 ! init number of tasks per thread

    do iperm = 1, nRows
      iRow = PERM(iperm)
      
      ! reuse elements of Y, i.e. keep the accumulator in the cache
      do Yind = iy(iRow), iy(iRow + 1) - 1
       if (Yind > shapeX(3)) stop __LINE__ ! DEBUG
        ithread = minloc(self%ntasks, dim=1) - 1 ! find out which thread id has so far received the least number of tasks 
        ! it is important that only one threads works on one Yind as different threads may run on different cores or ...
        ! ... even in different NUMA (non-uniform memory access) domains so that the update of Y(:,:,Yind) would lead to false sharing

        jCol = jy(Yind) ! update   matrix block element Y_full[iRow,jCol]

        beta = zero
        do Aind = ia(iRow), ia(iRow + 1) - 1
          if (Aind > shapeA(3)) stop __LINE__ ! DEBUG
          kCol = ja(Aind) ! for each matrix block element A_full[iRow,kCol]
#define   kRow kCol
          Xind = BSR_entry_exists(ix, jx, row=kRow, col=jCol) ! find out if X_full[kRow,jCol] exists
          if (Xind > 0) then ! yes
            if (Xind > shapeX(3)) stop __LINE__ ! DEBUG
            ! now: Y[:,:,Yind] += beta * A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
            !                    M     N     K          A                       B                             C
        !!! call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), lmsa, X(:,1,Xind), leadDim, beta, Y(:,1,Yind), leadDim)

            self%ntasks(ithread) = self%ntasks(ithread) + 1 ! schedule this task on thread i
            if (i01 > 0) self%task_AoSoA(0:3,self%ntasks(ithread),ithread) = [Yind, Aind, Xind, beta] ! store task in 2nd iteration

            self%nFlop = self%nFlop + (8_8*lmsd)*(nRHS*lmsd)

            beta = one
          endif ! X_full[kRow,jCol] exists
#undef    kRow
        enddo ! Aind

      enddo ! Yind
    enddo ! iRow

  enddo ! i01

#undef iy
#undef jy

    !=============================================================================================
    ! decide which kernel routine to be used for block-times-block operations
    !=============================================================================================
    selectcase (lmsa)
    !!! auto suggest hand implemented routines
!   case (    4) ; self%kernel =  4 ; name = "kernel4x4xN"
!   case (   16) ; self%kernel = 16 ; name = "kernel16x16xN"
    case default ; self%kernel =  0 ; name = "zgemm"
    endselect ! lda
#ifdef  KERNEL
!    write(*, '(9(a,i0))') "Warning! kernel is fixed at compile time: use ",KERNEL," while suggested kernel was ",self%kernel," (0=zgemm)"
    self%kernel = KERNEL ! overwrite variable in plan for correct display
#endif
    !=============================================================================================

  endsubroutine ! bsr_times_bsr


  subroutine bsr_times_bsr_planned(self, Y, A, X)
    type(bsrMultPlan), intent(in) :: self
    double complex, intent(out) :: Y(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    double complex, intent(in)  :: A(:,:,:) ! dim(lmsa,lmsd,A%nnzb)
    double complex, intent(in)  :: X(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    ! for the result Y we assume the same BSR structure as X

    external :: zgemm ! from BLAS

    ! local variables
    integer :: leadDim, lmsa, lmsd, nRHS, kernel
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    ! private variables
    integer :: Yind, Aind, Xind
    integer :: itask, ithread
    double complex :: beta

    if (any(shape(X) /= shape(Y))) stop __LINE__

    kernel = self%kernel ! must be lowercase !

    leadDim = size(X, 1)
    lmsa = size(A, 1)
    lmsd = size(A, 2)
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    nRHS = size(X, 1)
    if (nRHS /= leadDim) stop __LINE__
#else
    if (lmsa /= leadDim) stop __LINE__
    nRHS = size(X, 2)
#endif

!!! The thread parallelization forsees schedule(static, 1) as the
!!!       mapping to the individual threads has been done before.
!$omp parallel
!$omp do private(Yind, Aind, Xind, beta, itask, ithread) schedule(static, 1)
    do ithread = 0, self%nthreads - 1
      do itask = 1, self%ntasks(ithread)

        Yind = self%task_AoSoA(0,itask,ithread)
        Aind = self%task_AoSoA(1,itask,ithread)
        Xind = self%task_AoSoA(2,itask,ithread)
#ifdef DEBUG
        if (Yind > size(Y, 3)) stop __LINE__
        if (Aind > size(A, 3)) stop __LINE__
        if (Xind > size(X, 3)) stop __LINE__
#endif

        selectcase (KERNEL) ! this must be in ALLCAPS as we want to replace it by the preprocessor if we compile a production version
        case (4)
            if (self%task_AoSoA(3,itask,ithread) == 0) Y(:,:,Yind) = ZERO
            call kernel4x4xN(nRHS, A(:,:,Aind), X(:,:,Xind), Y(:,:,Yind))
        case (16)
            if (self%task_AoSoA(3,itask,ithread) == 0) Y(:,:,Yind) = ZERO
            call kernel16x16xN(nRHS, A(:,:,Aind), X(:,:,Xind), Y(:,:,Yind))
        case default
            beta = self%task_AoSoA(3,itask,ithread) * ONE

#ifdef  TRANSPOSE_TO_ROW_MAJOR
!           ! now: Y[:,:,Yind] += beta * A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
!           !                    M     N     K          A                     B                        C
            call zgemm('n', 'n', nRHS, lmsd, lmsd, one, X(:,1,Xind), leadDim, A(:,1,Aind), lmsa, beta, Y(:,1,Yind), leadDim)
#else
            ! now: Y[:,:,Yind] += beta * A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
            !                    M     N     K          A                  B                           C
            call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), lmsa, X(:,1,Xind), leadDim, beta, Y(:,1,Yind), leadDim)
#endif

        endselect ! kernel

      enddo ! itask
    enddo ! ithread
!$omp end do
!$omp end parallel

  endsubroutine ! bsr_times_bsr
  
  
  subroutine bsr_times_bsr_spontaneous(Y, ia, ja, A, ix, jx, X, nFlop)
    double complex, intent(out) :: Y(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    integer, intent(in) :: ia(:) !> dim(nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)    ! column indices
    double complex, intent(in)  :: A(:,:,:) ! dim(lmsa,lmsd,A%nnzb)
    integer, intent(in) :: ix(:) !> dim(nRows + 1) !  start indices
    integer, intent(in) :: jx(:) !> dim(X%nnzb)    ! column indices
    double complex, intent(in)  :: X(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    integer(kind=8), intent(inout) :: nFlop
    ! for the result Y we assume the same BSR structure as X
#define iy ix
#define jy jx

    external :: zgemm ! from BLAS

    ! local variables
    integer :: leadDim, lmsa, lmsd, nRHS, nRows
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    ! private variables
    integer :: Yind, Aind, Xind
    integer :: iRow, jCol, kCol
    double complex :: beta

#ifdef  TRANSPOSE_TO_ROW_MAJOR
    stop __FILE__ ! not prepared for this
#endif

    if (any(shape(X) /= shape(Y))) stop __LINE__

    lmsa = size(A, 1)
    leadDim = size(X, 1)

    lmsd = size(A, 2)
    nRHS = size(X, 2)

    if (lmsa /= leadDim) stop __LINE__ ! not prepared for TRANSPOSE_TO_ROW_MAJOR

    nRows = size(ix) - 1
    if (nRows /= size(ia) - 1) stop __LINE__

!$omp parallel
!$omp do private(iRow, Yind, jCol, Aind, kCol, Xind, beta) reduction(+:nFlop)
    do iRow = 1, nRows
      ! reuse elements of Y, i.e. keep the accumulator in the cache
      do Yind = iy(iRow), iy(iRow + 1) - 1
#ifdef DEBUG
       if (Yind > size(Y, 3)) stop __LINE__
#endif

        jCol = jy(Yind) ! update   matrix block element Y_full[iRow,jCol]

        beta = zero
        do Aind = ia(iRow), ia(iRow + 1) - 1
#ifdef DEBUG
          if (Aind > size(A, 3)) stop __LINE__
#endif
          kCol = ja(Aind) ! for each matrix block element A_full[iRow,kCol]
#define   kRow kCol          
          Xind = BSR_entry_exists(ix, jx, row=kRow, col=jCol) ! find out if X_full[kRow,jCol] exists
          if (Xind > 0) then ! yes
#ifdef DEBUG
            if (Xind > size(X, 3)) stop __LINE__
#endif
            ! now: Y[:,:,Yind] += beta * A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
            !                    M     N     K          A                       B                             C
            call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), lmsa, X(:,1,Xind), leadDim, beta, Y(:,1,Yind), leadDim)

            nFlop = nFlop + (8_8*lmsd)*(nRHS*lmsd)

            beta = one
          endif ! X_full[kRow,jCol] exists
#undef    kRow
        enddo ! Aind

      enddo ! Yind
    enddo ! iRow
!$omp end do
!$omp end parallel
#undef iy
#undef jy

  endsubroutine ! bsr_times_bsr

  integer function BSR_entry_exists(RowStart, Colindex, row, col) result(Ind)
    integer, intent(in) :: RowStart(:), ColIndex(:), row, col

    ! ToDo: make sure the ColIndex list is sorted ascendingly, then bisection search will be faster
    do Ind = RowStart(row), RowStart(row + 1) - 1
      if (ColIndex(Ind) == col) return ! Ind
    enddo ! Ind
    Ind = -1 ! not found

  endfunction ! exists

#ifdef DEBUG
  integer function show_BSR_structure(unit, ia, ja, ix, jx) result(ios)
    integer, intent(in) :: unit ! output unit to write to
    integer, intent(in) :: ia(:), ja(:), ix(:), jx(:)
    
    integer :: nRows, iRow, Aind, Xind, jCol
    character(len=192) :: line
    nRows = size(ia) - 1
    do iRow = 1, nRows
      line = "" ! clear line
      do Aind = ia(iRow), ia(iRow + 1) - 1 ; jCol = ja(Aind)                     ; if(jCol > len(line)) cycle
        line(jCol:jCol) = char(48 + mod(Aind, 10))
      enddo ! Aind                                                 ! add space
      do Xind = ix(iRow), ix(iRow + 1) - 1 ; jCol = jx(Xind)       + (nRows + 4) ; if(jCol > len(line)) cycle
        line(jCol:jCol) = char(48 + mod(Xind, 10))
      enddo ! Xind
#define gId
      write(unit, fmt='(i6,9a)', iostat=ios) gId(iRow),"    ",trim(line)
    enddo ! iRow
  endfunction ! show
#endif

  subroutine kernel16x16xN(nRHS, A, X, Y)
#define NM 15
    integer, intent(in)           :: nRHS
    double complex, intent(in)    :: X(0:NM,nRHS), A(0:NM,0:NM)
    double complex, intent(inout) :: Y(0:NM,nRHS)
    
    integer :: i, k
    do i = 1, nRHS
      do k = 0, NM
        Y(0:NM,i) = Y(0:NM,i) + A(0:NM,k)*X(k,i)
      enddo ! k
    enddo ! i
#undef NM
  endsubroutine ! kernel

  subroutine kernel4x4xN(nRHS, A, X, Y)
#define NM 3
    integer, intent(in)           :: nRHS
    double complex, intent(in)    :: X(0:NM,nRHS), A(0:NM,0:NM)
    double complex, intent(inout) :: Y(0:NM,nRHS)
    
    integer :: i, k
    do i = 1, nRHS
      do k = 0, NM
        Y(0:NM,i) = Y(0:NM,i) + A(0:NM,k)*X(k,i)
      enddo ! k
    enddo ! i
#undef NM
  endsubroutine ! kernel

endmodule ! bsrmm_mod

#ifdef TESTMAIN_bsrmm
!+ testmain_bsrmm

!!!
!!!>  ifort -warn -openmp -openmp-report -check all -check-bounds -traceback -O0 -g -mkl -D TESTMAIN_bsrmm bsrmm_mod.F90  && ./a.out 128 32 4
!!!>  ifort -openmp -mkl  -D TESTMAIN_bsrmm bsrmm_mod.F90  && time ./a.out 1024 64 16
!!!>  gfortran -ffree-line-length-0 -g -D TESTMAIN_bsrmm bsrmm_mod.F90 -lblas && ./a.out 1024 123 2
!!!
#define complex_data_t double complex


program test_bsrmm
  use bsrmm_mod !, only:
implicit none

  type BlockSparseRow
    integer :: fastBlockDim = 0 !< == (lmax+1)^2
    integer :: slowBlockDim = 0 !< == (lmax+1)^2
    integer :: mb = 0 !< number of block rows, slow block dim
    integer :: nb = 0 !< number of block cols, fast block dim
    integer :: nnzb = 0 !< number of non-zero blocks
    integer(kind=4), allocatable :: bsrRowPtr(:) !< dim(mb+1)
    integer(kind=4), allocatable :: bsrColInd(:) !< dim(nnzb)
#define  _bsrEndPtr(i)  %bsrRowPtr(i+1)-1
  endtype ! BlockSparseRow

#define create createBSR_from_full

  character(len=8) :: CLarg(0:3), method
  integer, parameter :: ShowR=0, ShowH=0, ShowG=0, Hfill=16
  integer(kind=8) :: nFlop = 0
  integer :: ilen, ios, iarg, mb, nb, kb, M, N, K, nn, mm, kk, bm, bn, bk, Rind, nerror(19)=0, fi, si, bs=2 ! BlockSize
  double precision, parameter :: Gfill=0.5, pointG=.5d0**8, pointH=.5d0**11, Giga = 0.5d0**30
  double precision :: elem, tick, tock, timediff, GiFlop, GiByte
  complex_data_t :: celem
  external :: zgemm ! BLAS matrix matrix multiplication
  complex_data_t, parameter :: one = 1.0, zero = 0.0
  complex_data_t, allocatable :: Hfull(:,:),  Gfull(:,:),  Rfull(:,:),  Rfill(:,:), t4(:,:,:,:)
  complex_data_t, allocatable :: Hval(:,:,:), Gval(:,:,:), Rval(:,:,:)
  logical(kind=1), allocatable :: Hnz(:,:), Gnz(:,:)
  type(BlockSparseRow) :: H, G!,R==G operators
  type(bsrMultPlan) :: plan

#define Wtime omp_get_wtime
  double precision, external :: Wtime

  do iarg = 0, ubound(CLarg, 1)
    call get_command_argument(iarg, CLarg(iarg), ilen, ios)
  enddo ! iarg
  read(unit=CLarg(1), fmt=*) mb ! number of target blocks
  read(unit=CLarg(2), fmt=*) nb ! number of RHS blocks
  read(unit=CLarg(3), fmt=*) bs ! block size

  kb = mb ! H is a square operator

  M = mb*bs
  N = nb*bs
  K = kb*bs

  tick = Wtime() ! start time

  allocate(Hfull(M,K), Hnz(mb,kb), Gfull(K,N), Gnz(kb,nb), Rfull(M,N)) ! H*G=R

  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for allocation  ',timediff,' sec'

  Hfull = zero ; Hnz = .false.
  do bk = 1, size(Hfull, 2)/bs ! fill H (Hamiltonian)
    do while(count(Hnz(:,bk)) < min(Hfill, size(Hfull, 1)/bs))
      call random_number(elem) ; bm = ceiling(size(Hfull, 1)*elem/bs)
!     if(ShowH>0) write(*, fmt="(A,9I6)") 'rand H',bk,bm !!! DEBUG
      do   si = 1, bs ; kk = bs*bk - bs + si
        do fi = 1, bs ; mm = bs*bm - bs + fi
          Hfull(mm,kk) = dcmplx(mm + pointH*kk, mm + pointH*kk)
        enddo ! fi
      enddo ! si
      Hnz(bm,bk) = .true.
    enddo ! while
    if(ShowH>0) write(*, fmt="(A,I4,99F9.4)") 'H',bk,Hfull(1:min(size(Hfull, 1), 8),bk)
  enddo ! bk
  write(*, fmt="(9(A,F0.3))") 'BSRtest: H  ',count(Hnz)*bs*bs*100./size(Hfull)," % =  ",count(Hnz)/1024.," ki of  ",size(Hfull)/(bs*1024.*bs*1024.)," Mi"

  Gfull = zero
  do bn = 1, size(Gfull, 2)/bs ! fill G (Green function)
    do while(count(Gnz(:,bn)) < Gfill*size(Gfull, 1)/bs)
      call random_number(elem) ; bm = ceiling(size(Gfull, 1)*elem/bs)
!     if(ShowG>0) write(*, fmt="(A,9I6)") 'rand G',bn,bm !!! DEBUG
      do   si = 1, bs ; nn = bs*bn - bs + si
        do fi = 1, bs ; mm = bs*bm - bs + fi
          Gfull(mm,nn) = dcmplx(mm + pointG*nn, -nn + pointG*mm)
        enddo ! fi
      enddo ! si
      Gnz(bm,bn) = .true.
    enddo ! while
    if(ShowG>0) write(*, fmt="(A,I4,99F9.4)") 'G',bn,Gfull(1:min(size(Gfull, 1), 8),bn)
  enddo ! nn
  write(*, fmt="(9(A,F0.3))") 'BSRtest: G  ',count(Gnz)*bs*bs*100./size(Gfull)," % =  ",count(Gnz)/1024.," ki of  ",size(Gfull)/(bs*1024.*bs)," ki"

  deallocate(Hnz, Gnz)

  Rfull = zero

  tock = Wtime() ; timediff = tock - tick ; tick = tock  
  write(*, fmt="(9(A,F0.3))") 'time for array filling  ',timediff,' sec'

  GiByte = (M*16d0*K + K*16d0*N + M*16d0*N)*Giga ! Total memory of all three dense matrices
  GiFlop = M*8d0*K*Giga*N
  write(*,"(9(A,F0.6))") "Dense  matrix-matrix  multiply: ",GiFlop,' GiFlop, ',GiByte,' GiByte'
  ! build reference A an M x K matrix, B a K x N matrix and C an M x N matrix using the BLAS routine
  !                    M  N  K       A         B               C         ! GEMM:  C(m,n) += A(m,k) * B(k,n) ! Fortran style
  call zgemm('n', 'n', M, N, K, one, Hfull, M, Gfull, K, zero, Rfull, M) ! here:  R(m,n) += H(m,k) * G(k,n) ! Fortran style
                                                                         ! or:   R[n][m] += G[n][k] * H[k][m]  !  C - style
  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for zgemm  ',timediff,' sec'
  if (timediff > 0.) write(*, fmt="(9(A,F0.3))") 'performance for zgemm    ',GiFlop/timediff,' GiFlop/sec'

#ifdef FULL_DEBUG
!+ full_debug

  nerror = 0
  do nn = 1, size(Rfull, 2)
    if(ShowR>0) write(*, fmt="(A,I4,99F9.4)") 'R',nn,Rfull(1:min(size(Rfull, 1), 8),nn)/K
    do mm = 1, size(Rfull, 1)
      celem = dot_product(Hfull(mm,:), Gfull(:,nn)) ! simple evaluation of a single element of a matrix-matrix product, but VERY SLOW
      call compare(celem, Rfull(mm,nn), nerror)
    enddo ! mm
  enddo ! nn
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(size(Gfull)*.01)

  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for dot_product  ',timediff,' sec'

!- full_debug
#endif 

  allocate(t4(bs,mb,bs,kb))
  t4 = reshape(Hfull, [bs,mb,bs,kb])
  call create(H, t4, bsrVal=Hval) ! reshape to dim(fast,nb=ncols,slow,mb=nrows)
  deallocate(t4, stat=ios)
  write(*,"(A,9('  ',i0))") "BlockSparseRow H: ",[H%fastBlockDim, H%slowBlockDim, H%nnzb]
  allocate(t4(bs,mb,bs,kb))
  t4 = reshape(Gfull, [bs,kb,bs,nb])
  call create(G, t4, bsrVal=Gval) ! reshape to dim(fast,nb=ncols,slow,mb=nrows)
  deallocate(t4, stat=ios)
  write(*,"(A,9('  ',i0))") "BlockSparseRow G: ",[G%fastBlockDim, G%slowBlockDim, G%nnzb]

  !! use the sparse structure of G for R
#define R G
  allocate(Rval(bs,bs,R%nnzb)) ; Rval = 0

  GiByte = 16d0*(size(Hval) + size(Gval) + size(Rval))*Giga ! Total memory of all three sparse matrices (index lists not included)

  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for BSR creation  ',timediff,' sec'

  call bsr_times_bsr(Rval, H%bsrRowPtr, H%bsrColInd, Hval, G%bsrRowPtr, G%bsrColInd, Gval, nFlop)
  GiFlop = nFlop*Giga
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop, ',GiByte,' GiByte'

  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for spontaneous BSR x BSR  ',timediff,' sec'
  if (timediff > 0.) write(*, fmt="(9(A,F0.3))") 'performance for spontan  ',GiFlop/timediff,' GiFlop/sec'

  nerror = 0
  do bm = 1, R%mb
    do Rind = R%bsrRowPtr(bm), R _bsrEndPtr(bm) ; bn = R%bsrColInd(Rind)
      do   si = 1, bs ; nn = bs*bn - bs + si 
        do fi = 1, bs ; mm = bs*bm - bs + fi
          call compare(Rval(fi,si,Rind), Rfull(mm,nn), nerror)
        enddo ! fi
      enddo ! si
    enddo ! Rind
  enddo ! bm
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(size(Gfull)*.01)

  !============================
  Rval = 0

  tick = Wtime() ! start time

  call bsr_times_bsr(plan, H%bsrRowPtr, H%bsrColInd, shape(Hval), G%bsrRowPtr, G%bsrColInd, shape(Gval))
  GiFlop = plan%nFlop*Giga
  write(*,"(9(A,F0.6))") "BlockSparseRow matrix multiply: ",GiFlop,' GiFlop (planned)'

  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for BSR x BSR planning  ',timediff,' sec'

  call bsr_times_bsr(plan, Rval, Hval, Gval)

  tock = Wtime() ; timediff = tock - tick ; tick = tock
  write(*, fmt="(9(A,F0.3))") 'time for BSR x BSR planned  ',timediff,' sec'
  if (timediff > 0.) write(*, fmt="(9(A,F0.3))") 'performance for planned  ',GiFlop/timediff,' GiFlop/sec'

  nerror = 0
  do bm = 1, R%mb
    do Rind = R%bsrRowPtr(bm), R _bsrEndPtr(bm) ; bn = R%bsrColInd(Rind)
      do   si = 1, bs ; nn = bs*bn - bs + si 
        do fi = 1, bs ; mm = bs*bm - bs + fi
          call compare(Rval(fi,si,Rind), Rfull(mm,nn), nerror)
        enddo ! fi
      enddo ! si
    enddo ! Rind
  enddo ! bm
  write(*, fmt="(A,99(' ',F0.1))") " errors", nerror/(size(Gfull)*.01)

  call destroy(plan)

  deallocate(Rval, Hval, Gval, stat=ios)

! allocate(Rfill(M,N))
!   !! test the multiplication with dense matrices  --> ToDo
! deallocate(Rfill)

#undef R
!   call destroy(H)
!   call destroy(G)

  contains

  subroutine compare(elem, eref, nerror)
    complex_data_t, intent(in) :: elem, eref
    integer, intent(inout) :: nerror(:)
    integer :: ip
    do ip = lbound(nerror, 1), ubound(nerror, 1)
      if (abs(elem - eref) > .1d0**ip) nerror(ip) = nerror(ip) + 1
    enddo ! ip
  endsubroutine ! compare

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
    write(*, '(9(A,F0.3))') 'create BSR: filling factor ',self%nnzb/(.01*max(1, size(nz))),' %'

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

! #ifdef BSRX
!     allocate(self%bsrEndPtr(self%mb), stat=ist)
!     self%bsrEndPtr(:) = self%bsrRowPtr(2:) - 1
! #endif
    deallocate(nz, stat=ist)
  endsubroutine ! create

endprogram ! test_bsrmm

!- testmain_bsrmm
#endif