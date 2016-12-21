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

!$omp parallel
!$omp do private(Yind, Aind, Xind, beta, itask, ithread) schedule(static, 1)
    do ithread = 0, self%nthreads - 1
      do itask = 1, self%ntasks(ithread)

        Yind = self%task_AoSoA(0,itask,ithread)
        Aind = self%task_AoSoA(1,itask,ithread)
        Xind = self%task_AoSoA(2,itask,ithread)
        beta = self%task_AoSoA(3,itask,ithread) * ONE
#ifdef DEBUG
        if (Yind > size(Y, 3)) stop __LINE__
        if (Aind > size(A, 3)) stop __LINE__
        if (Xind > size(X, 3)) stop __LINE__
#endif

!         selectcase (KERNEL) ! this must be in ALLCAPS as we want to replace it by the preprocessor if we compile a production version
!         case (4)
!           if (self%task_AoSoA(3,itask,ithread) == 0) Y(:,:,Yind) = ZERO
!           call kernel4x4xN(nRHS, A(:,:,Aind), X(:,:,Xind), Y(:,:,Yind))
!         case (16)
!           if (self%task_AoSoA(3,itask,ithread) == 0) Y(:,:,Yind) = ZERO
!           call kernel16x16xN(nRHS, A(:,:,Aind), X(:,:,Xind), Y(:,:,Yind))
!         case default

#ifdef  TRANSPOSE_TO_ROW_MAJOR
!           ! now: Y[:,:,Yind] += beta * A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
!           !                  M     N     K          A                     B                        C
          call zgemm('n', 'n', nRHS, lmsd, lmsd, one, X(:,1,Xind), leadDim, A(:,1,Aind), lmsa, beta, Y(:,1,Yind), leadDim)
#else          
          ! now: Y[:,:,Yind] += beta * A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
          !                    M     N     K          A                  B                           C
          call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), lmsa, X(:,1,Xind), leadDim, beta, Y(:,1,Yind), leadDim)
#endif

!         endselect ! kernel

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

endmodule ! vbrmv_mat_mod

