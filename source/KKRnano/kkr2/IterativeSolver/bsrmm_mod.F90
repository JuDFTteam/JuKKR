!> Multiplication of two BSR-matrices (block sparse row matrices)

module bsrmm_mod
  implicit none
  private
  public :: bsr_times_bsr, bsrMultPlan, destroy
  
  type bsrMultPlan
    integer :: nthreads = 1 ! default=1
    integer(kind=8) :: mtasks = 0 !< == maxval(ntasks)
    integer, allocatable :: ntasks(:) !< dim(nthreads)
    integer(kind=4), allocatable :: task(:,:,:) !< dim(0:3,mtasks,nthreads)
    integer(kind=8) :: nFlop = 0, nByte = 0 ! stats
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
  endsubroutine

  subroutine create_bsr_time_bsr_plan(self, ia, ja, shapeA, ix, jx, shapeX)
    type(bsrMultPlan), intent(inout) :: self
    integer, intent(in) :: ia(:) !> dim(A%nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)      ! column indices
    integer, intent(in) :: shapeA(:) ! == [lmsa, lmsd, A%nnzb, ...]
    integer, intent(in) :: ix(:) !> dim(X%nRows + 1) !  start indices
    integer, intent(in) :: jx(:) !> dim(X%nnzb)      ! column indices
    integer, intent(in) :: shapeX(:) ! == [lmsa, nRHS, X%nnzb]
    ! for the result Y we assume the same BSR structure as X
#define iy ix
#define jy jx

    external :: ZGEMM ! from BLAS

    ! local variables
    integer :: leadDim_Y, leadDim_X, leadDim_A, lmsd, nRHS, XnRows
    ! private variables
    integer :: Yind, Aind, Xind
    integer :: iRow, jCol, kCol
    integer :: ist, i01, ibeta
    integer(kind=8) :: nallops, nFlop, nByte
    
    leadDim_A = shapeA(1)
    leadDim_X = shapeX(1)
    leadDim_Y = shapeX(1)
    
    lmsd = shapeA(2)
    nRHS = shapeX(2)

    if (leadDim_A /= leadDim_X) stop __LINE__

    
    call destroy(self) ! deallocate everything
    allocate(self%task(0:3,0:0,0:0))

    XnRows = size(ix) - 1
#define YnRows XnRows
    if (XnRows /= size(ia) - 1) stop __LINE__
    
  do i01 = 0, 1 
    nallops = 0
    nFlop = 0
    nByte = 0
  
    do iRow = 1, YnRows
      ! reuse elements of Y, i.e. keep the accumulator in the cache
      do Yind = iy(iRow), iy(iRow + 1) - 1
#ifdef DEBUG
       if (Yind > size(Y, 3)) stop __LINE__
#endif
      
        jCol = jy(Yind) ! update   matrix block element Y_full[iRow,jCol]

        ! beta = zero
        ibeta = 0
        do Aind = ia(iRow), ia(iRow + 1) - 1
#ifdef DEBUG
          if (Aind > size(A, 3)) stop __LINE__
#endif
          kCol = ja(Aind) ! for each matrix block element A_full[iRow,kCol]
#define kRow kCol          
          Xind = BSR_entry_exists(ix, jx, row=kRow, col=jCol) ! find out if X_full[kRow,jCol] exists
          if (Xind > 0) then ! yes
#ifdef DEBUG
            if (Xind > size(X, 3)) stop __LINE__
#endif
            ! now: Y[:,:,Yind] += A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
            !                    M     N     K          A                       B                             C
!           call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), leadDim_A, X(:,1,Xind), leadDim_X, beta, Y(:,1,Yind), leadDim_Y)

            nallops = nallops + 1
            nFlop = nFlop + (8_8*lmsd)*(nRHS*lmsd)

            self%task(0:3,i01*nallops,i01*1) = [Yind, Aind, Xind, ibeta] ! store task

            ! beta = one
            ibeta = 1
          endif ! X_full[kRow,jCol] exists
#undef kRow
        enddo ! Aind

      enddo ! Yind
    enddo ! iRow
#undef YnRows

    if (i01 == 0) then
      call destroy(self) ! deallocate everything
      ! after 1st iteration
      self%nthreads = 1
      self%mtasks = nallops
      allocate(self%ntasks(self%nthreads), self%task(0:3,self%mtasks,self%nthreads), stat=ist)
      if (ist /= 0) stop __LINE__ ! allocation failed
      self%ntasks(1) = self%mtasks ! one thread takes all tasks
      self%nFlop = nFlop
      self%nByte = nByte
    else
      ! consistency checks in 2nd iteration
      if (self%nFlop /= nFlop)  stop __LINE__
      if (self%mtasks /= nallops) stop __LINE__
    endif

  enddo ! i01

#undef iy
#undef jy

  endsubroutine ! bsr_times_bsr
  
  
  subroutine bsr_times_bsr_planned(self, Y, A, X)
    type(bsrMultPlan), intent(in) :: self
    double complex, intent(out) :: Y(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    double complex, intent(in)  :: A(:,:,:) ! dim(lmsa,lmsd,A%nnzb)
    double complex, intent(in)  :: X(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    ! for the result Y we assume the same BSR structure as X

    ! local variables
    integer :: leadDim_Y, leadDim_X, leadDim_A, lmsd, nRHS
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    ! private variables
    integer :: Yind, Aind, Xind
    integer :: itask, ithread
    double complex :: beta

    leadDim_A = size(A, 1)
    leadDim_X = size(X, 1)
    leadDim_Y = size(Y, 1)
    
    lmsd = size(A, 2)
    nRHS = size(X, 2)

    if (leadDim_Y /= leadDim_X) stop __LINE__
    if (leadDim_A /= leadDim_X) stop __LINE__
    if (size(X, 2) /= size(Y, 2)) stop __LINE__
    if (size(X, 3) /= size(Y, 3)) stop __LINE__

!$omp parallel
!$omp do private(Yind, Aind, Xind, beta, itask, ithread) schedule(static, 1)
    do ithread = 1, self%nthreads
      do itask = 1, self%ntasks(ithread)
      
       Yind = self%task(0,itask,ithread)
       Aind = self%task(1,itask,ithread)
       Xind = self%task(2,itask,ithread)
       beta = self%task(3,itask,ithread) * ONE
#ifdef DEBUG
       if (Yind > size(Y, 3)) stop __LINE__
       if (Aind > size(A, 3)) stop __LINE__
       if (Xind > size(X, 3)) stop __LINE__
#endif
  
        ! now: Y[:,:,Yind] += A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
        !                    M     N     K          A                       B                             C
        call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), leadDim_A, X(:,1,Xind), leadDim_X, beta, Y(:,1,Yind), leadDim_Y)

      enddo ! itask
    enddo ! ithread
!$omp end do
!$omp end parallel

  endsubroutine ! bsr_times_bsr
  
  
  subroutine bsr_times_bsr_spontaneous(Y, ia, ja, A, ix, jx, X, nFlop)
    double complex, intent(out) :: Y(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    integer, intent(in) :: ia(:) !> dim(A%nRows + 1) !  start indices
    integer, intent(in) :: ja(:) !> dim(A%nnzb)      ! column indices
    double complex, intent(in)  :: A(:,:,:) ! dim(lmsa,lmsd,A%nnzb)
    integer, intent(in) :: ix(:) !> dim(X%nRows + 1) !  start indices
    integer, intent(in) :: jx(:) !> dim(X%nnzb)      ! column indices
    double complex, intent(in)  :: X(:,:,:) ! dim(lmsa,nRHS,X%nnzb)
    integer(kind=8), intent(inout) :: nFlop
    ! for the result Y we assume the same BSR structure as X
#define iy ix
#define jy jx

    ! local variables
    integer :: leadDim_Y, leadDim_X, leadDim_A, lmsd, nRHS, XnRows
    double complex, parameter :: ZERO=(0.d0, 0.d0), ONE=(1.d0, 0.d0)
    ! private variables
    integer :: Yind, Aind, Xind
    integer :: iRow, jCol, kCol
    double complex :: beta

    leadDim_A = size(A, 1)
    leadDim_X = size(X, 1)
    leadDim_Y = size(Y, 1)
    
    lmsd = size(A, 2)
    nRHS = size(X, 2)

    if (leadDim_Y /= leadDim_X) stop __LINE__
    if (leadDim_A /= leadDim_X) stop __LINE__
    if (size(X, 2) /= size(Y, 2)) stop __LINE__
    if (size(X, 3) /= size(Y, 3)) stop __LINE__

    XnRows = size(ix) - 1
#define YnRows XnRows
    if (XnRows /= size(ia) - 1) stop __LINE__

!$omp parallel
!$omp do private(iRow, Yind, jCol, Aind, kCol, Xind, beta) reduction(+:nFlop)
    do iRow = 1, YnRows
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
#define kRow kCol          
          Xind = BSR_entry_exists(ix, jx, row=kRow, col=jCol) ! find out if X_full[kRow,jCol] exists
          if (Xind > 0) then ! yes
#ifdef DEBUG
            if (Xind > size(X, 3)) stop __LINE__
#endif
            ! now: Y[:,:,Yind] += A[:,:,Aind] .times. X[:,:,Xind] ! GEMM:  C(m,n) += sum( A(m,:) * B(:,n) )
            !                    M     N     K          A                       B                             C
            call zgemm('n', 'n', lmsd, nRHS, lmsd, one, A(:,1,Aind), leadDim_A, X(:,1,Xind), leadDim_X, beta, Y(:,1,Yind), leadDim_Y)

            nFlop = nFlop + (8_8*lmsd)*(nRHS*lmsd)

            beta = one
          endif ! X_full[kRow,jCol] exists
#undef kRow
        enddo ! Aind

      enddo ! Yind
    enddo ! iRow
#undef YnRows
!$omp end do
!$omp end parallel
#undef iy
#undef jy

  endsubroutine ! bsr_times_bsr

  integer function BSR_entry_exists(RowStart, Colindex, row, col) result(Ind)
    integer, intent(in) :: RowStart(:), ColIndex(:), row, col

    ! ToDo: the ColIndex list should be sorted ascendingly, so bisection search will be faster
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
      write(unit, fmt='(i6,9a)', iostat=ios) gId(iRow),"    ",trim(line)
    enddo ! iRow
  endfunction ! show
#endif
  
endmodule ! vbrmv_mat_mod
