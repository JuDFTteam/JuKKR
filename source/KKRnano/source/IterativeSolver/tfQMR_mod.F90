!> tfQMR solver

module tfQMR_mod
#include "../DebugHelpers/logging_macros.h"
  use Logging_mod, only:    ! import no name, just mention it for the module dependency 
  implicit none
  private
  public :: solve

! #define EXPORT_tfQMR_PROBLEM

#ifdef  BENCHMARK_tfQMR
  public :: benchmark_tfQMR
#endif

  double complex, parameter, private :: CONE=(1.d0, 0.d0), ZERO=(0.d0, 0.d0)

  interface solve
    module procedure solve_with_tfQMR
  endinterface

#define column_index_t integer  

  contains

  !***********************************************************************
  ! v1 = mat_X
  ! v2 = mat_B
  !> @param op             coefficient matrix/operator
  !> @param initial_zero   true - use 0 as initial guess, false: provide own initial guess in mat_X
  !> @param nCols           number of right-hand sides = number of columns of B
  !> @param nrow           number of row elements of matrices mat_X, mat_B
  subroutine solve_with_tfQMR(op, mat_X, mat_B, tolerance, nRHSs, nCols, initial_zero, precond, use_precond, vecs, kernel_timer, &
                   iterations_needed, largest_residual, nFlops) ! optional output args
    USE_LOGGING_MOD
    use TimerMpi_mod, only: TimerMpi
    use SolverStats_mod, only: SolverStats
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, multiply

    type(KKROperator), intent(in) :: op
    double precision, intent(in) :: tolerance
    integer, intent(in) :: nRHSs ! number of Right Hand Sides
    integer, intent(in) :: nCols ! number of block columns
    double complex, intent(inout) :: mat_X(:,:,:) !< dim(X%fastBlockDim,X%slowBlockDim,X%nnzb)
    double complex, intent(in)    :: mat_B(:,:,:) !< dim(B%fastBlockDim,B%slowBlockDim,B%nnzb)

    logical, intent(in) :: initial_zero
    type(BCPOperator) :: precond
    logical, intent(in) :: use_precond
    double complex, intent(inout) :: vecs(:,:,:,3:) !< workspace dim(X%fastBlockDim,X%slowBlockDim,X%nnzb,3:9+1)
#define v3 vecs(:,:,:,3)
#define v4 vecs(:,:,:,4)
#define v5 vecs(:,:,:,5)
#define v6 vecs(:,:,:,6)
#define v7 vecs(:,:,:,7)
#define v8 vecs(:,:,:,8)
#define v9 vecs(:,:,:,9)
#define vP vecs(:,:,:,10)
    !!  vP is only accessed when preconditioning is active
    type(TimerMpi), intent(inout) :: kernel_timer

    ! optional args
    integer,          intent(out), optional :: iterations_needed
    double precision, intent(out), optional :: largest_residual
    integer(kind=8), intent(inout), optional :: nFlops
    
    ! locals

    integer, parameter :: MaxIterations = 2000 ! limit of max. 2000 iterations
    integer :: iteration, probe, iRHS, icol
    
    ! small, local arrays with dimension(nRHSs,nCols)
    double complex,   dimension(nRHSs,nCols) :: ZTMP, RHO, ETA, BETA, mALPHA ! -alpha
    double precision, dimension(nRHSs,nCols) :: RUB, DTMP, COSI, TAU, VAR, RESN, R0, N2B ! norm of right-hand side
    integer :: tfQMR_status(nRHSs,nCols) ! 0 = not converged, negative = breakdown, 1 = converged
    integer :: converged_at(nRHSs,nCols) ! stores iteration where calculation converged, 0 = never converged
    logical :: isDone

    !------------ convergence parameters-----------
    double precision :: max_residual, target_upper_bound, max_upper_bound
    double precision, parameter :: TEST_FACTOR = 100.d0, EPSILON_DP = tiny(0.d0)

    !------------- diagnostic variables -------------
    integer :: sparse_mult_count, res_probe_count ! count number of residual calculations
    integer(kind=8) :: mFlops
    
    double complex, allocatable :: v3aux(:,:,:) ! ToDo: remove again
    integer :: ibaux

#define ColIndices op%bsr_X%ColIndex
#define MINUS(Y, B) call subset_add(Y, B, -1.d0, op%B_subset_of_X, mFlops)  
#define PLUS(Y,  B) call subset_add(Y, B,  1.d0, op%B_subset_of_X, mFlops)  

    mFlops = 0
    tfQMR_status = 0
    converged_at = 0

    target_upper_bound = tolerance * TEST_FACTOR

    isDone = .false.

    sparse_mult_count = 0
    res_probe_count = 0
    
    if (initial_zero) then

      ! set x0 to 0
      mat_X = ZERO
   
      ! v5 = v2 ; r0 = B - A*x0 = B ! no need to multiply A here
      v5 = ZERO
      PLUS(v5, mat_B) ! v5 = mat_B ! add RightHandSide

    else

      ! v5 = A*v1
      call apply_precond_and_matrix(v5, op, mat_X, precond, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

      ! v5 = v2 - v5 ; r0 = b - A*x0
      MINUS(v5, mat_B) ! v5 = v5 - mat_B ! subtract RightHandSide
      v5 = -v5 ! correct for sign change at setup

    endif

    ! R0 = norm(v5)
    call col_norms(R0, v5, ColIndices, mFlops)

    ! use norm of B for convergence criterion - use it for residual normalisation instead of B-AX0 contrary to original tfQMR

    ! N2B = norm(v2)
    v4 = ZERO
    MINUS(v4, mat_B) ! v4 = v4 - mat_B ! subtract RightHandSide
    call col_norms(N2B, v4, ColIndices, mFlops) ! col_norms(N2B, mat_B)

    where (abs(N2B) < EPSILON_DP) N2B = 1.d0  ! where N2B = 0 use absolute residual

    ! Supply auxiliary start vector r*
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    allocate(v3aux(size(vecs, 2),size(vecs, 1),size(vecs, 3)))
    call ZRANDN(size(v3aux), v3aux, 1) ! fill v3aux with numbers in [0, 1]
    do ibaux = 1, size(vecs, 3)
      vecs(:,:,ibaux,3) = transpose(v3aux(:,:,ibaux))
    enddo ! ibaux
    deallocate(v3aux)
#else
    call ZRANDN(size(v3), v3, 1) ! fill v3 with numbers in [0, 1]
#endif

    ! Initialize the variables.

    RESN = 1.d0
    RHO  = CONE
    VAR  = 0.d0
    ETA  = ZERO
    TAU  = R0 * R0

    v8 = ZERO
    v4 = ZERO
    v6 = ZERO

    probe = 1

    do iteration = 1, MaxIterations

      ! ZTMP = v3*v5
      call col_dots(ZTMP, v3, v5, ColIndices, mFlops)

      where (abs(ZTMP) < EPSILON_DP .or. abs(RHO) < EPSILON_DP)
        ! severe breakdown
        BETA = ZERO
        RHO = ZERO
        tfQMR_status = -1
      elsewhere
        BETA = ZTMP / RHO
        RHO  = ZTMP
      endwhere

      ! v4 = beta*v4 + v8
      call col_xpay(v8, BETA, v4, ColIndices, mFlops)

      ! v6 = beta*v6 + v5
      call col_xpay(v5, BETA, v6, ColIndices, mFlops)


      ! v9 = A*v6
      call apply_precond_and_matrix(v9, op, v6, precond, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

      ! v4 = beta*v4 + v9
      call col_xpay(v9, BETA, v4, ColIndices, mFlops)

      ! ZTMP = v3*v4
      call col_dots(ZTMP, v3, v4, ColIndices, mFlops)

      where (abs(ZTMP) > EPSILON_DP .and. abs(RHO) > EPSILON_DP)
        mALPHA = -RHO / ZTMP
        ZTMP = VAR * ETA / (-mALPHA)
      elsewhere
        ! severe breakdown
        mALPHA = ZERO
        ZTMP = ZERO
        tfQMR_status = -1
      endwhere

      ! v7 = ZTMP*v7 + v6
      call col_xpay(v6, ZTMP, v7, ColIndices, mFlops)

      ! v5 = v5 - alpha*v9
      call col_axpy(mALPHA, v9, v5, ColIndices, mFlops)

      ! DTMP = norm(v5)
      call col_norms(DTMP, v5, ColIndices, mFlops)


      DTMP = DTMP * DTMP
      where (abs(TAU) > EPSILON_DP)
        VAR  = DTMP / TAU
        COSI  = 1.d0 / ( 1.d0 + VAR )
        ZTMP = VAR * COSI
      elsewhere
        ! early convergence or breakdown(stagnation)
        COSI = 0.d0
        VAR = 0.d0
        ZTMP = CONE
        tfQMR_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = -mALPHA * COSI

      ! do not modify brokedown components
      where (tfQMR_status < 0)
        ETA = ZERO
      endwhere

      ! v1 = v1 + eta*v7
      call col_axpy(ETA, v7, mat_X, ColIndices, mFlops)

      ! v6 = v6 - alpha*v4
      call col_axpy(mALPHA, v4, v6, ColIndices, mFlops)

      ! v7 = ZTMP*v7 + v6
      call col_xpay(v6, ZTMP, v7, ColIndices, mFlops)


      !=============================================================
      ! 2nd half-step
      !=============================================================

      ! v8 = A*v6
      call apply_precond_and_matrix(v8, op, v6, precond, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

      ! v5 = v5 - alpha*v8
      call col_axpy(mALPHA, v8, v5, ColIndices, mFlops)

      ! DTMP = norm(v5)
      call col_norms(DTMP, v5, ColIndices, mFlops)

      DTMP = DTMP * DTMP
      where (abs(TAU) > EPSILON_DP)
        VAR  = DTMP / TAU
        COSI  = 1.d0 / ( 1.d0 + VAR )
      elsewhere
        ! early convergence or breakdown
        VAR = 0.d0
        COSI = 0.d0
        tfQMR_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = -mALPHA * COSI

      ! do not modify brokedown components
      where (tfQMR_status < 0)
        ETA = ZERO
      endwhere

      ! v1 = v1 + eta*v7
      call col_axpy(ETA, v7, mat_X, ColIndices, mFlops)

      ! Residual upper bound calculation
      RUB = sqrt( (2*iteration + 1) * TAU) / N2B

      max_upper_bound = maxval(RUB)

      if (max_upper_bound <= target_upper_bound) then
        probe = iteration ! probe residual
      else
        probe = iteration + 1 ! do not probe residual
      endif

      if (iteration == MaxIterations) probe = iteration ! probe residual

      !check for complete breakdown
      isDone = .true.
      do icol = 1, nCols
        do iRHS = 1, nRHSs
          if (tfQMR_status(iRHS,icol) /= -1) isDone = .false.
        enddo ! iRHS
      enddo ! icol

      if (isDone) exit ! exit iteration loop


      if (MOD(iteration, probe) == 0) then

        res_probe_count = res_probe_count + 1

        ! in case of right-preconditioning
        !                  -1
        !         r = A * M  * y - b      otherwise   r = A * y - b
        !                  2
        ! has to be performed.

        ! v9 = A*v1
        call apply_precond_and_matrix(v9, op, mat_X, precond, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

        ! v9 = v2 - v9
        MINUS(v9, mat_B) ! v9 = v9 - mat_B ! flipped sign ! subtract RightHandSide

        ! RESN = norm(v9)
        call col_norms(RESN, v9, ColIndices, mFlops)

        RESN = RESN / N2B

        isDone = .true.
        do icol = 1, nCols
          do iRHS = 1, nRHSs
            if (RESN(iRHS,icol) > tolerance) then
              if (tfQMR_status(iRHS,icol) == 0) then
                ! if no breakdown has occured continue converging
                isDone = .false.
              endif
            else
              if (tfQMR_status(iRHS,icol) <= 0) then
                tfQMR_status(iRHS,icol) = 1
                converged_at(iRHS,icol) = iteration ! component converged
              endif
            endif
          enddo ! iRHS
        enddo ! icol

        max_residual = maxval(RESN)

        target_upper_bound = (max_upper_bound / max_residual) * tolerance

#ifdef DIAGNOSTICMMINVMOD
        write(*,*) iteration, minval(RESN), max_residual, max_upper_bound
#endif

        if (isDone) exit ! exit iteration loop

      endif ! iteration % probe == 0

    enddo ! iteration

    if (use_precond) then
      ! in case of right-preconditioning
      !              -1
      !         x = M  * y
      !              2
      ! has to be performed.
      call multiply(precond, mat_X, vP) ! warning: this part is not counted in mFlops and timers
      mat_X = vP
    endif

    ! ================================================================
    ! solution is in mat_X
    ! ================================================================

    WRITELOG(3,*) "number of sparse matrix multipl.: ", sparse_mult_count
    WRITELOG(3,*) "number of residual probes:        ", res_probe_count
    WRITELOG(3,*) "number of iterations:             ", iteration
    WRITELOG(3,*) "max. residual:             ", max_residual

    if (present(iterations_needed)) iterations_needed = iteration
    if (present(largest_residual)) largest_residual = max_residual
    if (present(nFlops)) nFlops = nFlops + mFlops

    do icol = 1, nCols
      do iRHS = 1, nRHSs
        if (converged_at(iRHS,icol) == 0) then
          select case (tfQMR_status(iRHS,icol))
            case (-1)    ; WRITELOG(3,*) "Component not converged (SEVERE breakdown): ", iRHS
            case (-2)    ; WRITELOG(3,*) "Component not converged (stagnated): ", iRHS
            case default ; WRITELOG(3,*) "Component not converged: ", iRHS
          endselect
        endif
      enddo ! iRHS
    enddo ! icol

    WRITELOG(3,*) tfQMR_status
    WRITELOG(3,*) converged_at
    WRITELOG(3,*) RESN

#ifdef EXPORT_tfQMR_PROBLEM
  call dump_tfQMR_problem(op, tolerance, nRHSs, nCols, mat_X, mat_B)
#endif

  contains

    !------------------------------------------------------------------------------
    !> Applies the preconditioner (optional), then the sparse matrix on 'mat' and puts result into 'mat_out'.
    !>
    !> mat_out = A P mat_in
    !> preconditioner is used only when use_precond=.true.
    subroutine apply_precond_and_matrix(mat_out, op, mat_in, precond, temp)!, use_precond, mFlops, kernel_timer, sparse_mult_count)
      use KKROperator_mod, only: KKROperator, multiply
      use BCPOperator_mod, only: BCPOperator, multiply
      use TimerMpi_mod, only: startTimer, stopTimer!, TimerMpi
      double complex, intent(out)   :: mat_out(:,:,:)
      type(KKROperator), intent(in) :: op
      double complex, intent(in)    :: mat_in(:,:,:)
      type(BCPOperator), intent(in) :: precond
      double complex, intent(out)   :: temp(:,:,:)
!     logical, intent(in) :: use_precond
!     integer(kind=8), intent(inout) :: mFlops
!     type(TimerMpi), intent(inout) :: kernel_timer
!     integer, intent(inout) :: sparse_mult_count

      if (use_precond) then
      
        call multiply(precond, mat_in, temp) ! act with preconditioner on mat_in
        
        call startTimer(kernel_timer)
        call multiply(op, temp, mat_out, mFlops) ! act with op on temp
        call  stopTimer(kernel_timer)
        
      else
      
        call startTimer(kernel_timer)
        call multiply(op, mat_in, mat_out, mFlops) ! act with op on mat_in
        call  stopTimer(kernel_timer)
        
      endif
      
      sparse_mult_count = sparse_mult_count + 1
      
    endsubroutine ! apply

  endsubroutine ! solve
  


  !------------------------------------------------------------------------------
  subroutine col_norms(norms, vector, colInd, mFlops)
    double precision, intent(out) :: norms(:,:)
    double complex, intent(in) :: vector(:,:,:)
    column_index_t, intent(in) :: colInd(:)
    integer(kind=8), intent(inout) :: mFlops

    integer :: col, block, jCol, nrow
    double precision, external :: DZNRM2 ! BLAS double complex 2-norm
    nrow = size(vector, 1)
#ifndef NDEBUG
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    if (size(norms, 1) /= size(vector, 1)) stop 'col_norms: strange! TRANSPOSE_TO_ROW_MAJOR'
#else    
    if (size(norms, 1) /= size(vector, 2)) stop 'col_norms: strange!'
#endif
    if (size(colInd) /= size(vector, 3)) stop 'col_norms: strange! (colInd vs. v)'
#endif
    norms = 0.d0
    
!$omp parallel
    !$omp do private(col, block, jCol) reduction(+:norms)
    do block = 1, size(vector, 3)
      jCol = colInd(block)
      do col = 1, size(vector, 2)
#ifdef  TRANSPOSE_TO_ROW_MAJOR
        norms(:,jCol) = norms(:,jCol) + dreal(vector(:,col,block))**2 + aimag(vector(:,col,block))**2
#else
        norms(col,jCol) = norms(col,jCol) + DZNRM2(nrow, vector(:,col,block), 1)**2 ! this is wrong without the ^2 as DZNRM2 is not a linear function
#endif
      enddo ! col
    enddo ! block
    !$omp end do

    !$omp barrier

    !$omp do private(jCol)
    do jCol = 1, size(norms, 2)
      norms(:,jCol) = sqrt(norms(:,jCol)) ! ToDo: discuss if we need to take the square root at all
    enddo ! jCol
    !$omp end do

!$omp end parallel

    mFlops = mFlops + 4_8*size(vector)
  endsubroutine ! norms

  !------------------------------------------------------------------------------
  subroutine col_dots(dots, vector, wektor, colInd, mFlops)
    double complex, intent(out) :: dots(:,:)
    double complex, intent(in) :: vector(:,:,:)
    double complex, intent(in) :: wektor(:,:,:)
    column_index_t, intent(in) :: colInd(:)
    integer(kind=8), intent(inout) :: mFlops

    integer :: col, block, jCol, nrow
    double complex, external :: ZDOTU ! BLAS double complex inner product, Unconjugated
    nrow = size(vector, 1)
#ifndef NDEBUG
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    if (size(dots, 1) /= size(vector, 1)) stop 'col_dots: strange! TRANSPOSE_TO_ROW_MAJOR'
#else    
    if (size(dots, 1) /= size(vector, 2)) stop 'col_dots: strange! (v)'
#endif
    if (any(shape(vector) /= shape(wektor))) stop 'col_dots: strange! (v vs. w)'
    if (size(colInd) /= size(vector, 3)) stop 'col_dots: strange! (colInd vs. v)'
#endif
    dots = ZERO

!$omp parallel
    !$omp do private(col, block, jCol) reduction(+:dots)
    do block = 1, size(vector, 3)
      jCol = colInd(block)
      do col = 1, size(vector, 2)
#ifdef  TRANSPOSE_TO_ROW_MAJOR
        dots(:,jCol) = dots(:,jCol) + vector(:,col,block)*wektor(:,col,block) ! row-major
#else
        dots(col,jCol) = dots(col,jCol) + ZDOTU(nrow, vector(:,col,block), 1, wektor(:,col,block), 1)
#endif
      enddo ! col
    enddo ! block
    !$omp end do
!$omp end parallel

    mFlops = mFlops + 8_8*size(vector)
  endsubroutine ! dots

  !------------------------------------------------------------------------------
  subroutine col_axpy(factors, xvector, yvector, colInd, mFlops)
    double complex, intent(in) :: factors(:,:)
    double complex, intent(in) :: xvector(:,:,:)
    double complex, intent(inout) :: yvector(:,:,:)
    column_index_t, intent(in) :: colInd(:)
    integer(kind=8), intent(inout) :: mFlops

    integer :: col, block, jCol
#ifndef NDEBUG
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    if (size(factors, 1) /= size(yvector, 1)) stop 'col_axpy: strange! TRANSPOSE_TO_ROW_MAJOR'
#else    
    if (size(factors, 1) /= size(yvector, 2)) stop 'col_axpy: strange! (y)'
#endif
    if (any(shape(xvector) /= shape(yvector))) stop 'col_axpy: strange! (y vs. x)'
    if (size(colInd) /= size(xvector, 3)) stop 'col_axpy: strange! (colInd vs. x)'
#endif

#ifdef  TRANSPOSE_TO_ROW_MAJOR
#define COL :
#endif

!$omp parallel
    !$omp do private(col, block, jCol)
    do block = 1, size(yvector, 3)
      jCol = colInd(block)
      do col = 1, size(yvector, 2)
        yvector(:,col,block) = yvector(:,col,block) + factors(COL,jCol) * xvector(:,col,block)
      enddo ! col
    enddo ! block
    !$omp end do
!$omp end parallel

    mFlops = mFlops + 8_8*size(yvector)
  endsubroutine ! y := a*x+y

  !------------------------------------------------------------------------------
  subroutine col_xpay(xvector, factors, yvector, colInd, mFlops)
    double complex, intent(in) :: factors(:,:)
    double complex, intent(in) :: xvector(:,:,:)
    double complex, intent(inout) :: yvector(:,:,:)
    column_index_t, intent(in) :: colInd(:)
    integer(kind=8), intent(inout) :: mFlops

    integer :: col, block, jCol
#ifndef NDEBUG
#ifdef  TRANSPOSE_TO_ROW_MAJOR
    if (size(factors, 1) /= size(yvector, 1)) stop 'col_xpay: strange! TRANSPOSE_TO_ROW_MAJOR'
#else    
    if (size(factors, 1) /= size(yvector, 2)) stop 'col_xpay: strange! (y)'
#endif
    if (any(shape(xvector) /= shape(yvector))) stop 'col_xpay: strange! (y vs. x)'
    if (size(colInd) /= size(xvector, 3)) stop 'col_xpay: strange! (colInd vs. x)'
#endif

!$omp parallel
    !$omp do private(col, block, jCol)
    do block = 1, size(yvector, 3)
      jCol = colInd(block)
      do col = 1, size(yvector, 2)
        yvector(:,col,block) = xvector(:,col,block) + factors(COL,jCol) * yvector(:,col,block)
      enddo ! col
    enddo ! block
    !$omp end do
!$omp end parallel

#ifdef COL
#undef COL
#endif
    
    mFlops = mFlops + 8_8*size(yvector)
  endsubroutine ! y := x+a*y

  !------------------------------------------------------------------------------
  subroutine subset_add(yvector, bvector, factor, subInd, mFlops)
    double complex, intent(inout) :: yvector(:,:,:)
    double complex, intent(in)    :: bvector(:,:,:)
    double precision, intent(in) :: factor ! can be extended to double complex if that is necessary
    integer, intent(in) :: subInd(:) !> precomputed sub-index list
    integer(kind=8), intent(inout) :: mFlops

    integer :: block, Yind
#ifndef NDEBUG
    if (size(bvector, 1) /= size(yvector, 1)) stop 'subset_add: strange! (1)'
    if (size(bvector, 2) /= size(yvector, 2)) stop 'subset_add: strange! (2)'
    if (size(bvector, 3)  > size(yvector, 3)) stop 'subset_add: strange! (y vs. b)'
    if (size(subInd) /= size(bvector, 3)) stop 'subset_add: strange! (b vs. subInd)'
#endif

!$omp parallel
    !$omp do private(block, Yind)
    do block = 1, size(bvector, 3)
      Yind = subInd(block)
      yvector(:,:,Yind) = yvector(:,:,Yind) + factor * bvector(:,:,block)
    enddo ! block
    !$omp end do
!$omp end parallel
    
    mFlops = mFlops + 4_8*size(bvector)
  endsubroutine ! y := y + a*b with different shapes of Y and B

#ifdef  EXPORT_tfQMR_PROBLEM
  subroutine dump_tfQMR_problem(op, tolerance, nRHSs, nCols, mat_X, mat_B)
    use KKROperator_mod, only: KKROperator
    type(KKROperator), intent(in) :: op
    double precision, intent(in) :: tolerance
    integer, intent(in) :: nRHSs ! number of Right Hand Sides
    integer, intent(in) :: nCols ! number of block columns
    double complex, intent(in)    :: mat_X(:,:,:) !< dim(X%fastBlockDim,X%slowBlockDim,X%nnzb)
    double complex, intent(in)    :: mat_B(:,:,:) !< dim(B%fastBlockDim,B%slowBlockDim,B%nnzb)

    integer, parameter :: iou = 753
    integer, save :: filenumber = 0
    character(len=256) :: filename
    integer :: ios, rank
    include 'mpif.h'
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ios) ! uses the world rank

    filenumber = filenumber + 1
    write(unit=filename, fmt="(9(a,i0))") "tfqmr_problem.",filenumber,"-rank.",rank
    open(iou, file=filename, action="write", status="old", iostat=ios)
    if (ios /= 0) then
      write(*,"(9a)") __FILE__," file '",trim(filename),"' does not exist!" ! use touch to create them
      return ! failed
    else
      write(*,"(9a)") __FILE__," export tfQMR problem to file '",trim(filename),"'."
    endif

    write(iou, "(9(a,i0))")   "nRHSs   ",nRHSs,"   ! number of Right Hand Sides"
    write(iou, "(9(a,i0))")   "nCols   ",nCols,"   ! number of block columns"
    write(iou, "(a,es12.3)")  "tolerance  ",tolerance

    write(iou, "(a)") ! blank line
    write(iou, "(9(a,i0))")       "bsr_A%nCols   ",op%bsr_A%nCols,"   ! number of columns"
    write(iou, "(9(a,i0))")   "size(bsr_A%RowStart) ",size(op%bsr_A%RowStart),"   ! should match nRows+1= ",op%bsr_A%nRows+1 ! should match
    write(iou, "(10(' ',i0))") op%bsr_A%RowStart
    write(iou, "(9(a,i0))")   "size(bsr_A%ColIndex) ",size(op%bsr_A%ColIndex),"   ! should match nnzb= ",op%bsr_A%nnzb ! should match
    write(iou, "(10(' ',i0))") op%bsr_A%ColIndex

    write(iou, "(a)") ! blank line
    write(iou, "(9(a,i0))")       "bsr_B%nCols   ",op%bsr_B%nCols,"   ! number of columns"
    write(iou, "(9(a,i0))")   "size(bsr_B%RowStart)   ",size(op%bsr_B%RowStart),"   ! should match nRows+1= ",op%bsr_B%nRows+1 ! should match
    write(iou, "(10(' ',i0))") op%bsr_B%RowStart
    write(iou, "(9(a,i0))")   "size(bsr_B%ColIndex)   ",size(op%bsr_B%ColIndex),"   ! should match nnzb= ",op%bsr_B%nnzb ! should match
    write(iou, "(10(' ',i0))") op%bsr_B%ColIndex

    write(iou, "(a)") ! blank line
    write(iou, "(9(a,i0))")       "bsr_X%nCols   ",op%bsr_X%nCols,"   ! number of columns"
    write(iou, "(9(a,i0))")   "size(bsr_X%RowStart)   ",size(op%bsr_X%RowStart),"   ! should match nRows+1= ",op%bsr_X%nRows+1 ! should match
    write(iou, "(10(' ',i0))") op%bsr_X%RowStart
    write(iou, "(9(a,i0))")   "size(bsr_X%ColIndex)   ",size(op%bsr_X%ColIndex),"   ! should match nnzb= ",op%bsr_X%nnzb ! should match
    write(iou, "(10(' ',i0))") op%bsr_X%ColIndex

    write(iou, "(a)") ! blank line
    write(iou, "(a,9(' ',i0))") "shape(mat_A)   ",shape(op%mat_A) ! may have a forth entry: 1+Lly
    write(iou, *) op%mat_A

    write(iou, "(a)") ! blank line
    write(iou, "(a,9(' ',i0))") "shape(mat_B)   ",shape(mat_B)
    write(iou, *) mat_B

    write(iou, "(a)") ! blank line
    write(iou, "(a,9(' ',i0))") "shape(mat_X)   ",shape(mat_X)
    write(iou, *) mat_X

    close(iou, iostat=ios)
  endsubroutine ! dump
#endif

#ifdef  BENCHMARK_tfQMR
  subroutine benchmark_tfQMR()
    use TimerMpi_mod, only: TimerMpi, outTimeStats, startTimer, stopTimer, createTimer
    use KKROperator_mod, only: KKROperator, destroy
    use BCPOperator_mod, only: BCPOperator
    use SparseMatrixDescription_mod, only: subset!, SparseMatrixDescription

    type(TimerMpi)    :: kernel_timer, solver_timer
    integer(kind=8)   :: nFlops
    type(KKROperator) :: op
    type(BCPOperator) :: precond
    double precision  :: tolerance, largest_residual
    integer           :: nRHSs, nCols, iterations_needed
    double complex, allocatable :: vecs(:,:,:,:) !< workspace dim(X%fastBlockDim,X%slowBlockDim,X%nnzb,3:9+1)
    integer           :: d(3), ist, iter, nRepetitions = 16
    logical, parameter :: initial_zero=.true., use_precond=.false.

    call load_tfQMR_problem(op, tolerance, nRHSs, nCols, filenumber=0)

    allocate(op%B_subset_of_X(op%bsr_B%nnzb))
    ist = subset(set=op%bsr_X, sub=op%bsr_B, list=op%B_subset_of_X)
    if (ist /= 0) stop 'KKROperator: creation of subset list failed!'

    d = shape(op%mat_X)
    allocate(vecs(d(1),d(2),d(3),3:11)) ! allocate workspace
    vecs(:,:,:,11) = op%mat_X ! copy reference solution from the problem file

    call createTimer(solver_timer); call createTimer(kernel_timer) ! init timers
    write(*, "(9(a,i0))") "tfQMR problem loaded, start ",nRepetitions," repetitions."
    nFlops = 0
    do iter = 1, nRepetitions
      call startTimer(solver_timer)
      !======================================================================
      call solve_with_tfQMR(op, op%mat_X, op%mat_B, tolerance, nRHSs, nCols, &
               initial_zero, precond, use_precond, vecs, kernel_timer, &
               iterations_needed, largest_residual, nFlops)
      !======================================================================
      call stopTimer(solver_timer)
    enddo ! iter
    write(*, "(9(a,es9.1))") "tfQMR quality: largest deviation is",maxval(abs(op%mat_X - vecs(:,:,:,11)))
    write(*, "(9(a,f0.6))") "tfQMR performed ",nFlops*1e-12," TFlop."
    call outTimeStats(solver_timer, name=" tfQMR solver", unit=6)
    call outTimeStats(kernel_timer, name="   KKR kernel", unit=6)

    call destroy(op); deallocate(vecs, stat=ist) ! clean up
  endsubroutine ! benchmark_tfQMR

  subroutine load_tfQMR_problem(op, tolerance, nRHSs, nCols, filenumber)
    use KKROperator_mod, only: KKROperator
    type(KKROperator), intent(inout) :: op
    double precision, intent(out) :: tolerance
    integer, intent(out) :: nRHSs ! number of Right Hand Sides
    integer, intent(out) :: nCols ! number of block columns

    integer, intent(in) :: filenumber
    integer, parameter :: iou = 753
    character(len=256) :: filename, line
    character(len=32) :: word
    integer :: ios, n1, n2, n3

    write(unit=filename, fmt="(9(a,i0))") "tfqmr_problem.",filenumber
    open(iou, file=filename, action="read", iostat=ios)
    if (ios /= 0) then
      write(*,'(9a)') __FILE__," Failed to open ",trim(filename)
      stop
    endif

      read(iou, fmt='(a)', iostat=ios) line ! read 1st line
    do while(ios == 0)
! write(*,"(9a)") trim(filename),": read line ",trim(line) ! DEBUG
      if(line /= "") then

		read(unit=line, fmt=*) word, n1
! write(*,"(9a)") trim(filename),": read word ",trim(word) ! DEBUG
		selectcase(word)
		case("nRHSs")      ; nRHSs = n1
		case("nCols")      ; nCols = n1
		case("tolerance")  ; read(unit=line, fmt=*) word, tolerance

#define ABX bsr_A
		case("bsr_A%nCols"); op%ABX%nCols = n1
		case("size(bsr_A%RowStart)"); op%ABX%nRows = n1-1; allocate(op%ABX%RowStart(n1)); read(iou, fmt=*) op%ABX%RowStart
		case("size(bsr_A%ColIndex)"); op%ABX%nnzb  = n1;   allocate(op%ABX%ColIndex(n1)); read(iou, fmt=*) op%ABX%ColIndex
#define ABX bsr_B
		case("bsr_B%nCols"); op%ABX%nCols = n1
		case("size(bsr_B%RowStart)"); op%ABX%nRows = n1-1; allocate(op%ABX%RowStart(n1)); read(iou, fmt=*) op%ABX%RowStart
		case("size(bsr_B%ColIndex)"); op%ABX%nnzb  = n1;   allocate(op%ABX%ColIndex(n1)); read(iou, fmt=*) op%ABX%ColIndex
#define ABX bsr_X
		case("bsr_X%nCols"); op%ABX%nCols = n1
		case("size(bsr_X%RowStart)"); op%ABX%nRows = n1-1; allocate(op%ABX%RowStart(n1)); read(iou, fmt=*) op%ABX%RowStart
		case("size(bsr_X%ColIndex)"); op%ABX%nnzb  = n1;   allocate(op%ABX%ColIndex(n1)); read(iou, fmt=*) op%ABX%ColIndex
#undef  ABX

		case("shape(mat_A)"); read(unit=line, fmt=*) word, n1, n2, n3; allocate(op%mat_A(n1,n2,n3,0:0)) ; read(iou, fmt=*) op%mat_A
		case("shape(mat_B)"); read(unit=line, fmt=*) word, n1, n2, n3; allocate(op%mat_B(n1,n2,n3))     ; read(iou, fmt=*) op%mat_B
		case("shape(mat_X)"); read(unit=line, fmt=*) word, n1, n2, n3; allocate(op%mat_X(n1,n2,n3))     ; read(iou, fmt=*) op%mat_X

		case default; write(*,*) __FILE__," found unrekognized '",word,"' in ",filename
		endselect ! word
! write(*,"(9a)") trim(filename),": found ",trim(word) ! DEBUG
      else
! write(*,"(9a)") trim(filename),": found empty line" ! DEBUG
      endif ! line not empty

      read(iou, fmt='(a)', iostat=ios) line ! read next line
    enddo ! while

    close(iou, iostat=ios)
  endsubroutine ! load
#endif

endmodule ! tfQMR_mod
