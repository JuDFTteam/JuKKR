!> TFQMR solver


module TFQMR_mod
#include "../DebugHelpers/logging_macros.h"
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  implicit none
  private
  public :: solve

  double complex, parameter, private :: CONE=(1.d0, 0.d0), ZERO=(0.d0, 0.d0)

  interface solve
    module procedure solve_with_TFQMR
  endinterface

#define column_index_t integer  
  
  contains

  !***********************************************************************
  ! v1 = mat_X
  ! v2 = mat_B
  !> @param op             coefficient matrix/operator
  !> @param initial_zero   true - use 0 as initial guess, false: provide own initial guess in mat_X
  !> @param ncol           number of right-hand sides = number of columns of B
  !> @param nrow           number of row elements of matrices mat_X, mat_B
  subroutine solve_with_TFQMR(op, mat_X, mat_B, tolerance, ncol, nRHSs, initial_zero, precond, use_precond, vecs, kernel_timer, &
                   iterations_needed, largest_residual, nFlops) ! optional output args
    USE_LOGGING_MOD
    use TimerMpi_mod, only: TimerMpi
    use SolverStats_mod, only: SolverStats
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, multiply

    type(KKROperator), intent(in) :: op
    double precision, intent(in) :: tolerance
    integer, intent(in) :: ncol, nRHSs
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

    integer, parameter :: MaxIterations = 200 ! limit of max. 2000 iterations
    integer :: iteration, probe, icol, iRHS

    ! small, local arrays with dimension(ncol,nRHSs)
    double complex,   dimension(ncol,nRHSs) :: ZTMP, RHO, ETA, BETA, mALPHA ! -alpha
    double precision, dimension(ncol,nRHSs) :: RUB, DTMP, COSI, TAU, VAR, RESN, R0, N2B ! norm of right-hand side
    integer :: tfqmr_status(ncol,nRHSs) ! 0 = not converged, negative = breakdown, 1 = converged
    integer :: converged_at(ncol,nRHSs) ! stores iteration where calculation converged, 0 = never converged
    logical :: isDone

    !------------ convergence parameters-----------
    double precision :: max_residual, target_upper_bound, max_upper_bound
    double precision, parameter :: TEST_FACTOR = 100.d0, EPSILON_DP = tiny(0.d0)

    !------------- diagnostic variables -------------
    integer :: sparse_mult_count, res_probe_count ! count number of residual calculations
    integer(kind=8) :: mFlops

#define ColIndices op%bsr_X%ColIndex

    mFlops = 0
    tfqmr_status = 0
    converged_at = 0

    target_upper_bound = tolerance * TEST_FACTOR

    isDone = .false.

    sparse_mult_count = 0
    res_probe_count = 0
    
    if (initial_zero) then

      ! set x0 to 0
      mat_X = ZERO
   
      ! v5 = v2 ; r0 = B - A*x0 = B ! no need to multiply A here
      v5 = -mat_B ! subtract RightHandSide

    else

      ! v5 = A*v1
      call apply_precond_and_matrix(op, precond, mat_X, v5, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

      ! v5 = v2 - v5 ; r0 = b - A*x0
      v5 = v5 - mat_B ! subtract RightHandSide

    endif

    v5 = -v5 ! correct for sign change at setup

    ! R0 = norm(v5)
    call col_norms(R0, v5, ColIndices, mFlops)

    ! use norm of B for convergence criterion - use it for residual normalisation instead of B-AX0 contrary to original TFQMR

    ! N2B = norm(v2)
    v4 = ZERO
    v4 = v4 - mat_B ! subtract RightHandSide
    call col_norms(N2B, v4, ColIndices, mFlops) ! col_norms(N2B, mat_B)

    where (abs(N2B) < EPSILON_DP) N2B = 1.d0  ! where N2B = 0 use absolute residual

    ! Supply auxiliary start vector r*
    call ZRANDN(size(v3), v3, 1) ! fill v3 with numbers in [0, 1]

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
        tfqmr_status = -1
      elsewhere
        BETA = ZTMP / RHO
        RHO  = ZTMP
      endwhere

      ! v4 = beta*v4 + v8
      call col_xpay(v8, BETA, v4, ColIndices, mFlops)

      ! v6 = beta*v6 + v5
      call col_xpay(v5, BETA, v6, ColIndices, mFlops)


      ! v9 = A*v6
      call apply_precond_and_matrix(op, precond, v6, v9, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

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
        tfqmr_status = -1
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
        tfqmr_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = -mALPHA * COSI

      ! do not modify brokedown components
      where (tfqmr_status < 0)
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
      call apply_precond_and_matrix(op, precond, v6, v8, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

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
        tfqmr_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = -mALPHA * COSI

      ! do not modify brokedown components
      where (tfqmr_status < 0)
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
      do iRHS = 1, nRHSs
        do icol = 1, ncol
          if (tfqmr_status(icol,iRHS) /= -1) isDone = .false.
        enddo ! icol
      enddo ! iRHS

      if (isDone) exit ! exit iteration loop


      if (MOD(iteration, probe) == 0) then

        res_probe_count = res_probe_count + 1

        ! in case of right-preconditioning
        !                  -1
        !         r = A * M  * y - b      otherwise   r = A * y - b
        !                  2
        ! has to be performed.

        ! v9 = A*v1
        call apply_precond_and_matrix(op, precond, mat_X, v9, vP)!, use_precond, mFlops, kernel_timer, sparse_mult_count)

        ! v9 = v2 - v9
        v9 = v9 - mat_B ! flipped sign ! subtract RightHandSide

        ! RESN = norm(v9)
        call col_norms(RESN, v9, ColIndices, mFlops)

        RESN = RESN / N2B

        isDone = .true.
        do iRHS = 1, nRHSs
          do icol = 1, ncol
            if (RESN(icol,iRHS) > tolerance) then
              if (tfqmr_status(icol,iRHS) == 0) then
                ! if no breakdown has occured continue converging
                isDone = .false.
              endif
            else
              if (tfqmr_status(icol,iRHS) <= 0) then
                tfqmr_status(icol,iRHS) = 1
                converged_at(icol,iRHS) = iteration ! component converged
              endif
            endif
          enddo ! icol
        enddo ! iRHS

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

    do iRHS = 1, nRHSs
      do icol = 1, ncol
        if (converged_at(icol,iRHS) == 0) then
          select case (tfqmr_status(icol,iRHS))
            case (-1)    ; WRITELOG(3,*) "Component not converged (SEVERE breakdown): ", icol
            case (-2)    ; WRITELOG(3,*) "Component not converged (stagnated): ", icol
            case default ; WRITELOG(3,*) "Component not converged: ", icol
          endselect
        endif
      enddo ! icol
    enddo ! iRHS

    WRITELOG(3,*) tfqmr_status
    WRITELOG(3,*) converged_at
    WRITELOG(3,*) RESN

  contains

    !------------------------------------------------------------------------------
    !> Applies the preconditioner (optional), then the sparse matrix on 'mat' and puts result into 'mat_out'.
    !>
    !> mat_out = A P mat_in
    !> preconditioner is used only when use_precond=.true.
    subroutine apply_precond_and_matrix(op, precond, mat_in, mat_out, temp)!, use_precond, mFlops, kernel_timer, sparse_mult_count)
      use KKROperator_mod, only: KKROperator, multiply
      use BCPOperator_mod, only: BCPOperator, multiply
      use TimerMpi_mod, only: startTimer, stopTimer!, TimerMpi

      type(KKROperator), intent(in) :: op
      type(BCPOperator), intent(in) :: precond

      double complex, intent(in)  :: mat_in(:,:,:)
      double complex, intent(out) :: mat_out(:,:,:)
      double complex, intent(out) :: temp(:,:,:)
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

    integer :: col, nrow, block, jCol
    double precision, external :: DZNRM2 ! BLAS double complex 2-norm

    nrow = size(vector, 1)
#ifndef NDEBUG
    if (size(norms, 1) /= size(vector, 2)) stop 'col_norms: strange!'
    if (size(colInd) /= size(vector, 3)) stop 'col_norms: strange! (colInd vs. v)'
#endif    
    norms = 0.d0
    
    !$omp do private(col, block, jCol) reduction(+:norms)
    do block = 1, size(vector, 3)
      jCol = colInd(block)
      do col = 1, size(vector, 2)
        norms(col,jCol) = norms(col,jCol) + DZNRM2(nrow, vector(:,col,block), 1)
      enddo ! col
    enddo ! block
    !$omp end do

    mFlops = mFlops + 4*size(vector)
  endsubroutine ! norms

  !------------------------------------------------------------------------------
  subroutine col_dots(dots, vector, wektor, colInd, mFlops)
    double complex, intent(out) :: dots(:,:)
    double complex, intent(in) :: vector(:,:,:)
    double complex, intent(in) :: wektor(:,:,:)
    column_index_t, intent(in) :: colInd(:)
    integer(kind=8), intent(inout) :: mFlops

    integer :: col, nrow, block, jCol
    double complex, external :: ZDOTU ! BLAS double complex inner product, Unconjugated

    nrow = size(vector, 1)
#ifndef NDEBUG
    if (size(dots, 1) /= size(vector, 2)) stop 'col_dots: strange! (v)'
    if (size(dots, 1) /= size(wektor, 2)) stop 'col_dots: strange! (w)'
    if (any(shape(vector) /= shape(wektor))) stop 'col_dots: strange! (v vs. w)'
    if (size(colInd) /= size(vector, 3)) stop 'col_dots: strange! (colInd vs. v)'
#endif
    dots = ZERO

    !$omp do private(col, block, jCol) reduction(+:dots)
    do block = 1, size(vector, 3)
      jCol = colInd(block)
      do col = 1, size(vector, 2)
        dots(col,jCol) = dots(col,jCol) + ZDOTU(nrow, vector(:,col,block), 1, wektor(:,col,block), 1)
      enddo ! col
    enddo ! block
    !$omp end do

    mFlops = mFlops + 8*size(vector)
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
    if (size(factors, 1) /= size(yvector, 2)) stop 'col_axpy: strange! (y)'
    if (any(shape(xvector) /= shape(yvector))) stop 'col_axpy: strange! (y vs. x)'
    if (size(colInd) /= size(xvector, 3)) stop 'col_axpy: strange! (colInd vs. x)'
#endif

    !$omp do private(col, block, jCol)
    do block = 1, size(yvector, 3)
      jCol = colInd(block)
      do col = 1, size(yvector, 2)
        yvector(:,col,block) = factors(col,jCol) * xvector(:,col,block) + yvector(:,col,block)
      enddo ! col
    enddo ! block
    !$omp end do

    mFlops = mFlops + 8*size(yvector)
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
    if (size(factors, 1) /= size(yvector, 2)) stop 'col_xpay: strange! (y)'
    if (any(shape(xvector) /= shape(yvector))) stop 'col_xpay: strange! (y vs. x)'
    if (size(colInd) /= size(xvector, 3)) stop 'col_xpay: strange! (colInd vs. x)'
#endif

    !$omp do private(col, block, jCol)
    do block = 1, size(yvector, 3)
      jCol = colInd(block)
      do col = 1, size(yvector, 2)
        yvector(:,col,block) = xvector(:,col,block) + factors(col,jCol) * yvector(:,col,block)
      enddo ! col
    enddo ! block
    !$omp end do
    
    mFlops = mFlops + 8*size(yvector)
  endsubroutine ! y := x+a*y

endmodule ! TFQMR_mod
