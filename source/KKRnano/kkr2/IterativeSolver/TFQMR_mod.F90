!> TFQMR solver

#include "../DebugHelpers/logging_macros.h"

module TFQMR_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  implicit none
  private
  public :: solve

  double complex, parameter, private :: CONE = (1.d0, 0.d0), ZERO=(0.d0, 0.d0)

  interface solve
    module procedure solve_with_TFQMR
  endinterface
  
  contains

  !***********************************************************************
  ! v1 = mat_X
  ! v2 = mat_B
  !> @param op             coefficient matrix/operator
  !> @param initial_zero   true - use 0 as initial guess, false: provide own initial guess in mat_X
  !> @param ncol           number of right-hand sides = number of columns of B
  !> @param nrow           number of row elements of matrices mat_X, mat_B
  subroutine solve_with_TFQMR(op, mat_X, mat_B, tolerance, ncol, nrow, initial_zero, precond, use_precond, vecs, &
                   iterations_needed, largest_residual) ! optional output args
    USE_LOGGING_MOD
    use SolverStats_mod, only: SolverStats
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, multiply

    type(KKROperator), intent(in) :: op
    double precision, intent(in) :: tolerance
    integer, intent(in) :: nrow, ncol
    double complex, intent(inout) :: mat_X(nrow,ncol)
    double complex, intent(in)    :: mat_B(nrow,ncol)

    logical, intent(in) :: initial_zero
    type(BCPOperator) :: precond
    logical, intent(in) :: use_precond
    double complex, intent(inout) :: vecs(:,:,3:) ! workspace
#define v3 vecs(:,:,3)
#define v4 vecs(:,:,4)
#define v5 vecs(:,:,5)
#define v6 vecs(:,:,6)
#define v7 vecs(:,:,7)
#define v8 vecs(:,:,8)
#define v9 vecs(:,:,9)
#define vP vecs(:,:,10)
    !!  vP is only accessed when preconditioning is active

    ! optional args
    integer,          intent(out), optional :: iterations_needed
    double precision, intent(out), optional :: largest_residual
    
    ! locals

    integer, parameter :: MaxIterations = 2000 ! limit of max. 2000 iterations
    integer :: iteration, probe, icol

    !     small, local arrays with dimension ncol
    double complex,   dimension(ncol) :: ZTMP, RHO, ETA, BETA, mALPHA ! -alpha
    double precision, dimension(ncol) :: residual_upper_bound, DTMP, COSI, TAU, VAR, RESN, R0, N2B  ! norm of right-hand side

    integer :: tfqmr_status(ncol) ! 0 = not converged, negative = breakdown, 1 = converged
    integer :: converged_at(ncol) ! stores iteration where calculation converged, 0 = never converged
    logical :: isDone

    !------------ convergence parameters-----------
    double precision :: max_residual, target_upper_bound, max_upper_bound
    double precision, parameter :: TEST_FACTOR = 100.d0, EPSILON_DP = tiny(0.d0)

    !------------- diagnostic variables -------------
    integer :: sparse_mult_count, res_probe_count ! count number of residual calculations
    
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
      call apply_precond_and_matrix(op, precond, mat_X, v5, vP, use_precond)

      sparse_mult_count = sparse_mult_count + 1

      ! v5 = v2 - v5 ; r0 = b - A*x0
      v5 = v5 - mat_B ! subtract RightHandSide

    endif

    v5 = -v5 ! correct for sign change at setup
    
    ! R0 = norm(v5)
    call col_norms(R0, v5)

    ! use norm of B for convergence criterion - use it for residual normalisation instead of B-AX0 contrary to original TFQMR

    ! N2B = norm(v2)
    v4 = ZERO
    v4 = v4 - mat_B ! subtract RightHandSide
    call col_norms(N2B, v4) ! col_norms(N2B, mat_B)

    where (abs(N2B) < EPSILON_DP) N2B = 1.d0  ! where N2B = 0 use absolute residual

    ! Supply auxiliary start vector r*
    call ZRANDN(nrow*ncol, v3, 1) ! fill v3 with numbers in [0, 1]

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
      call col_dots(ZTMP, v3, v5)

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
      call col_xpay(v8, BETA, v4)

      ! v6 = beta*v6 + v5
      call col_xpay(v5, BETA, v6)


      ! v9 = A*v6
      call apply_precond_and_matrix(op, precond, v6, v9, vP, use_precond)

      sparse_mult_count = sparse_mult_count + 1

      ! v4 = beta*v4 + v9
      call col_xpay(v9, BETA, v4)

      ! ZTMP = v3*v4
      call col_dots(ZTMP, v3, v4)

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
      call col_xpay(v6, ZTMP, v7)

      ! v5 = v5 - alpha*v9
      call col_axpy(mALPHA, v9, v5)

      ! DTMP = norm(v5)
      call col_norms(DTMP, v5)


      DTMP = DTMP * DTMP
      where (abs(TAU) > EPSILON_DP)
        VAR  = DTMP / TAU
        COSI  = 1.d0 / ( 1.d0 + VAR )
        ZTMP = VAR * COSI
      elsewhere
        ! early convergence or breakdown(stagnation)
        COSI = ZERO
        VAR = ZERO
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
      call col_axpy(ETA, v7, mat_X)

      ! v6 = v6 - alpha*v4
      call col_axpy(mALPHA, v4, v6)

      ! v7 = ZTMP*v7 + v6
      call col_xpay(v6, ZTMP, v7)


      !=============================================================
      ! 2nd half-step
      !=============================================================

      ! v8 = A*v6
      call apply_precond_and_matrix(op, precond, v6, v8, vP, use_precond)

      sparse_mult_count = sparse_mult_count + 1

      ! v5 = v5 - alpha*v8
      call col_axpy(mALPHA, v8, v5)

      ! DTMP = norm(v5)
      call col_norms(DTMP, v5)

      DTMP = DTMP * DTMP
      where (abs(TAU) > EPSILON_DP)
        VAR  = DTMP / TAU
        COSI  = 1.d0 / ( 1.d0 + VAR )
      elsewhere
        ! early convergence or breakdown
        VAR = ZERO
        COSI = ZERO
        tfqmr_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = -mALPHA * COSI

      ! do not modify brokedown components
      where (tfqmr_status < 0)
        ETA = ZERO
      endwhere

      ! v1 = v1 + eta*v7
      call col_axpy(ETA, v7, mat_X)

      ! Residual upper bound calculation
      residual_upper_bound = sqrt( (2*iteration + 1) * TAU) / N2B

      max_upper_bound = maxval(residual_upper_bound)

      if (max_upper_bound <= target_upper_bound) then
        probe = iteration ! probe residual
      else
        probe = iteration + 1 ! do not probe residual
      endif

      if (iteration == MaxIterations) probe = iteration ! probe residual

      !check for complete breakdown
      isDone = .true.
      do icol = 1, ncol
        if (tfqmr_status(icol) /= -1) then
          isDone = .false.
        endif
      enddo

      if (isDone) exit ! exit iteration loop


      if (MOD(iteration, probe) == 0) then

        res_probe_count = res_probe_count + 1

        ! in case of right-preconditioning
        !                  -1
        !         r = A * M  * y - b      otherwise   r = A * y - b
        !                  2
        ! has to be performed.

        ! v9 = A*v1
        call apply_precond_and_matrix(op, precond, mat_X, v9, vP, use_precond)

        sparse_mult_count = sparse_mult_count + 1

        ! v9 = v2 - v9
        v9 = v9 - mat_B ! flipped sign ! subtract RightHandSide

        ! RESN = norm(v9)
        call col_norms(RESN, v9)

        RESN = RESN / N2B

        isDone = .true.
        do icol = 1, ncol
          if (RESN(icol) > tolerance) then
            if (tfqmr_status(icol) == 0) then
              ! if no breakdown has occured continue converging
              isDone = .false.
            endif
          else
            if (tfqmr_status(icol) <= 0) then
              tfqmr_status(icol) = 1
              converged_at(icol) = iteration ! component converged
            endif
          endif
        enddo

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
      call multiply(precond, mat_X, vP)
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

    do icol = 1, ncol
      if (converged_at(icol) == 0) then
        select case (tfqmr_status(icol))
          case (-1)    ; WRITELOG(3,*) "Component not converged (SEVERE breakdown): ", icol
          case (-2)    ; WRITELOG(3,*) "Component not converged (stagnated): ", icol
          case default ; WRITELOG(3,*) "Component not converged: ", icol
        endselect
      endif
    enddo ! icol

    WRITELOG(3,*) tfqmr_status
    WRITELOG(3,*) converged_at
    WRITELOG(3,*) RESN

  endsubroutine ! solve

  !------------------------------------------------------------------------------
  !> Applies the preconditioner (optional), then the sparse matrix on 'mat' and puts result
  !> into 'mat_out'.
  !>
  !> mat_out = A P mat_in
  !> preconditioner is used only when use_precond=.true.
  subroutine apply_precond_and_matrix(op, precond, mat_in, mat_out, temp, use_precond)
    use KKROperator_mod, only: KKROperator, multiply
    use BCPOperator_mod, only: BCPOperator, multiply

    type(KKROperator), intent(in) :: op
    type(BCPOperator), intent(in) :: precond

    double complex, intent(in)  :: mat_in(:,:)
    double complex, intent(out) :: mat_out(:,:)
    double complex, intent(out) :: temp(:,:)
    logical, intent(in) :: use_precond

    if (use_precond) then
      call multiply(precond, mat_in, temp)
      call multiply(op, temp, mat_out)
    else
      call multiply(op, mat_in, mat_out)
    endif

  endsubroutine ! apply

  !------------------------------------------------------------------------------
  subroutine col_norms(norms, vectors)
    double precision, intent(out) :: norms(:)
    double complex, intent(in) :: vectors(:,:)

    integer :: col, ncol, nrow
    double precision, external :: DZNRM2 ! BLAS

    ncol = size(norms)
    nrow = size(vectors, 1)

    !$omp do private(col)
    do col = 1, ncol
      norms(col) = DZNRM2(nrow, vectors(:,col), 1)
    enddo ! col
    !$omp end do
    
  endsubroutine ! norms

  !------------------------------------------------------------------------------
  subroutine col_dots(dots, vectorsv, vectorsw)
    double complex, intent(out) :: dots(:)
    double complex, intent(in) :: vectorsv(:,:)
    double complex, intent(in) :: vectorsw(:,:)

    integer :: col, ncol, nrow
    double complex, external :: ZDOTU ! BLAS

    ncol = size(vectorsv, 2)
    nrow = size(vectorsv, 1)

    !$omp do private(col)
    do col = 1, ncol
      dots(col) = ZDOTU(nrow, vectorsv(:,col), 1, vectorsw(:,col), 1)
    enddo ! col
    !$omp end do

  endsubroutine ! dots


  !------------------------------------------------------------------------------
  subroutine col_axpy(factors, xvector, yvector)
    double complex, intent(in) :: factors(:)
    double complex, intent(in) :: xvector(:,:)
    double complex, intent(inout) :: yvector(:,:)

    integer :: col, ncol, nrow

    ncol = size(factors)
    nrow = size(xvector, 1)

    !$omp do private(col)
    do col = 1, ncol
      yvector(:,col) = factors(col) * xvector(:,col) + yvector(:,col)
    enddo ! col
    !$omp end do

  endsubroutine ! y := a*x+y

  !------------------------------------------------------------------------------
  subroutine col_xpay(xvector, factors, yvector)
    double complex, intent(in) :: factors(:)
    double complex, intent(in) :: xvector(:,:)
    double complex, intent(inout) :: yvector(:,:)

    integer :: col, ncol, nrow

    ncol = size(factors)
    nrow = size(xvector, 1)

    !$omp do private(col)
    do col = 1, ncol
      yvector(:,col) = xvector(:,col) + factors(col) * yvector(:,col)
    enddo ! col
    !$omp end do
    
  endsubroutine ! y := x+a*y

endmodule ! TFQMR_mod
