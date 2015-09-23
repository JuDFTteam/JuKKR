!> TFQMR solver

#include "../DebugHelpers/logging_macros.h"

#define DOTPRODUCT(VDOTW, V, W) call col_dots(VDOTW, V, W)
#define COLUMNNORMS(NORMS, VECTORS) call col_norms(NORMS, VECTORS)

! YVECTOR = FACTORS*XVECTOR + YVECTOR
#define COLUMN_AXPY(FACTORS, XVECTOR, YVECTOR) call col_axpy(FACTORS, XVECTOR, YVECTOR)
! YVECTOR = -FACTORS*XVECTOR + YVECTOR
#define COLUMN_MAXPY(FACTORS, XVECTOR, YVECTOR) call col_maxpy(FACTORS, XVECTOR, YVECTOR)

! YVECTOR =         XVECTOR + FACTORS*YVECTOR
#define COLUMN_XPAY(FACTORS, XVECTOR, YVECTOR) call col_xpay(FACTORS, XVECTOR, YVECTOR)

module mminvmod_oop_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  implicit none
  private
  public :: mminvmod_oop

  contains

  ! :-)
#define THREE 1
#define FOUR 2
#define FIVE 3
#define SIX 4
#define SEVEN 5
#define EIGHT 6
#define NINE 7

  !***********************************************************************
  ! v1 = mat_X
  ! v2 = mat_B
  !> @param op             coefficient matrix/operator
  !> @param initial_zero   true - use 0 as initial guess, false: provide own initial guess in mat_X
  !> @param num_columns    number of right-hand sides = number of columns of B
  !> @param NLEN           number of row elements of matrices mat_X, mat_B
  subroutine mminvmod_oop(op, mat_X, mat_B, TOL, num_columns, NLEN, initial_zero, stats, precond, use_precond, VECS, temp)
    USE_LOGGING_MOD
    use SolverStats_mod, only: SolverStats
    use OperatorT_mod, only: OperatorT

    class(OperatorT), intent(in) :: op

    double precision, intent(in) :: TOL
    integer, intent(in) :: num_columns
    INTEGER, intent(in) :: NLEN

    double complex, intent(inout) :: mat_X(NLEN,num_columns)
    double complex, intent(in) :: mat_B(NLEN,num_columns)

    logical, intent(in) :: initial_zero
    type(SolverStats), intent(inout) :: stats
    class(OperatorT) :: precond
    logical, intent(in) :: use_precond
    double complex, intent(inout) :: VECS(:,:,:) ! workspace
    double complex, intent(inout) :: temp(:,:) ! workspace only used when use_precond is true

    !----------------- local variables --------------------------------------------

    double complex, parameter :: CONE  = (1.d0, 0.d0)
    double complex, parameter :: CZERO = (0.d0, 0.d0)

    integer :: NLIM

    integer :: IT
    integer :: PROBE
    integer :: ind

    ! local arrays ..

    !     small, local arrays with dimension num_columns
    double complex   :: ZTMP(num_columns)
    double precision :: DTMP(num_columns)

    double complex :: RHO(num_columns)
    double complex :: ETA(num_columns)
    double complex :: BETA(num_columns)
    double complex :: ALPHA(num_columns)

    double precision :: R0(num_columns)
    double precision :: N2B(num_columns)  ! norm of right-hand side
    double precision :: RESN(num_columns)
    double precision :: VAR(num_columns)
    double precision :: TAU(num_columns)
    double precision :: COSI(num_columns)
    double precision :: residual_upper_bound(num_columns)

    ! 0 = not converged, negative = breakdown
    ! 1 = converged
    integer :: tfqmr_status(num_columns)

    ! stores iteration where calculation converged
    ! 0 = never converged
    integer :: converged_at(num_columns)

    logical :: isDone

    !------------ convergence parameters-----------
    double precision :: max_residual
    double precision :: target_upper_bound
    double precision :: max_upper_bound
    double precision, parameter :: TEST_FACTOR = 100.d0
    double precision :: EPSILON_DP

    !------------- diagnostic variables -------------
    integer :: sparse_mult_count
    integer :: res_probe_count     ! count number of residual calculations

    !=======================================================================
    ! INITIALIZATION
    !=======================================================================

    EPSILON_DP = epsilon(0.d0)
    tfqmr_status = 0
    converged_at = 0

    target_upper_bound = TOL * TEST_FACTOR

    isDone = .false.

    NLIM = 2000  ! limit of max. 2000 iterations

    sparse_mult_count = 0
    res_probe_count = 0

    if (initial_zero) then

      ! v5 = v2   r0 = B - AX0 = B
      VECS(:,:,FIVE) = mat_B

      ! set x0 to 0
      mat_X = CZERO

    else

      !==============================================================================
      ! V9 = A*V1
      !==============================================================================

      call apply_precond_and_matrix(op, precond, mat_X, VECS(:,:,NINE), temp, use_precond)

      sparse_mult_count = sparse_mult_count + 1

      !r0 = b - Ax0 = v2 - v9
      VECS(:,:,FIVE) = mat_B - VECS(:,:,NINE)

    endif

    ! R0 = norm(v5)
    COLUMNNORMS(R0, VECS(:,:,FIVE))

    ! use norm of B for convergence criterion - use it for residual normalisation
    ! instead of B-AX0 contrary to original TFQMR

    ! N2B = norm(v2)
    COLUMNNORMS(N2B, mat_B)

    where (abs(N2B) < EPSILON_DP) N2B = 1.d0  ! where N2B = 0 use absolute residual

    ! Supply auxiliary start vector r*
    call ZRANDN (NLEN*num_columns,VECS(:,:,THREE),1)

    !     Initialize the variables.

    RESN = 1.d0
    RHO  = CONE
    VAR  = 0.d0
    ETA  = CZERO
    TAU  = R0 * R0

    VECS(:,:,EIGHT) = CZERO
    VECS(:,:,FOUR) = CZERO
    VECS(:,:,SIX) = CZERO

    PROBE= 1

    !============================================================================
    !============================================================================
    ! ITERATION
    do IT=1, NLIM
      !============================================================================
      !============================================================================

      ! ZTMP = v3*v5
      DOTPRODUCT(ZTMP, VECS(:,:,THREE), VECS(:,:,FIVE))

      where (abs(ZTMP) < EPSILON_DP .or. abs(RHO) < EPSILON_DP)
        ! severe breakdown
        BETA = CZERO
        RHO = CZERO
        tfqmr_status = -1
      elsewhere
        BETA = ZTMP / RHO
        RHO  = ZTMP
      endwhere

      ! v4 = beta*v4 + v8
      COLUMN_XPAY(BETA, VECS(:,:,EIGHT), VECS(:,:,FOUR))

      ! v6 = beta*v6 + v5
      COLUMN_XPAY(BETA, VECS(:,:,FIVE), VECS(:,:,SIX))


      !====================================================================
      ! V9 = A*V6
      !====================================================================

      call apply_precond_and_matrix(op, precond, VECS(:,:,SIX), VECS(:,:,NINE), temp, use_precond)

      sparse_mult_count = sparse_mult_count + 1

      !     VECS(:,:,6) input vector to be multiplied by A = smat
      !     VECS(:,:,9) result


      ! v4 = beta*v4 + v9
      COLUMN_XPAY(BETA, VECS(:,:,NINE), VECS(:,:,FOUR))

      ! ZTMP = v3*v4
      DOTPRODUCT(ZTMP, VECS(:,:,THREE), VECS(:,:,FOUR))

      where (abs(ZTMP) > EPSILON_DP .and. abs(RHO) > EPSILON_DP)
        ALPHA = RHO / ZTMP
        ZTMP = VAR * ETA / ALPHA
      elsewhere
        ! severe breakdown
        ALPHA = CZERO
        ZTMP = CZERO
        tfqmr_status = -1
      endwhere

      ! v7 = ZTMP*v7 + v6
      COLUMN_XPAY(ZTMP, VECS(:,:,SIX), VECS(:,:,SEVEN))

      ! v5 = v5 - alpha*v9
      COLUMN_MAXPY(ALPHA, VECS(:,:,NINE), VECS(:,:,FIVE))

      ! DTMP = norm(v5)
      COLUMNNORMS(DTMP, VECS(:,:,FIVE))


      DTMP = DTMP * DTMP
      where (abs(TAU) > EPSILON_DP)
        VAR  = DTMP / TAU
        COSI  = 1.d0 / ( 1.d0 + VAR )
        ZTMP = VAR * COSI
      elsewhere
        ! early convergence or breakdown(stagnation)
        COSI = CZERO
        VAR = CZERO
        ZTMP = CONE
        tfqmr_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = ALPHA * COSI

      ! do not modify brokedown components
      where (tfqmr_status < 0)
        ETA = CZERO
      endwhere

      ! v1 = v1 + eta*v7
      COLUMN_AXPY(ETA, VECS(:,:,SEVEN), mat_X)

      ! v6 = v6 - alpha*v4
      COLUMN_MAXPY(ALPHA, VECS(:,:,FOUR), VECS(:,:,SIX))

      !ZTMP = VAR*COSI ! moved!

      ! v7 = ZTMP*v7 + v6
      COLUMN_XPAY(ZTMP, VECS(:,:,SIX), VECS(:,:,SEVEN))


      !=============================================================
      ! 2nd half-step
      !=============================================================

      !=========================================
      ! V8 = A*V6
      !=========================================

      call apply_precond_and_matrix(op, precond, VECS(:,:,SIX), VECS(:,:,EIGHT), temp, use_precond)

      sparse_mult_count = sparse_mult_count + 1

      !     VECS(:,:,6) input vector to be multiplied by A = GLLH1
      !     VECS(:,:,8) result


      ! v5 = v5 - alpha*v8
      COLUMN_MAXPY(ALPHA, VECS(:,:,EIGHT), VECS(:,:,FIVE))

      ! DTMP = norm(v5)
      COLUMNNORMS(DTMP, VECS(:,:,FIVE))

      DTMP = DTMP * DTMP
      where (abs(TAU) > EPSILON_DP)
        VAR  = DTMP / TAU
        COSI  = 1.d0 / ( 1.d0 + VAR )
      elsewhere
        ! early convergence or breakdown
        VAR = CZERO
        COSI = CZERO
        tfqmr_status = -2
      endwhere
      TAU  = DTMP * COSI
      ETA  = ALPHA * COSI

      ! do not modify brokedown components
      where (tfqmr_status < 0)
        ETA = CZERO
      endwhere

      ! v1 = v1 + eta*v7
      COLUMN_AXPY(ETA, VECS(:,:,SEVEN), mat_X)

      ! Residual upper bound calculation
      residual_upper_bound = sqrt( (2*IT + 1) * TAU) / N2B

      max_upper_bound = maxval(residual_upper_bound)

      if (max_upper_bound <= target_upper_bound) then
        PROBE = IT ! probe residual
      else
        PROBE = IT+1 ! don't probe residual
      endif

      if (IT == NLIM) then
        PROBE = IT
      endif

      !check for complete breakdown
      isDone = .true.
      do ind=1,num_columns
        if (tfqmr_status(ind) /= -1) then
          isDone = .false.
        endif
      enddo

      if (isDone) exit ! exit iteration loop


      ! >>>>>>>>>>>
      if (MOD(IT,PROBE) == 0) then

        res_probe_count = res_probe_count + 1

        ! in case of right-preconditioning
        !                  -1
        !         r = A * M  * y - b      otherwise   r = A * y - b
        !                  2
        ! >>>>>>>>>>>
        ! has to be performed ..

        !=========================================
        ! V9 = A*V1
        !=========================================

        call apply_precond_and_matrix(op, precond, mat_X, VECS(:,:,NINE), temp, use_precond)

        sparse_mult_count = sparse_mult_count + 1

        !     VECS(:,:,1) input vector to be multiplied by A = GLLH1
        !     VECS(:,:,9) result
        !--------------

        ! v9 = v2 - v9
        VECS(:,:,NINE) = mat_B - VECS(:,:,NINE)

        ! RESN = norm(v9)
        COLUMNNORMS(RESN, VECS(:,:,NINE))

        RESN = RESN / N2B

        isDone = .true.
        do ind=1,num_columns
          if (RESN(ind) > TOL) then
            if (tfqmr_status(ind) == 0) then
              ! if no breakdown has occured continue converging
              isDone = .false.
            endif
          else
            if (tfqmr_status(ind) <= 0) then
              tfqmr_status(ind) = 1
              converged_at(ind) = IT ! component converged
            endif
          endif
        enddo

        max_residual = maxval(RESN)

        target_upper_bound = (max_upper_bound / max_residual) * TOL

#ifdef DIAGNOSTICMMINVMOD
        write(*,*) IT, minval(RESN), max_residual, max_upper_bound
#endif

        if (isDone) exit ! exit iteration loop

      ! <<<<<<<<<<<<
      endif
    ! <<<<<<<<<<<<

    !============================================================================
    !============================================================================
    enddo
   ! ITERATION
   !============================================================================
   !============================================================================
67 continue
      !     Done.

      ! >>>>>>>>>>>

      ! in case of right-preconditioning
      !              -1
      !         x = M  * y
      !              2
      ! has to be performed ..

   if (use_precond) then
     call precond%apply(mat_X, temp)
     mat_X = temp
   endif

    ! ================================================================
    ! solution is in mat_X
    ! ================================================================

   WRITELOG(3,*) "number of sparse matrix multipl.: ", sparse_mult_count
   WRITELOG(3,*) "number of residual probes:        ", res_probe_count
   WRITELOG(3,*) "number of iterations:             ", IT
   WRITELOG(3,*) "max. residual:             ", max_residual

   stats%iterations = IT
   stats%max_residual = max_residual

   do ind=1, num_columns
     if (converged_at(ind) == 0) then
       select case (tfqmr_status(ind))
         case (-1)    ; WRITELOG(3,*) "Component not converged (SEVERE breakdown): ", ind
         case (-2)    ; WRITELOG(3,*) "Component not converged (stagnated): ", ind
         case default ; WRITELOG(3,*) "Component not converged: ", ind
       endselect
     endif
   enddo

   WRITELOG(3,*) tfqmr_status
   WRITELOG(3,*) converged_at
   WRITELOG(3,*) RESN

 endsubroutine

 !------------------------------------------------------------------------------
 !> Applies the preconditioner (optional), then the sparse matrix on 'mat' and puts result
 !> into 'mat_out'.
 !>
 !> mat_out = A P mat_in
 !> preconditioner is used only when use_precond=.true.
 subroutine apply_precond_and_matrix(op, precond, mat, mat_out, temp, use_precond)
   use OperatorT_mod, only: OperatorT
   class(OperatorT), intent(in) :: op, precond
   double complex, intent(in)   :: mat(:,:)
   double complex, intent(out)  :: mat_out(:,:)
   double complex, intent(out)  :: temp(:,:)
   logical, intent(in) :: use_precond

   if (use_precond) then
     call precond%apply(mat, temp)
     call op%apply(temp, mat_out)
   else
     call op%apply(mat, mat_out)
   endif

 endsubroutine

 !------------------------------------------------------------------------------
 subroutine col_AXPY(factors, xvector, yvector)
   double complex, intent(in) :: factors(:)
   double complex, intent(in) :: xvector(:,:)
   double complex, intent(inout) :: yvector(:,:)

   integer :: col, ncol, nlen
   double complex, parameter :: CONE = (1.d0, 0.d0)

   ncol = size(factors)
   nlen = size(xvector,1)

   do col = 1, ncol
     call zaxpby(nlen, yvector(:,col), factors(col), xvector(:,col), CONE, yvector(:,col))
   enddo
 endsubroutine

 !------------------------------------------------------------------------------
 subroutine col_MAXPY(factors, xvector, yvector)
   double complex, intent(in) :: factors(:) 
   double complex, intent(in) :: xvector(:,:)
   double complex, intent(inout) :: yvector(:,:)

   integer :: col, ncol, nlen
   double complex, parameter :: CONE = (1.d0, 0.d0)

   ncol = size(factors)
   nlen = size(xvector,1)

   do col = 1, ncol
     call zaxpby(nlen, yvector(:,col), -factors(col), xvector(:,col), CONE, yvector(:,col))
   enddo ! col
 endsubroutine

 !------------------------------------------------------------------------------
 subroutine col_XPAY(factors, xvector, yvector)
   double complex, intent(in) :: factors(:)
   double complex, intent(in) :: xvector(:,:)
   double complex, intent(inout) :: yvector(:,:)

   integer :: col, ncol, nlen
   double complex, parameter :: CONE = (1.d0, 0.d0)

   ncol = size(factors)
   nlen = size(xvector,1)

   do col = 1, ncol
     call zaxpby(nlen, yvector(:,col), CONE, xvector(:,col), factors(col), yvector(:,col))
   enddo ! col
 endsubroutine

 !------------------------------------------------------------------------------
 subroutine col_norms(norms, vectors)
   double precision, intent(out) :: norms(:)
   double complex, intent(in) :: vectors(:,:)

   integer :: col, ncol, nlen
   double precision, external :: DZNRM2

   ncol = size(norms)
   nlen = size(vectors,1)

   do col = 1, ncol
     norms(col) = DZNRM2(nlen,vectors(:,col),1)
   enddo ! col
 endsubroutine

 !------------------------------------------------------------------------------
 subroutine col_dots(dots, vectorsv, vectorsw)
   double complex, intent(out) :: dots(:)
   double complex, intent(in) :: vectorsv(:,:)
   double complex, intent(in) :: vectorsw(:,:)

   integer :: col, ncol, nlen
   double complex, external :: ZDOTU

   ncol = size(vectorsv,2)
   nlen = size(vectorsv,1)

   do col = 1, ncol
     dots(col) = ZDOTU(nlen, vectorsv(:,col), 1, vectorsw(:,col), 1)
   enddo ! col
 endsubroutine

endmodule
