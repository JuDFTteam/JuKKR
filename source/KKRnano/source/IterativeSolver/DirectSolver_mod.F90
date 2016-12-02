!> Direct solver.
!>

module DirectSolver_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: DirectSolver, solve, destroy !, create
  
  type :: DirectSolver
    double complex, allocatable :: full_A(:,:) !< dim(n,n)
    double complex, allocatable :: full_X(:,:) !< dim(n,nRHSs)
  endtype

  interface solve
    module procedure solve_with_solver
  endinterface

  interface destroy
    module procedure destroy_solver
  endinterface

!   interface create
!     module procedure create_solver
!   endinterface

  contains

!   subroutine create_solver(self)
!     type(DirectSolver), intent(inout) :: self
!   endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Solves problem for right hand side mat_B, solution in mat_X.
  !>
  !> The workspace is allocated on demand and stays allocated.
  !> Deallocate with call IterativeSolver%destroy
  subroutine solve_with_solver(self, op, mat_X)
    use KKROperator_mod, only: KKROperator
    use fillKKRMatrix_mod, only: convertBSRToFullMatrix, convertFullMatrixToBSR
    type(DirectSolver), intent(inout) :: self
    type(KKROperator), intent(in) :: op
    double complex, intent(out) :: mat_X(:,:,:)

    integer :: n, nRHSs, ist, nd
    
    ! adopt the dimensions of mat_X
    n =     op%bsr_A%nRows*size(op%mat_A, 2)
    nRHSs = op%bsr_X%nCols*size(op%mat_X, 2)

    if (any(shape(self%full_A) /= [n,n]) .or. any(shape(self%full_X) /= [n,nRHSs])) then
      ! resize
      call destroy(self)
      nd = n ! memory alignment could be done here, see if convert routines are ready to handle that
      allocate(self%full_A(nd,n), self%full_X(nd,nRHSs), stat=ist)
      if (ist /= 0) die_here("failed to allocate dense matrices with"+(nd*.5**26*(n+nRHSs))+"GiByte!")
    endif ! needs resize

    call convertBSRToFullMatrix(self%full_A, op%bsr_A, op%mat_A(:,:,:,0))
    
!   TESTARRAYLOG(3, full_A)

    call convertBSRToFullMatrix(self%full_X, op%bsr_B, op%mat_B) ! convert op%mat_B to full_B

    ist = solveFull(n, self%full_A, self%full_X) ! on entry, full_X contains mat_B, compute the direct solution using LAPACK
    if (ist /= 0) die_here("failed to directly invert a matrix of dim"+n+"with"+nRHSs+"right hand sides!")

    call convertFullMatrixToBSR(mat_X, op%bsr_X, self%full_X) ! convert back full_X to op%mat_X

  endsubroutine ! solve

  
  !----------------------------------------------------------------------------
  !> Solution of a system of linear equations with multiple right hand sides,
  !> using standard dense matrix LAPACK routines.
  integer function solveFull(n, full_A, full_X) result(info)
    integer, intent(in)           :: n
    double complex, intent(inout) :: full_A(:,:)
    double complex, intent(inout) :: full_X(:,:) ! on entry this contains full_B, on exit the solution

    integer, allocatable :: ipvt(:)
    integer :: lda, nRHSs
    external :: zgetrf, zgetrs ! LAPACK

    lda  =  size(full_A, 1)
    assert( size(full_A, 2) == n )
    nRHSs = size(full_X, 2)
    assert( size(full_X, 1) == lda ) ! must match the dims of A
    
    allocate(ipvt(lda))

    ! factorization
    call zgetrf(n, n, full_A, lda, ipvt, info) ! LU-factorize
!   if (info /= 0) ! ToDo: warn

    call zgetrs('n', n, nRHSs, full_A, lda, ipvt, full_X, lda, info) ! solve the system of linear equations
!   if (info /= 0) ! ToDo: warn

    deallocate(ipvt, stat=info)
  endfunction ! solveFull
  
  !----------------------------------------------------------------------------
  elemental subroutine destroy_solver(self)
    type(DirectSolver), intent(inout) :: self
    
    integer :: ist
    deallocate(self%full_A, self%full_X, stat=ist) ! ignore status
  endsubroutine ! destroy

endmodule ! DirectSolver_mod

