module TFQMRSolver_mod
  use Solver_mod
  use OperatorT_mod
  implicit none

  type TFQMRSolver
    private
    class(OperatorT), pointer :: op => null()
    class(OperatorT), pointer :: prec => null()
    contains
      procedure :: init => init_solver
      procedure :: init_precond => init_precond_solver
      procedure :: solve => solve_with_solver
  end type

  contains

  subroutine init_solver(self, op)
    class(TFQMRSolver) :: self
    class(OperatorT), target :: op
    self%op => op
  end subroutine

  subroutine init_precond_solver(self, prec)
    class(TFQMRSolver) :: self
    class(OperatorT), target :: prec
    self%prec => prec
  end subroutine

  subroutine solve_with_solver(self, mat_X, mat_B)
    class(TFQMRSolver) :: self
    double complex, intent(inout) :: mat_X(:,:)
    double complex, intent(in)    :: mat_B(:,:)
    !call self%op%apply
    !if (associated(self%prec)) call self%prec%apply
  end subroutine
end module

