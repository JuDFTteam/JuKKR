module Solver_mod
  use Operator_mod
  use Precond_mod
  implicit none

  type Solver
    class(Operator), pointer :: op => null()
    class(Precond), pointer :: prec => null()
    contains
      procedure :: init => init_solver
      procedure :: init_precond => init_precond_solver
      procedure :: solve => solve_with_solver
  end type

  contains

  subroutine init_solver(self, op)
    class(Solver) :: self
    class(Operator), target :: op
    self%op => op
  end subroutine

  subroutine init_precond_solver(self, prec)
    class(Solver) :: self
    class(Precond), target :: prec
    self%prec => prec
  end subroutine

  subroutine solve_with_solver(self)
    class(Solver) :: self
    call self%op%apply
    if (associated(self%prec)) call self%prec%apply
  end subroutine
end module

