module Solver_mod
  use OperatorT_mod
  implicit none

  type Solver
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
    class(Solver) :: self
    class(OperatorT), target :: op
    self%op => op
  end subroutine

  subroutine init_precond_solver(self, prec)
    class(Solver) :: self
    class(OperatorT), target :: prec
    self%prec => prec
  end subroutine

  subroutine solve_with_solver(self)
    class(Solver) :: self
    !call self%op%apply
    !if (associated(self%prec)) call self%prec%apply
  end subroutine
end module

