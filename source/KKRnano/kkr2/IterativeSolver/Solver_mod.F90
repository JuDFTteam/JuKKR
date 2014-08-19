module Solver_mod
  implicit none

  type, abstract :: Solver
    contains
      procedure (solve_interface), deferred :: solve
  end type

  interface
    !> Specify a right hand side mat_B - get result in mat_X
    subroutine solve_interface(self, mat_X, mat_B)
      import Solver
      class(Solver) :: self
      double complex, intent(inout) :: mat_X(:,:)
      double complex, intent(inout) :: mat_B(:,:)
    end subroutine
  end interface

end module

