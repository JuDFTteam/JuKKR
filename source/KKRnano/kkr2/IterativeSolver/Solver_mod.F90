!> A simple linear equation solver class.
!>
!> Extendfrom this class to implement a solver.
!>
!> @author Elias Rabel

module Solver_mod
  implicit none
  private
  public :: Solver

  type, abstract :: Solver
    contains
      procedure (solve_interface), deferred :: solve
  endtype

  interface
    !> Specify a right hand side mat_B - get result in mat_X
    subroutine solve_interface(self, mat_X, mat_B)
      import Solver
      class(Solver) :: self
      double complex, intent(inout) :: mat_X(:,:)
      double complex, intent(inout) :: mat_B(:,:)
    endsubroutine
  endinterface

endmodule

