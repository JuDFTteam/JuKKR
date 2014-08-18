module Operator_mod
  implicit none

  type, abstract :: Operator
    contains
    procedure(apply), deferred :: apply
  end type

  interface
    !----------------------------------------------------------------------------
    !> Applies Operator on mat_X and returns result in mat_AX
    subroutine apply(self, mat_X, mat_AX)
      import Operator
      class(Operator) :: self
      double complex, intent(in)  :: mat_X
      double complex, intent(out) :: mat_AX
    end subroutine
  end interface

end module
